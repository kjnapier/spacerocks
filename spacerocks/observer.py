from .spacerock import SpaceRock
from .units import Units

from astropy import units as u
import numpy as np

import pandas as pd

from skyfield.api import Topos, Loader
# Load in planets for ephemeride calculation.
load = Loader('./Skyfield-Data', expire=False, verbose=False)
planets = load('de440.bsp')
earth = planets['earth']
sun = planets['sun']


from skyfield.api import wgs84
from skyfield.data import iers

ts = load.timescale()

import pkg_resources

DATA_PATH = pkg_resources.resource_filename('spacerocks', 'data/observatories.csv')

observatories = pd.read_csv(DATA_PATH)


class Observer:

    def __init__(self, obscode):
        self.obscode = obscode

    def at(self, epoch):
        unique_times = np.unique(alltimes)
        t = ts.tdb(jd=unique_times)
        earth = planets['earth']

        #observer = Observer(epoch=self.epoch)

        if obscode == 'ssb':
            x_observer = np.zeros(len(alltimes)) * u.au
            y_observer = np.zeros(len(alltimes)) * u.au
            z_observer = np.zeros(len(alltimes)) * u.au

            vx_observer = np.zeros(len(alltimes)) * u.au / u.day
            vy_observer = np.zeros(len(alltimes)) * u.au / u.day
            vz_observer = np.zeros(len(alltimes)) * u.au / u.day

        elif obscode == 'sun':
            ss = sun.at(t)

            # xx, yy, zz = ss.position.au * u.au # earth ICRS position
            # vxx, vyy, vzz = ss.velocity.au_per_d * u.au / u.day # earth ICRS position

            xx, yy, zz = ss.ecliptic_xyz().au * u.au # earth ICRS position
            vxx, vyy, vzz = ss.ecliptic_velocity().au_per_d * u.au / u.day # earth ICRS position


            sxs = {t:x.value for t, x in zip(unique_times, xx)}
            sys = {t:y.value for t, y in zip(unique_times, yy)}
            szs = {t:z.value for t, z in zip(unique_times, zz)}
            svxs = {t:vx.value for t, vx in zip(unique_times, vxx)}
            svys = {t:vy.value for t, vy in zip(unique_times, vyy)}
            svzs = {t:vz.value for t, vz in zip(unique_times, vzz)}

            x_observer = np.array([sxs[t] for t in alltimes]) * u.au
            y_observer = np.array([sys[t] for t in alltimes]) * u.au
            z_observer = np.array([szs[t] for t in alltimes]) * u.au

            vx_observer = np.array([svxs[t] for t in alltimes]) * u.au / u.day
            vy_observer = np.array([svys[t] for t in alltimes]) * u.au / u.day
            vz_observer = np.array([svzs[t] for t in alltimes]) * u.au / u.day

        # Only used for the topocentric calculation.
        else:
            if (obscode != 500) and (obscode != '500'):
                earth += wgs84.latlon(obslat, obslon, obselev)

            ee = earth.at(t)

            # xx, yy, zz = ee.position.au * u.au # earth ICRS position
            # vxx, vyy, vzz = ee.velocity.au_per_d * u.au / u.day # earth ICRS position

            xx, yy, zz = ee.ecliptic_xyz().au * u.au # earth ICRS position
            vxx, vyy, vzz = ee.ecliptic_velocity().au_per_d * u.au / u.day # earth ICRS position

            exs = {t:x.value for t, x in zip(unique_times, xx)}
            eys = {t:y.value for t, y in zip(unique_times, yy)}
            ezs = {t:z.value for t, z in zip(unique_times, zz)}
            evxs = {t:vx.value for t, vx in zip(unique_times, vxx)}
            evys = {t:vy.value for t, vy in zip(unique_times, vyy)}
            evzs = {t:vz.value for t, vz in zip(unique_times, vzz)}

            x_observer = np.array([exs[t] for t in alltimes]) * u.au
            y_observer = np.array([eys[t] for t in alltimes]) * u.au
            z_observer = np.array([ezs[t] for t in alltimes]) * u.au

            vx_observer = np.array([evxs[t] for t in alltimes]) * u.au / u.day
            vy_observer = np.array([evys[t] for t in alltimes]) * u.au / u.day
            vz_observer = np.array([evzs[t] for t in alltimes]) * u.au / u.day


def observer(obscodes, epochs):

    x = []
    y = []
    z = []
    vx = []
    vy = []
    vz = []

    for obscode, epoch in zip(obscodes, epochs):

        print(obscode, epoch)

        earth = planets['earth']

        if obscode != 500:
            obscode = str(obscode).zfill(3)
            obs = observatories[observatories.obscode == obscode]
            obslat = obs.lat.values
            obslon = obs.lon.values
            obselev = obs.elevation.values

            earth += wgs84.latlon(obslat, obslon, obselev)
                    

        print(epoch.tdb.jd)
        
        t = ts.tdb(jd=epoch.tdb.jd)
        
        print(t)
        ee = earth.at(t)
        print('.')
        xx, yy, zz = ee.ecliptic_xyz().au 
        vxx, vyy, vzz = ee.ecliptic_velocity().au_per_d 

        # ephemeris_time = spice.str2et(['JD{} TDB'.format(epoch.tdb.jd)])
        # state, _ = spice.spkezr(obscode, ephemeris_time, 'J2000', 'none', '0')
        # state = np.vstack(state).T
        #xx, yy, zz, vxx, vyy, vzz = state
        x.append(xx)
        y.append(yy)
        z.append(zz)
        vx.append(vxx)
        vy.append(vyy)
        vz.append(vzz)

    units = Units()
    units.distance = u.km
    units.speed = u.km / u.s
    units.timescale = 'tdb'

    rocks = SpaceRock(x=x,
                      y=y,
                      z=z,
                      vx=vx,
                      vy=vy,
                      vz=vz,
                      units=units,
                      epoch=epochs.tdb.jd,
                      frame='barycentric')

    return rocks
