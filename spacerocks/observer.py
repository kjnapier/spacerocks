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
