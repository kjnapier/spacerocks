import os
import warnings

import healpy as hp
import pandas as pd
import numpy as np
from skyfield.api import Topos, Loader
from astropy import units as u
from astropy.coordinates import Angle, SkyCoord

from numpy import sin, cos, arctan2, arcsin, array, isscalar, sqrt, pi

from .constants import epsilon, mu_bary, c
#from .skyfuncs import Transformations
from .convenience import Convenience
from .vector import Vector

# Read in the observatory codes file and rehash as a dataframe.
observatories = pd.read_csv(os.path.join(os.path.dirname(__file__),
                            'data',
                            'observatories.csv'))

# Load in planets for ephemeride calculation.
load = Loader('./Skyfield-Data', expire=False, verbose=False)
#ts = load.timescale()
planets = load('de423.bsp')

from skyfield.api import wgs84

from skyfield.data import iers

url = load.build_url('finals2000A.all')
with load.open(url) as f:
    finals_data = iers.parse_x_y_dut1_from_finals_all(f)

ts = load.timescale()
iers.install_polar_motion_table(ts, finals_data)

class Observe(Convenience):

    def __init__(self, rocks, obscode=None, NSIDE=None):

        self.mu = mu_bary
        if obscode is not None:
            self.__class__.obscode = str(obscode).zfill(3)
            obs = observatories[observatories.obscode == self.__class__.obscode]
            self.__class__.obslat = obs.lat.values
            self.__class__.obslon = obs.lon.values
            self.__class__.obselev = obs.elevation.values
        else:
            self.__class__.obscode = 500


        #self.tel_position, self.tel_velocity = self.xyz_to_tel(rocks)
        self.xT, self.yT, self.zT, self.vxT, self.vyT, self.vzT = self.xyz_to_tel(rocks)


        self.name = rocks.name
        self.epoch = rocks.epoch

        self.dec = Angle(np.arcsin(self.zT / sqrt(self.xT**2 + self.yT**2 + self.zT**2)), u.rad)
        self.ra = Angle(np.arctan2(self.yT, self.xT), u.rad).wrap_at(2 * np.pi * u.rad)
        self.dec_rate = (-self.zT * (self.xT * self.vxT + self.yT * self.vyT) + ((self.xT**2 + self.yT**2) * self.vzT)) \
                / (sqrt(self.xT**2 + self.yT**2) * (self.xT**2 + self.yT**2 + self.zT**2)) * u.rad
        self.ra_rate = - cos(self.dec) * (self.yT * self.vxT - self.xT * self.vyT) / (self.xT**2 + self.yT**2) * u.rad
        self.delta = sqrt(self.xT**2 + self.yT**2 + self.zT**2)
        earth_dis = 1 * u.au
        self.phase_angle = Angle(np.arccos(-(earth_dis**2 - rocks.r**2 - self.delta**2)/(2 * rocks.r * self.delta)), u.rad)
        self.elong = Angle(np.arccos(-(rocks.r**2 - self.delta**2 - earth_dis**2)/(2 * self.delta * earth_dis)), u.rad)

    @property
    def mag(self):
        return self._mag

    @mag.setter
    def mag(self, rocks):
        self._mag = self.__estimate_mag(rocks)


    def xyz_to_tel(self, rocks):
        '''
        Transform from barycentric ecliptic Cartesian coordinates to
        telescope-centric coordinates.

        Routine corrects iteratively for light travel time.
        '''

        rocks.to_bary()

        #t = ts.tdb(jd=rocks.epoch.tdb.jd)
        alltimes = rocks.epoch.tdb.jd
        unique_times = np.unique(alltimes)
        t = ts.tdb(jd=unique_times)
        earth = planets['earth']
        earth += wgs84.latlon(self.__class__.obslat, self.__class__.obslon, elevation_m=self.__class__.obselev)
        ee = earth.at(t)
        #x_earth, y_earth, z_earth = ee.ecliptic_xyz().au * u.au # earth ICRS position
        #vx_earth, vy_earth, vz_earth = ee.ecliptic_velocity().au_per_d * u.au / u.day # earth ICRS position

        xx, yy, zz = ee.position.au * u.au # earth ICRS position
        vxx, vyy, vzz = ee.velocity.au_per_d * u.au / u.day # earth ICRS position

        exs = {t:x for t, x in zip(unique_times, xx)}
        eys = {t:y for t, y in zip(unique_times, yy)}
        ezs = {t:z for t, z in zip(unique_times, zz)}
        evxs = {t:vx for t, vx in zip(unique_times, vxx)}
        evys = {t:vy for t, vy in zip(unique_times, vyy)}
        evzs = {t:vz for t, vz in zip(unique_times, vzz)}

        x_earth = np.array([exs[t] for t in alltimes])
        y_earth = np.array([eys[t] for t in alltimes])
        z_earth = np.array([ezs[t] for t in alltimes])

        vx_earth = np.array([evxs[t] for t in alltimes])
        vy_earth = np.array([evys[t] for t in alltimes])
        vz_earth = np.array([evzs[t] for t in alltimes])

        x0, y0, z0 = rocks.x, rocks.y, rocks.z
        vx0, vy0, vz0 = rocks.vx, rocks.vy, rocks.vz

        for idx in range(10):

            dx = x0 - x_earth
            dy = y0 * np.cos(epsilon) - z0 * np.sin(epsilon) - y_earth
            dz = y0 * np.sin(epsilon) + z0 * np.cos(epsilon) - z_earth
            dvx = vx0 - vx_earth
            dvy = vy0 * np.cos(epsilon) - vz0 * np.sin(epsilon) - vy_earth
            dvz = vy0 * np.sin(epsilon) + vz0 * np.cos(epsilon) - vz_earth

            if idx < 9:

                delta = sqrt(dx**2 + dy**2 + dz**2)
                ltt = delta / c
                M = rocks.M - (ltt * rocks.n)
                x0, y0, z0, vx0, vy0, vz0 = self.kep_to_xyz_temp(rocks.a,
                                                                 rocks.e,
                                                                 rocks.inc,
                                                                 rocks.arg,
                                                                 rocks.node,
                                                                 M)

        #return Vector(xT, yT, zT), Vector(vxT, vyT, vzT)
        #return xT, yT, zT, vxT, vyT, vzT
        return dx, dy, dz, dvx, dvy, dvz

    def kep_to_xyz_temp(self, a, e, inc, arg, node, M):
        '''
        Just compute the xyz position of an object. Used for iterative equatorial
        calculation.
        '''
        # compute eccentric anomaly E
        e = e
        M = array(M.rad)
        M[M > pi] -= 2 * pi
        alpha = (3 * pi**2 + 1.6 * (pi**2 - pi * abs(M))/(1 + e))/(pi**2 - 6)
        d = 3 * (1 - e) + alpha * e
        q = 2 * alpha * d * (1 - e) - M**2
        r = 3 * alpha * d * (d - 1 + e) * M + M**3
        w = (abs(r) + sqrt(q**3 + r**2))**(2/3)
        E1 = (2 * r * w / (w**2 + w*q + q**2) + M)/d
        f2 = e * sin(E1)
        f3 = e * cos(E1)
        f0 = E1 - f2 - M
        f1 = 1 - f3
        d3 = -f0 / (f1 - f0 * f2 / (2 * f1))
        d4 = -f0 / (f1 + f2 * d3 / 2 + d3**2 * f3/6)
        d5 = -f0 / (f1 + d4*f2/2 + d4**2*f3/6 - d4**3*f2/24)
        E = E1 + d5
        E[E < 0] += 2 * pi
        #E = E % (2 * pi)

        # compute true anomaly ν
        true_anomaly = 2 * np.arctan2((1 + e)**0.5*np.sin(E/2), (1 - e)**0.5*np.cos(E/2))

        # compute the distance to the central body r
        r = a * (1 - e * np.cos(E))

        # obtain the position vector o
        o = Vector(r * cos(true_anomaly), r * sin(true_anomaly), np.zeros_like(true_anomaly))
        ov = Vector((mu_bary * a)**0.5 / r * (-np.sin(E))/ u.rad, (mu_bary * a)**0.5 / r * ((1 - e**2)**0.5 * np.cos(E))/ u.rad, np.zeros(len(true_anomaly))/ u.rad)

        # Rotate o to the inertial frame
        position = o.euler_rotation(arg, inc, node) #* u.au
        velocity = ov.euler_rotation(arg, inc, node) #* u.au / u.day

        return position.x, position.y, position.z, velocity.x, velocity.y, velocity.z
