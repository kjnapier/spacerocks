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
ts = load.timescale()
planets = load('de423.bsp')

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
        self.ra_rate = -(self.yT * self.vxT - self.xT * self.vyT) / (self.xT**2 + self.yT**2) * u.rad
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

        t = ts.tt(jd=rocks.epoch.tt.jd)
        earth = planets['earth']

        # Only used for the topocentric calculation.
        #if self.__class__.obscode != 500:
        #    earth += Topos(latitude_degrees=self.__class__.obslat,
        #                   longitude_degrees=self.__class__.obslon,
        #                   elevation_m=self.__class__.obselev) # topocentric calculation

        ee = earth.at(t)
        #x_earth, y_earth, z_earth = ee.ecliptic_xyz().au * u.au # earth ICRS position
        #vx_earth, vy_earth, vz_earth = ee.ecliptic_velocity().au_per_d * u.au / u.day # earth ICRS position

        x_earth, y_earth, z_earth = ee.position.au * u.au # earth ICRS position
        vx_earth, vy_earth, vz_earth = ee.velocity.au_per_d * u.au / u.day # earth ICRS position

        x0, y0, z0 = rocks.x, rocks.y, rocks.z
        vx0, vy0, vz0 = rocks.vx, rocks.vy, rocks.vz

        for idx in range(10):
            # transfer ecliptic to ICRS and shift to Geocentric (topocentric)
            x = x0 - x_earth
            y = y0 * np.cos(epsilon) - z0 * np.sin(epsilon) - y_earth
            z = y0 * np.sin(epsilon) + z0 * np.cos(epsilon) - z_earth
            vx = vx0 - vx_earth
            vy = vy0 * np.cos(epsilon) - vz0 * np.sin(epsilon) - vy_earth
            vz = vy0 * np.sin(epsilon) + vz0 * np.cos(epsilon) - vz_earth
            delta = sqrt(x**2 + y**2 + z**2)
            ltt = delta / c
            M = rocks.M - ltt * (mu_bary / rocks.a**3)**0.5
            if idx < 9:
                x0, y0, z0, vx0, vy0, vz0 = self.kep_to_xyz_temp(rocks.a, rocks.e, rocks.inc,
                                                                 rocks.arg, rocks.node, M)

        #for idx in range(5):

        #    # transfer ecliptic to ICRS and shift to Geocentric (topocentric)
        #    #xT = x0 - x_earth
        #    #yT = y0 * cos(epsilon) - z0 * sin(epsilon) - y_earth
        #    #zT = y0 * sin(epsilon) + z0 * cos(epsilon) - z_earth
        #    #vxT = vx0 - vx_earth
        #    #vyT = vy0 * cos(epsilon) - vz0 * sin(epsilon) - vy_earth
        #    #vzT = vy0 * sin(epsilon) + vz0 * cos(epsilon) - vz_earth

        #    xT = x0 - x_earth
        #    yT = y0 - y_earth * cos(-epsilon) - z_earth * sin(epsilon)
        #    zT = z0 - y_earth * sin(-epsilon) - z_earth * cos(epsilon)
        #    vxT = vx0 - vx_earth
        #    vyT = vy0 - vy_earth * cos(-epsilon) - vz_earth * sin(epsilon)
        #    vzT = vz0 - vy_earth * sin(-epsilon) - vz_earth * cos(epsilon)


        #    if idx < 4:

        #        delta = sqrt(xT**2 + yT**2 + zT**2)
        #        ltt = delta / c
        #        M = rocks.M - (ltt * rocks.n)
        #        x0, y0, z0, vx0, vy0, vz0 = self.kep_to_xyz_temp(rocks.a,
        #                                                         rocks.e,
        #                                                         rocks.inc,
        #                                                         rocks.arg,
        #                                                         rocks.node,
        #                                                         M)

        #return Vector(xT, yT, zT), Vector(vxT, vyT, vzT)
        #return xT, yT, zT, vxT, vyT, vzT
        return x, y, z, vx, vy, vz

    def kep_to_xyz_temp(self, a, e, inc, arg, node, M):
        '''
        Just compute the xyz position of an object. Used for iterative equatorial
        calculation.
        '''
        # compute eccentric anomaly E
        e = e.value
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
