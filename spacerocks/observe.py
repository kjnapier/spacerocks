import os
import warnings

import healpy as hp
import pandas as pd
import numpy as np
from skyfield.api import Topos, Loader
from astropy import units as u
from astropy.coordinates import Angle, SkyCoord

from numpy import sin, cos, arctan2, arcsin, array, isscalar

from .constants import epsilon, mu_bary, c
#from .skyfuncs import Transformations
from .linalg3d import norm
from .convenience import Convenience

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


        self.xT, self.yT, self.zT, self.vxT, self.vyT, self.vzT = self.xyz_to_tel(rocks)


        self.name = rocks.name
        self.epoch = rocks.epoch


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

        if rocks.frame == 'heliocentric':
            rocks.to_bary()

        t = ts.tt(jd=rocks.epoch.tt.jd)
        earth = planets['earth']

        # Only used for the topocentric calculation.
        if self.__class__.obscode != 500:
            earth += Topos(latitude_degrees=self.__class__.obslat,
                           longitude_degrees=self.__class__.obslon,
                           elevation_m=self.__class__.obselev) # topocentric calculation

        ee = earth.at(t)
        x_earth, y_earth, z_earth = ee.position.au * u.au # earth ICRS position
        vx_earth, vy_earth, vz_earth = ee.velocity.au_per_d * u.au / u.day # earth ICRS position

        for idx in range(5):

            # transfer ecliptic to ICRS and shift to Geocentric (topocentric)
            xT = rocks.x - x_earth
            yT = rocks.y * cos(epsilon) - rocks.z * sin(epsilon) - y_earth
            zT = rocks.y * sin(epsilon) + rocks.z * cos(epsilon) - z_earth
            vxT = rocks.vx - vx_earth
            vyT = rocks.vy * cos(epsilon) - rocks.vz * sin(epsilon) - vy_earth
            vzT = rocks.vy * sin(epsilon) + rocks.vz * cos(epsilon) - vz_earth

            if idx < 4:

                delta = norm([xT, yT, zT])
                ltt = delta / c
                M = rocks.M - ltt * rocks.n
                x0, y0, z0, vx0, vy0, vz0 = self.kep_to_xyz_temp(rocks.a,
                                                                 rocks.e,
                                                                 rocks.inc,
                                                                 rocks.arg,
                                                                 rocks.node,
                                                                 M)

        return xT, yT, zT, vxT, vyT, vzT
