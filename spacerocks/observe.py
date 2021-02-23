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
from .transformations import Transformations
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

class Observe(Transformations, Convenience):

    def __init__(self, rocks, obscode=None, NSIDE=None):

        if obscode is not None:
            self.__class__.obscode = str(obscode).zfill(3)
            obs = observatories[observatories.obscode == self.__class__.obscode]
            self.__class__.obslat = obs.lat.values
            self.__class__.obslon = obs.lon.values
            self.__class__.obselev = obs.elevation.values
        else:
            self.__class__.obscode = 500

        if NSIDE is not None:
            if isscalar(NSIDE):
                NSIDE = array([NSIDE])
            self.__class__.NSIDE = NSIDE
        else:
            self.__class__.NSIDE = None

        self.xyz_to_equa(rocks)

        if rocks.frame == 'barycentric':
            self.xyz_to_equa(rocks)

        elif rocks.frame == 'heliocentric':
            rocks.to_bary()
            self.xyz_to_equa(rocks)
            rocks.to_helio()


        if self.__class__.NSIDE is not None:
            for value in Observe.NSIDE:
                setattr(self, 'HPIX_{}'.format(value), self.radec_to_hpix(value))

        self.name = rocks.name
        self.epoch = rocks.epoch

        try:
            self.mag = self.estimate_mag(rocks)
        except:
            pass

    def radec_to_hpix(self, NSIDE):
        '''
        Convert (ra, dec) into healpix.
        '''
        return hp.pixelfunc.ang2pix(NSIDE, np.pi/2 - self.dec.radian, self.ra.radian, nest=True)

    @property
    def ecliptic_longitude(self):
        return arctan2((cos(epsilon)*cos(self.dec.rad)*sin(self.ra.rad) + sin(epsilon)*sin(self.dec.rad)), cos(self.dec.rad) * cos(self.ra.rad))

    @property
    def ecliptic_latitude(self):
        return arcsin(cos(epsilon)*sin(self.dec.rad) - sin(epsilon)*cos(self.dec.rad)*sin(self.ra.rad))
