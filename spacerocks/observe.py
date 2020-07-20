import os
import warnings

import healpy as hp
import numpy as np
import pandas as pd
from skyfield.api import Topos, Loader
from astropy import units as u
from astropy.coordinates import Angle, SkyCoord

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
            Observe.obscode = str(obscode).zfill(3)
            obs = observatories[observatories.obscode == Observe.obscode]
            Observe.obslat = obs.lat.values
            Observe.obslon = obs.lon.values
            Observe.obselev = obs.elevation.values
        else:
            Observe.obscode = 500

        if NSIDE is not None:
            if np.isscalar(NSIDE):
                NSIDE = np.array([NSIDE])
            Observe.NSIDE = NSIDE
        else:
            Observe.NSIDE = None

        self.xyz_to_equa(rocks)

        if rocks.frame == 'barycentric':
            self.xyz_to_equa(rocks)

        elif rocks.frame == 'heliocentric':
            rocks.to_bary()
            self.xyz_to_equa(rocks)
            rocks.to_helio()

        try:
            self.mag = self.estimate_mag(rocks)
        except:
            pass

        if Observe.NSIDE is not None:
            for value in Observe.NSIDE:
                setattr(self, 'HPIX_{}'.format(value), self.radec_to_hpix(value))

        self.name = rocks.name
        self.epoch = rocks.epoch

        try:
            self.mag = self.estimate_mag(rocks)
        except:
            pass
