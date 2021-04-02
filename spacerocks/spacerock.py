################################################################################
# SpaceRocks, version 0.7.5
#
# 0.7.2:
#     - changed the obsdate keyword to epoch
#
# 0.7.4:
#     - Major restructuring. Transformations, Convenience, Observe classes
#
# Author: Kevin Napier kjnapier@umich.edu
################################################################################

import sys
import os
import random
import copy

from math import pi

import numpy as np
import pandas as pd

from astropy import units as u
from astropy.coordinates import Angle, Distance
from astropy.time import Time

import dateutil
import matplotlib.pyplot as plt

from numpy import sin, cos, arctan2, sqrt

## Read in the observatory codes file and rehash as a dataframe.
#observatories = pd.read_csv(os.path.join(os.path.dirname(__file__),
#                            'data',
#                            'observatories.csv'))
#
#from skyfield.api import Topos, Loader
## Load in planets for ephemeride calculation.
#load = Loader('./Skyfield-Data', expire=False, verbose=False)
#ts = load.timescale()
#planets = load('de423.bsp')
#earth = planets['earth']

from .linalg3d import *
from .constants import *
from .orbitfuncs import OrbitFuncs
from .convenience import Convenience
from .observe import Observe
from .units import Units
#from .jacobians import *


class SpaceRock(OrbitFuncs, Convenience):

    def __init__(self, input_frame='barycentric', units=Units(), *args, **kwargs):

        # Case-insensitive keyword arguments.
        coords = self.detect_coords(kwargs)
        input_frame = input_frame.lower()

        # input -> arrays
        for idx, key in enumerate([*kwargs]):
            if not hasattr(kwargs.get(key), '__len__'):
                kwargs[key] = np.array([kwargs.get(key)])
            else:
                kwargs[key] = np.array(kwargs.get(key))

        SpaceRock.frame = input_frame
        if SpaceRock.frame == 'barycentric':
            SpaceRock.mu = mu_bary
        elif SpaceRock.frame == 'heliocentric':
            SpaceRock.mu = mu_helio

        if units.timeformat is None:
            self.epoch = self.detect_timescale(kwargs.get('epoch'), units.timescale)
        else:
            self.epoch = Time(kwargs.get('epoch'), format=units.timeformat, scale=units.timescale)

        self.t0 = Time(self.epoch.jd, format='jd', scale=units.timescale)

        if kwargs.get('name') is not None:
            self.name = kwargs.get('name').astype(str)
        else:
            # produces random, non-repeting integers between 0 and 1e10 - 1
            self.name = np.array(['{:010}'.format(value) for value in random.sample(range(int(1e10)), len(self.epoch))])



        if coords == 'kep':

            self.a = Distance(kwargs.get('a'), units.distance).to(u.au)
            self.e = kwargs.get('e') #* u.dimensionless_unscaled
            self.inc = Angle(kwargs.get('inc'), units.angle).to(u.rad)

            if kwargs.get('node') is not None:
                self.node = Angle(kwargs.get('node'), units.angle).to(u.rad)

            if kwargs.get('arg') is not None:
                self.arg = Angle(kwargs.get('arg'), units.angle).to(u.rad)

            if kwargs.get('varpi') is not None:
                self.varpi = Angle(kwargs.get('varpi'), units.angle).to(u.rad)

            if kwargs.get('t_peri') is not None:
                if units.timeformat is None:
                    self.t_peri = self.detect_timescale(kwargs.get('t_peri'), units.timescale)
                else:
                    self.t_peri = Time(kwargs.get('t_peri'), format=units.timeformat, scale=units.timescale)

            if kwargs.get('M') is not None:
                self.M = Angle(kwargs.get('M'), units.angle).to(u.rad)
                self.E = self._calc_E_from_M(self.e, self.M.rad)

            if kwargs.get('E') is not None:
                self.E = Angle(kwargs.get('E'), units.angle).to(u.rad)
                self.M = self.E * self.e * sin(self.E)

            if kwargs.get('true_anomaly') is not None:
                self.true_anomaly = Angle(kwargs.get('true_anomaly'), units.angle)
                self.E = Angle(2 * arctan2(sqrt(1-self.e) * sin(self.true_anomaly/2), sqrt(1+self.e) * cos(self.true_anomaly/2)), u.rad)

            if kwargs.get('true_longitude') is not None:
                self.true_longitude = Angle(kwargs.get('true_longitude'), units.angle).to(u.rad)
                self.true_anomaly = self.true_longitude - self.varpi

            if kwargs.get('mean_longitude') is not None:
                self.mean_longitude = Angle(kwargs.get('mean_longitude'), units.angle)
                self.M = Angle((self.mean_longitude.rad - self.varpi.rad) % (2 * pi), u.rad)


            self.kep_to_xyz()


        elif coords == 'xyz':

            self.x = Distance(kwargs.get('x'), units.distance, allow_negative=True).to(u.au)
            self.y = Distance(kwargs.get('y'), units.distance, allow_negative=True).to(u.au)
            self.z = Distance(kwargs.get('z'), units.distance, allow_negative=True).to(u.au)
            self.vx = (kwargs.get('vx') * units.speed).to(u.au / u.day)
            self.vy = (kwargs.get('vy') * units.speed).to(u.au / u.day)
            self.vz = (kwargs.get('vz') * units.speed).to(u.au / u.day)

            self.xyz_to_kep()
