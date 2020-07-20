################################################################################
# SpaceRocks, version 0.7.5
#
# 0.7.2:
#     - changed the obsdate keyword to epoch
#
# 0.7.4:
#     - Major restructuring. Addec Transformations, Convenience, Observe classes
#
# Author: Kevin Napier kjnapier@umich.edu
################################################################################

import healpy as hp
import sys
import os
import random
import copy
import warnings

import numpy as np
import pandas as pd

from astropy import units as u
from astropy.table import Table
from astropy.coordinates import Angle
from astropy.time import Time
from astropy.coordinates import SkyCoord

# Read in the observatory codes file and rehash as a dataframe.
observatories = pd.read_csv(os.path.join(os.path.dirname(__file__),
                            'data',
                            'observatories.csv'))

from skyfield.api import Topos, Loader
# Load in planets for ephemeride calculation.
load = Loader('./Skyfield-Data', expire=False, verbose=False)
ts = load.timescale()
planets = load('de423.bsp')
earth = planets['earth']

from .linalg3d import *
from .constants import *
from .transformations import Transformations
from .convenience import Convenience
#from .jacobians import *


class SpaceRock(Transformations, Convenience):

    def __init__(self, input_coordinates='keplerian', input_frame='barycentric',
                 input_angles='degrees', input_time_format='jd', calc_abg=False,
                 input_time_scale='utc', abg_obscode=500, *args, **kwargs):

        if abg_obscode != 500:
            SpaceRock.obscode = str(abg_obscode).zfill(3)
            obs = observatories[observatories.obscode == SpaceRock.obscode]
            SpaceRock.obslat = obs.lat.values
            SpaceRock.obslon = obs.lon.values
            SpaceRock.obselev = obs.elevation.values
        else:
            SpaceRock.obscode = 500

        # Case-insensitive keyword arguments.
        kwargs = {key.lower(): data for key, data in kwargs.items()}
        # keywords = ['a', 'e', 'inc', 'node', 'arg', 'm',
        #             'x', 'y', 'z', 'vx', 'vy', 'vz',
        #             'epoch', 't_peri', 'h', 'name', 'g']

        # if not all(key in keywords for key in [*kwargs]):
        #     raise ValueError('Keywords are limited to a, e, inc, node,\
        #                       arg, m, x, y, z, vx, vy, vz, epoch, t_peri,\
        #                       h, name, g')

        input_coordinates = input_coordinates.lower()
        input_frame = input_frame.lower()
        input_angles = input_angles.lower()
        calc_abg = calc_abg

        # scalar input -> arrays
        for idx, key in enumerate([*kwargs]):
            if np.isscalar(kwargs.get(key)):
                kwargs[key] = np.array([kwargs.get(key)])

        if kwargs.get('delta_h') is not None:
            self.delta_H = kwargs.get('delta_h')
        if kwargs.get('rotation_period') is not None:
            self.rotation_period = kwargs.get('rotation_period') * u.day
        if kwargs.get('phi0') is not None:
            self.phi_0 = kwargs.get('phi0')


        self.t0 = kwargs.get('epoch') * u.day

        if kwargs.get('name') is not None:
            self.name = kwargs.get('name')
        else:
            # produces random, non-repeting integers between 0 and 1e10 - 1
            self.name = np.array(['{:010}'.format(value) for value in random.sample(range(int(1e10)), len(self.t0))])

        self.epoch = Time(kwargs.get('epoch'), format=input_time_format, scale=input_time_scale)

        SpaceRock.frame = input_frame


        if SpaceRock.frame == 'barycentric':
            mu = mu_bary
        elif SpaceRock.frame == 'heliocentric':
            mu = mu_helio

        if input_angles == 'degrees':
            angle_unit = u.degree
        elif input_angles == 'radians':
            angle_unit = u.radian
        else:
            raise ValueError('The input_angles argument must be a string \
                              that reads either degrees or radians.')

        if input_coordinates == 'keplerian':

            self.a = kwargs.get('a') * u.au
            self.e = kwargs.get('e') * u.dimensionless_unscaled
            self.inc = Angle(kwargs.get('inc'), angle_unit).to(u.rad)
            self.node = Angle(kwargs.get('node'), angle_unit).to(u.rad)
            self.arg = Angle(kwargs.get('arg'), angle_unit).to(u.rad)

            if (kwargs.get('t_peri') is None) and (kwargs.get('m') is not None):
                self.M = Angle(kwargs.get('m'), angle_unit).rad * u.rad
                self.t_peri = self.calc_t_peri()

            elif (kwargs.get('m') is None) and (kwargs.get('t_peri') is not None):
                self.t_peri = Time(kwargs.get('t_peri'), format=input_time_format, scale=input_time_scale)
                self.M = np.sqrt(mu / self.a**3) * (self.epoch.jd * u.day - self.t_peri.jd * u.day)

            self.kep_to_xyz(mu)

            if calc_abg == True:
                if self.__class__.frame == 'heliocentric':
                    self.to_bary()
                    self.xyz_to_abg()
                    self.to_helio()
                else:
                    self.xyz_to_abg()


        elif input_coordinates == 'cartesian':

            self.x = kwargs.get('x') * u.au
            self.y = kwargs.get('y') * u.au
            self.z = kwargs.get('z') * u.au
            self.vx = kwargs.get('vx') * (u.au / u.day)
            self.vy = kwargs.get('vy') * (u.au / u.day)
            self.vz = kwargs.get('vz') * (u.au / u.day)

            self.xyz_to_kep(mu)

            if calc_abg == True:
                if self.__class__.frame == 'heliocentric':
                    self.to_bary()
                    self.xyz_to_abg()
                    self.to_helio()
                else:
                    self.xyz_to_abg()


        elif input_coordinates == 'abg':

            self.alpha = kwargs.get('alpha')
            self.beta = kwargs.get('beta')
            self.gamma = kwargs.get('gamma') / u.au
            self.alpha_dot = kwargs.get('alpha_dot') / u.day
            self.beta_dot = kwargs.get('beta_dot') / u.day
            self.gamma_dot = kwargs.get('gamma_dot') / u.day

            self.abg_to_xyz()
            self.xyz_to_kep(mu)


        self.varpi = Angle((self.arg + self.node).wrap_at(2 * np.pi * u.rad), u.rad)

        if kwargs.get('h') is not None:
            self.H = kwargs.get('h')
            self.G = kwargs.get('g')
            if self.G is None:
                self.G = np.repeat(0.15, len(self))

        # self.epoch = Time(self.epoch, format='jd', scale='utc')

    def plot_orbits(self):
        '''
        Plot the orbits of all your rocks.
        '''
        x = np.zeros([500, len(self.a)])
        y = np.zeros([500, len(self.a)])
        z = np.zeros([500, len(self.a)])
        for idx, M in enumerate(np.linspace(0, 2*np.pi, 500)):
            xx, yy, zz, _, _, _ = self.kep_to_xyz_temp(self.a, self.e, self.inc.rad,
                                                       self.arg.rad, self.node.rad,
                                                       np.repeat(M, len(self.a)) * u.rad)
            x[idx] = xx
            y[idx] = yy
            z[idx] = zz

        fig = plt.figure(figsize=(12, 12))
        ax = fig.add_subplot(111, projection='3d')
        for idx in range(len(self.a)):
            ax.plot(x.T[idx], y.T[idx], z.T[idx], color=np.random.choice(['#00274C', '#FFCB05']))

        ax.axis('off')
        ax.view_init(90, 0)

        return fig, ax
