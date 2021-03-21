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

import dateutil
import matplotlib.pyplot as plt

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
from .observe import Observe
from .units import Units
#from .jacobians import *


class SpaceRock(Transformations, Convenience):

    def __init__(self, input_frame='barycentric', units=Units(), *args, **kwargs):

        # Case-insensitive keyword arguments.
        kwargs = {key.lower(): data for key, data in kwargs.items()}
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
            mu = mu_bary
        elif SpaceRock.frame == 'heliocentric':
            mu = mu_helio

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

        if kwargs.get('delta_h') is not None:
            self.delta_H = kwargs.get('delta_h')
        if kwargs.get('rotation_period') is not None:
            self.rotation_period = kwargs.get('rotation_period') * units.rotation_period
        if kwargs.get('phi0') is not None:
            self.phi_0 = Angle(kwargs.get('phi0'), units.angle)


        if coords == 'kep':

            self.a = kwargs.get('a') * units.distance
            self.e = kwargs.get('e') * u.dimensionless_unscaled
            self.inc = Angle(kwargs.get('inc'), units.angle).to(u.rad)
            self.node = Angle(kwargs.get('node'), units.angle).to(u.rad)
            self.arg = Angle(kwargs.get('arg'), units.angle).to(u.rad)

            if (kwargs.get('t_peri') is None) and (kwargs.get('m') is not None):
                self.M = Angle(kwargs.get('m'), units.angle).to(u.rad)
                self.t_peri = self.calc_t_peri()

            elif (kwargs.get('m') is None) and (kwargs.get('t_peri') is not None):
                if units.timeformat is None:
                    self.t_peri = self.detect_timescale(kwargs.get('t_peri'), units.timescale)
                else:
                    self.t_peri = Time(kwargs.get('t_peri'), format=units.timeformat, scale=units.timescale)
                self.M = np.sqrt(mu / self.a**3) * (self.epoch.jd * u.day - self.t_peri.jd * u.day)

            self.kep_to_xyz(mu)


        elif coords == 'xyz':

            self.x = kwargs.get('x') * units.distance
            self.y = kwargs.get('y') * units.distance
            self.z = kwargs.get('z') * units.distance
            self.vx = kwargs.get('vx') * units.speed
            self.vy = kwargs.get('vy') * units.speed
            self.vz = kwargs.get('vz') * units.speed

            self.xyz_to_kep(mu)


        self.varpi = Angle((self.arg + self.node).wrap_at(2 * np.pi * u.rad), u.rad)
        self.true_longitude = Angle((self.true_anomaly + self.varpi).wrap_at(2 * np.pi * u.rad), u.rad)
        self.mean_longitude = Angle((self.M + self.varpi).wrap_at(2 * np.pi * u.rad), u.rad)

        if kwargs.get('h') is not None:
            self.H = kwargs.get('h')
            if kwargs.get('g') is not None:
                self.G = kwargs.get('g')
            else:
                self.G = np.repeat(0.15, len(self))

    @property
    def q(self):
        return self.a * (1 - self.e)


    def ecl_to_tel(self, x_ec, y_ec, z_ec, l0, b0):
        '''
        This is right
        '''
        sl = np.sin(l0)
        cl = np.cos(l0)
        sb = np.sin(b0)
        cb = np.cos(b0)

        xT = -sl*x_ec + cl*y_ec
        yT = -cl*sb*x_ec -sl*sb*y_ec + cb*z_ec
        zT = cl*cb*x_ec + sl*cb*y_ec + sb*z_ec

        return xT, yT, zT


    def calc_abg(self, obscode=500):
        '''
        [x] Transform equatorial {xyz, vxyz} (barycentric) to ecliptic {xyz, vxyz}.
        [x] Compute the (ra, dec) of the rock.
        [x] Get {xyz, vxyz} of earth in the ecliptic frame.
        [x] Transform {xyz, vxyz} from ecliptic to telescope-centric coordinates.
        [x] Use xyzT, vxyzT to compute abg coordinates.
        '''

        x_ec, y_ec, z_ec = self.x, self.y, self.z
        vx_ec, vy_ec, vz_ec = self.vx, self.vy, self.vz

        o = Observe(self, obscode=obscode)
        l = o.ecliptic_longitude
        b = o.ecliptic_latitude

        #t = ts.tt(jd=self.epoch.tt.jd)
        t = ts.ut1(jd=self.epoch.ut1.jd)
        earth = planets['earth']

        # Only used for the topocentric calculation.
        if obscode != 500:
            obscode = str(obscode).zfill(3)
            obs = observatories[observatories.obscode == obscode]
            earth += Topos(latitude_degrees=obs.lat.values, longitude_degrees=obs.lon.values, elevation_m=obs.elevation.values)

        ee = earth.at(t)

        x0, y0, z0 = ee.ecliptic_xyz().au * u.au
        vx0, vy0, vz0 = ee.ecliptic_velocity().au_per_d * u.au / u.d

        xT, yT, zT = self.ecl_to_tel(x_ec - x0, y_ec - y0, z_ec - z0, l, b)
        vxT, vyT, vzT = self.ecl_to_tel(vx_ec - vx0, vy_ec - vy0, vz_ec - vz0, l, b)

        self.gamma = 1/zT
        self.alpha = self.gamma * xT
        self.beta = self.gamma * yT
        self.dalpha = self.gamma * vxT
        self.dbeta = self.gamma * vyT
        self.dgamma = self.gamma * vzT

    def equ_to_ecl(self, x, y, z):
        '''
        This is right
        '''
        x_ec = x
        y_ec = np.cos(epsilon)*y + np.sin(epsilon)*z
        z_ec = -np.sin(epsilon)*y + np.cos(epsilon)*z
        return x_ec, y_ec, z_ec


#class SpaceRock(Transformations, Convenience):
#
#    def __init__(self, input_coordinates='keplerian', input_frame='barycentric', input_angles='degrees', *args, **kwargs):
#
#        if abg_obscode != 500:
#            SpaceRock.obscode = str(abg_obscode).zfill(3)
#            obs = observatories[observatories.obscode == SpaceRock.obscode]
#            SpaceRock.obslat = obs.lat.values
#            SpaceRock.obslon = obs.lon.values
#            SpaceRock.obselev = obs.elevation.values
#        else:
#            SpaceRock.obscode = 500
#
#        # Case-insensitive keyword arguments.
#        kwargs = {key.lower(): data for key, data in kwargs.items()}
#
#
#        input_coordinates = input_coordinates.lower()
#        input_frame = input_frame.lower()
#        input_angles = input_angles.lower()
#
#        # scalar input -> arrays
#        for idx, key in enumerate([*kwargs]):
#            if np.isscalar(kwargs.get(key)):
#                kwargs[key] = np.array([kwargs.get(key)])
#
#        if kwargs.get('delta_h') is not None:
#            self.delta_H = kwargs.get('delta_h')
#        if kwargs.get('rotation_period') is not None:
#            self.rotation_period = kwargs.get('rotation_period') * u.day
#        if kwargs.get('phi0') is not None:
#            self.phi_0 = kwargs.get('phi0')
#
#
#        self.t0 = kwargs.get('epoch') * u.day
#
#        if kwargs.get('name') is not None:
#            self.name = kwargs.get('name')
#        else:
#            # produces random, non-repeting integers between 0 and 1e10 - 1
#            self.name = np.array(['{:010}'.format(value) for value in random.sample(range(int(1e10)), len(self.t0))])
#
#        self.epoch = Time(kwargs.get('epoch'), format=input_time_format, scale=input_time_scale)
#
#        SpaceRock.frame = input_frame
#
#
#        if SpaceRock.frame == 'barycentric':
#            mu = mu_bary
#        elif SpaceRock.frame == 'heliocentric':
#            mu = mu_helio
#
#        if input_angles == 'degrees':
#            angle_unit = u.degree
#        elif input_angles == 'radians':
#            angle_unit = u.radian
#        else:
#            raise ValueError('The input_angles argument must be a string \
#                              that reads either degrees or radians.')
#
#        if input_coordinates == 'keplerian':
#
#            self.a = kwargs.get('a') * u.au
#            self.e = kwargs.get('e') * u.dimensionless_unscaled
#            self.inc = Angle(kwargs.get('inc'), angle_unit).to(u.rad)
#            self.node = Angle(kwargs.get('node'), angle_unit).to(u.rad)
#            self.arg = Angle(kwargs.get('arg'), angle_unit).to(u.rad)
#
#            if (kwargs.get('t_peri') is None) and (kwargs.get('m') is not None):
#                self.M = Angle(kwargs.get('m'), angle_unit).to(u.rad)
#                self.t_peri = self.calc_t_peri()
#
#            elif (kwargs.get('m') is None) and (kwargs.get('t_peri') is not None):
#                self.t_peri = Time(kwargs.get('t_peri'), format=input_time_format, scale=input_time_scale)
#                self.M = np.sqrt(mu / self.a**3) * (self.epoch.jd * u.day - self.t_peri.jd * u.day)
#
#            self.kep_to_xyz(mu)
#
#            if calc_abg == True:
#                if self.__class__.frame == 'heliocentric':
#                    self.to_bary()
#                    self.xyz_to_abg()
#                    self.to_helio()
#                else:
#                    self.xyz_to_abg()
#
#
#        elif input_coordinates == 'cartesian':
#
#            self.x = kwargs.get('x') * u.au
#            self.y = kwargs.get('y') * u.au
#            self.z = kwargs.get('z') * u.au
#            self.vx = kwargs.get('vx') * (u.au / u.day)
#            self.vy = kwargs.get('vy') * (u.au / u.day)
#            self.vz = kwargs.get('vz') * (u.au / u.day)
#
#            self.xyz_to_kep(mu)
#
#            if calc_abg == True:
#                if self.__class__.frame == 'heliocentric':
#                    self.to_bary()
#                    self.xyz_to_abg()
#                    self.to_helio()
#                else:
#                    self.xyz_to_abg()
#
#        self.varpi = Angle((self.arg + self.node).wrap_at(2 * np.pi * u.rad), u.rad)
#
#        if kwargs.get('h') is not None:
#            self.H = kwargs.get('h')
#            self.G = kwargs.get('g')
#            if self.G is None:
#                self.G = np.repeat(0.15, len(self))
#
#        # self.epoch = Time(self.epoch, format='jd', scale='utc')
#
    def plot_orbits(self, N, colors, alphas, ax, linewidths):
        '''
        Plot the orbits of all your rocks.
        '''
        x = np.zeros([N, len(self.a)])
        y = np.zeros([N, len(self.a)])
        z = np.zeros([N, len(self.a)])
        for idx, M in enumerate(np.linspace(0, 2*np.pi, N)):
            xx, yy, zz, _, _, _ = self.kep_to_xyz_temp(self.a, self.e, self.inc.rad,
                                                       self.arg.rad, self.node.rad,
                                                       np.repeat(M, len(self.a)) * u.rad)
            x[idx] = xx
            y[idx] = yy
            z[idx] = zz

        for idx in range(len(self.a)):
            ax.plot(x.T[idx], y.T[idx], color=colors[idx], alpha=alphas[idx], linewidth=linewidths[idx])
