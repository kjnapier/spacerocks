###############################################################################
# SpaceRocks, version 0.7.0
#
# Author: Kevin Napier kjnapier@umich.edu
################################################################################

import sys
import os
import random
import copy

import healpy as hp

import rebound
import reboundx
from reboundx import constants

import numpy as np
import pandas as pd
from numba import jit
import ephem
from skyfield.api import Topos, Loader

from astropy import units as u
from astropy.table import Table
from astropy.coordinates import Angle
from astropy.time import Time
from astropy.coordinates import SkyCoord

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from mpl_toolkits.mplot3d import Axes3D

try:
    import cartopy.crs as ccrs
except:
    pass

from .linalg3d import *
from .constants import *
#from .jacobians import *

# Read in the observatory codes file and rehash as a dataframe.
observatories = pd.read_csv(os.path.join(os.path.dirname(__file__),
                            'data',
                            'observatories.csv'))

# Load in planets for ephemeride calculation.
load = Loader('./Skyfield-Data', expire=False, verbose=False)
ts = load.timescale()
planets = load('de423.bsp')
sun = planets['sun']


class SpaceRock:

    def __init__(self, input_coordinates='keplerian', input_frame='barycentric',
                 input_angles='degrees', input_time_format='jd', calc_equa=True,
                 input_time_scale='utc', NSIDE=None, obscode=None,
                 uncertainties=None, *args, **kwargs):

        # Case-insensitive keyword arguments.
        kwargs = {key.lower(): data for key, data in kwargs.items()}
        keywords = ['a', 'e', 'inc', 'node', 'arg', 'm',
                    'x', 'y', 'z', 'vx', 'vy', 'vz',
                    'obsdate', 't_peri', 'h', 'name']

        if not all(key in keywords for key in [*kwargs]):
            raise ValueError('Keywords are limited to a, e, inc, node,\
                              arg, m, x, y, z, vx, vy, vz, obsdate, t_peri,\
                              H, name')

        input_coordinates = input_coordinates.lower()
        input_frame = input_frame.lower()
        input_angles = input_angles.lower()

        # scalar input -> arrays
        for idx, key in enumerate([*kwargs]):
            if np.isscalar(kwargs.get(key)):
                kwargs[key] = np.array([kwargs.get(key)])

        self.obsdate = Time(kwargs.get('obsdate'),
                        format=input_time_format,
                        scale=input_time_scale).jd * u.day

        self.H = kwargs.get('h')

        if NSIDE is not None:
            if np.isscalar(NSIDE):
                NSIDE = np.array([NSIDE])
            SpaceRock.NSIDE = NSIDE
        else:
            SpaceRock.NSIDE = None

        SpaceRock.calc_equa = calc_equa

        # Base attributes. Required by every object.
        attributes = ['a', 'e', 'inc', 'node', 'arg', 'M', 't_peri', 'varpi',
                      'x', 'y', 'z', 'vx', 'vy', 'vz', 'obsdate', 'name', 'r']

        if SpaceRock.calc_equa == True:

            for attr in ['ra', 'dec', 'skycoord', 'elong',
                         'delta', 'ltt', 'phase_angle']:
                attributes.append(attr)

            if self.H is not None:
                attributes.append('H')
                attributes.append('mag')

            if SpaceRock.NSIDE is not None:
                for value in SpaceRock.NSIDE:
                    attributes.append('HPIX_{}'.format(value))

        SpaceRock.attributes = attributes
        SpaceRock.frame = input_frame

        if obscode is not None:
            SpaceRock.obscode = str(obscode).zfill(3)
            obs = observatories[observatories.obscode == SpaceRock.obscode]
            SpaceRock.obslat = obs.lat.values
            SpaceRock.obslon = obs.lon.values
            SpaceRock.obselev = obs.elevation.values
        else:
            SpaceRock.obscode = None

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
                lp = self.M < np.pi * u.rad
                self.t_peri = np.zeros(len(self.obsdate))
                self.t_peri[lp] = self.obsdate.value[lp] - self.M.value[lp] \
                                  / np.sqrt(mu_bary.value / self.a.value[lp]**3)
                self.t_peri[~lp] = self.obsdate.value[~lp] \
                                   + (2*np.pi - self.M.value[~lp]) \
                                   / np.sqrt(mu_bary.value / self.a.value[~lp]**3)
                self.t_peri = Time(self.t_peri, format='jd', scale='utc')

            elif (kwargs.get('m') is None) and (kwargs.get('t_peri') is not None):
                self.t_peri = Time(kwargs.get('t_peri'),
                                  format=input_time_format,
                                  scale=input_time_scale)
                self.M = np.sqrt(mu / self.a**3) * (self.obsdate - self.t_peri.jd * u.day)

            # this looks redundant but it allows for broadcasring.
            # self.obsdate = self.t_peri + self.M / np.sqrt(mu / self.a**3)
            self.kep_to_xyz(mu)

            if SpaceRock.calc_equa == True:

                if SpaceRock.frame == 'barycentric':
                    self.xyz_to_equa()

                elif SpaceRock.frame == 'heliocentric':
                    self.to_bary()
                    self.xyz_to_equa()
                    self.to_helio()

                if self.H is not None:
                    self.mag = self.estimate_mag()

                if NSIDE is not None:
                    for value in SpaceRock.NSIDE:
                        setattr(SpaceRock,
                                'HPIX_{}'.format(value),
                                self.radec_to_hpix(value))


        elif input_coordinates == 'cartesian':

            self.x = kwargs.get('x') * u.au
            self.y = kwargs.get('y') * u.au
            self.z = kwargs.get('z') * u.au
            self.vx = kwargs.get('vx') * (u.au / u.day)
            self.vy = kwargs.get('vy') * (u.au / u.day)
            self.vz = kwargs.get('vz') * (u.au / u.day)

            self.xyz_to_kep(mu)
            lp = self.M < np.pi * u.rad
            self.t_peri[lp] = self.obsdate.jd[lp] * u.day - self.M[lp] \
                              / np.sqrt(mu_bary / self.a[lp]**3)
            self.t_peri[~lp] = self.obsdate.jd[~lp] * u.day \
                               + (2*np.pi * u.rad - self.M[~lp]) \
                               / np.sqrt(mu_bary / self.a[~lp]**3)
            self.t_peri = Time(self.t_peri, format='jd', scale='utc')

            # this looks redundant but it allows for broadcasring.
            # self.obsdate = self.t_peri + self.M / np.sqrt(mu / self.a**3)

            if SpaceRock.calc_equa == True:

                if SpaceRock.frame == 'barycentric':
                    self.xyz_to_equa()

                elif SpaceRock.frame == 'heliocentric':
                    self.to_bary()
                    self.xyz_to_equa()
                    self.to_helio()

                self.H = kwargs.get('h')
                if self.H is not None:
                    self.mag = self.estimate_mag()

                if NSIDE is not None:
                    for value in SpaceRock.NSIDE:
                        setattr(SpaceRock,
                                'HPIX_{}'.format(value),
                            self.radec_to_hpix(value))

        self.varpi = Angle((self.arg + self.node).wrap_at(2 * np.pi * u.rad), u.rad)

        if kwargs.get('name') is not None:
            self.name = kwargs.get('name')
        else:
            # produces random, non-repeting integers between 0 and 1e10 - 1
            self.name = ['{:010}'.format(value) for value in random.sample(range(int(1e10)), len(self.a))]

        self.obsdate = Time(self.obsdate, format='jd', scale='utc')

    def __len__(self):
        '''
        This method allows you to use the len() function on a SpaceRocks object.
        '''
        return len(self.name)


    def __getitem__(self, idx):
        '''
        This method allows you to index a SpaceRocks object.
        '''
        p = copy.copy(self)
        for attr in SpaceRock.attributes:
            setattr(p, attr, getattr(self, attr)[idx])

        return p


    def calc_E(self, e, M):
        '''
        This method solves Kepler's Equation with the algorithm desscribed here:
        https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19950021346.pdf
        '''
        M[M > np.pi] -= 2*np.pi
        α = (3 * np.pi**2 + 1.6 * (np.pi**2 - np.pi * abs(M))/(1 + e))/(np.pi**2 - 6)
        d = 3 * (1 - e) + α * e
        q = 2 * α * d * (1 - e) - M**2
        r = 3 * α * d * (d - 1 + e) * M + M**3
        w = (abs(r) + np.sqrt(q**3 + r**2))**(2/3)
        E1 = (2 * r * w / (w**2 + w*q + q**2) + M)/d
        f2 = e * np.sin(E1)
        f3 = e * np.cos(E1)
        f0 = E1 - f2 - M
        f1 = 1 - f3
        δ3 = -f0 / (f1 - f0 * f2 / (2 * f1))
        δ4 = -f0 / (f1 + f2 * δ3 / 2 + δ3**2 * f3/6)
        δ5 = -f0 / (f1 + δ4*f2/2 + δ4**2*f3/6 - δ4**3*f2/24)

        E = E1 + δ5
        E[E < 0] += 2 * np.pi
        return E


    def xyz_to_equa(self):
        '''
        Transform from barycentric Cartesian coordinates to equatorial
        coordinates. If you need very precise values, the method will use
        a topocentric correction to the Earth's position.
        See https://en.wikipedia.org/wiki/Horizontal_coordinate_system.
        This results in a significant slowdown to the code.
        '''

        t = ts.tt(jd=self.obsdate.value)
        earth = planets['earth']

        # Only used for the topocentric calculation.
        if SpaceRock.obscode is not None:
            earth += Topos(latitude_degrees=SpaceRock.obslat,
                           longitude_degrees=SpaceRock.obslon,
                           elevation_m=SpaceRock.obselev) # topocentric calculation

        x_earth, y_earth, z_earth = earth.at(t).position.au * u.au # earth ICRS position
        earth_dis = norm([x_earth, y_earth, z_earth])
        x0, y0, z0 = self.x, self.y, self.z
        for idx in range(10):
            # transfer ecliptic to ICRS and shift to Geocentric (topocentric)
            x = x0 - x_earth
            y = y0 * np.cos(epsilon) - z0 * np.sin(epsilon) - y_earth
            z = y0 * np.sin(epsilon) + z0 * np.cos(epsilon) - z_earth
            delta = norm([x, y, z])
            ltt = delta / c
            M = self.M - ltt * (mu_bary / self.a**3)**0.5
            if idx < 9:
                x0, y0, z0 = self.kep_to_xyz_pos(self.a, self.e, self.inc,
                                                 self.arg, self.node, M)

        # Cartesian to spherical coordinate
        self.delta = norm([x, y, z])
        self.ltt = self.delta / c
        self.dec = Angle(np.arcsin(z / norm([x, y, z])), u.rad)
        self.ra = Angle(np.arctan2(y, x), u.rad).wrap_at(2 * np.pi * u.rad)
        self.phase_angle = Angle(np.arccos(-(earth_dis**2 - self.r**2 - self.delta**2)/(2 * self.r* self.delta)), u.rad)
        self.elong = Angle(np.arccos(-(self.r**2 - self.delta**2 - earth_dis**2)/(2 * self.delta * earth_dis)), u.rad)
        self.skycoord = SkyCoord(self.ra, self.dec, frame='icrs')

        return self

    #def sky_error(self):

    #    kep_to_xyz_jac = kep_to_xyz_jacobian(self.x.value,
    #                                         self.y.value,
    #                                         self.z.value,
    #                                         self.vx.value,
    #                                         self.vy.value,
    #                                         self.vz.value,
    #                                         mu_bary.value)
    #    J = np.linalg.inv(kep_to_xyz_jac)
    #    cov_xyz = np.matmul(J, np.matmul(self.cov_kep, J.T))


    #    return self


    def radec_to_hpix(self, NSIDE):
        '''
        Convert (ra, dec) into healpix.
        '''
        return hp.pixelfunc.ang2pix(NSIDE, np.pi/2 - self.dec.radian, self.ra.radian, nest=True)


    def estimate_mag(self):
        '''
        Estimate the apparent magnitude of a TNO
        '''
        q = (self.r**2 + self.delta**2 - 1 * u.au**2)/(2 * self.r * self.delta)

        ## pyephem
        beta = np.zeros(len(q))
        beta[np.where(q <= 1)[0]] = np.pi * u.rad
        beta[np.where(q >= 1)[0]] = 0 * u.rad

        Psi_1 = np.exp(-3.33 * np.tan(beta/2)**0.63)
        Psi_2 = np.exp(-1.87 * np.tan(beta/2)**1.22)
        mag = self.H + 5 * np.log10(self.r * self.delta / u.au**2)

        not_zero = np.where((Psi_1 != 0) | (Psi_2 != 0))[0]
        mag[not_zero] -= 2.5 * np.log10((1 - G) * Psi_1[not_zero] + G * Psi_2[not_zero])

        return mag


    def kep_to_xyz(self, mu):
        '''
        Transform from Keplerian to cartesian coordinates. There is no analytic
        solution to solve Kepler's Equation M = E - eSin[e] for the eccentric
        anomaly (E), but I use a close approximation.
        '''
        # compute eccentric anomaly E
        E = self.calc_E(self.e.value, self.M.value) * u.rad

        # compute true anomaly v
        ν = 2 * np.arctan2((1 + self.e)**0.5*np.sin(E/2.), (1 - self.e)**0.5*np.cos(E/2.))

        # compute the distance to the central body r
        r = self.a * (1 - self.e*np.cos(E))

        # obtain the position o and velocity ov vector
        o = [r * np.cos(ν), r * np.sin(ν), np.zeros(len(ν))]
        ov = [(mu * self.a)**0.5 / r * (-np.sin(E)),
              (mu * self.a)**0.5 / r * ((1-self.e**2)**0.5 * np.cos(E)),
              np.zeros(len(ν))]

        # Rotate o and ov to the inertial frame
        self.x, self.y, self.z = euler_rotation(self.arg, self.inc, self.node, o) * u.au
        self.vx, self.vy, self.vz = euler_rotation(self.arg, self.inc, self.node, ov) * u.au / u.day
        self.r = norm([self.x, self.y, self.z])
        return self


    def xyz_to_kep(self, mu):
        '''
        Transform from Cartesian to Keplerian coordinates. The units in this method
        are a bit hacky due to unexpected behavior and known issues with astropy
        units. But it returns the correct values with the correct units.
        '''
        ## there are a lot of hacky units in here because astropy doesn't behave as expected.
        # compute the barycentric distance r
        self.r = norm([self.x, self.y, self.z])
        rrdot = dot([self.x, self.y, self.z], [self.vx, self.vy, self.vz])

        # compute the specific angular momentum h
        hx, hy, hz = cross([self.x, self.y, self.z], [self.vx, self.vy, self.vz])
        h = norm([hx, hy, hz])

        # compute eccentricity vector
        ### hacky units
        ex, ey, ez = u.au**3 * u.rad**2/ u.day**2 * np.array(cross([self.vx, self.vy, self.vz], \
                     [hx, hy, hz])) / mu  - [self.x, self.y, self.z]*u.au/self.r
        self.e = norm([ex, ey, ez])

        # compute vector n
        ### hacky units
        nx, ny, nz = -hy.value, hx.value, np.zeros(len(hz))
        n = norm([nx, ny, nz])

        # compute true anomaly ν, the angle between e and r
        ν = np.arccos(dot([ex, ey, ez], [self.x, self.y, self.z]) / (self.e*self.r))
        ν[rrdot < 0] = 2 * np.pi * u.rad - ν[rrdot < 0]

        # compute inclination
        self.inc = Angle(np.arccos(hz/h), u.rad)

        # compute eccentric anomaly E
        E = 2 * np.arctan2(np.sqrt(1-self.e) * np.sin(ν/2), np.sqrt(1+self.e) * np.cos(ν/2))

        # compute ascending node
        node = np.arccos(nx/n)
        node[ny < 0] = 2 * np.pi - node[ny < 0]
        self.node = Angle(node, u.rad)

        # compute argument of periapsis, the angle between e and n
        arg = np.arccos(dot([nx, ny, nz], [ex, ey, ez]) / (n*self.e))
        arg[ez < 0] = 2 * np.pi * u.rad - arg[ez < 0]
        self.arg = Angle(arg, u.rad)

        # compute mean anomaly
        M = E - self.e * np.sin(E) * u.rad
        M[M < 0] += 2 * np.pi * u.rad
        self.M = Angle(M, u.rad)

        # compute a
        self.a = 1 / (2 / self.r - norm([self.vx, self.vy, self.vz])**2 / mu * u.rad**2)

        return self


    def to_bary(self):
        '''
        Method to convert heliocentric coordinates to barycentric coordinates.
        '''
        if SpaceRock.frame == 'heliocentric':
            t = ts.tt(jd=self.obsdate.value)
            x_sun, y_sun, z_sun = sun.at(t).ecliptic_xyz().au * u.au
            vx_sun, vy_sun, vz_sun = sun.at(t).ecliptic_velocity().au_per_d * u.au / u.day
            # calculate the barycentric xyz postion
            self.x += x_sun
            self.y += y_sun
            self.z += z_sun
            self.vx += vx_sun
            self.vy += vy_sun
            self.vz += vz_sun

            # calculate barycentric keplerian elements
            self.xyz_to_kep(mu_bary)
            self.varpi = (self.arg + self.node).wrap_at(2 * np.pi * u.rad)
            lp = self.M < np.pi * u.rad
            self.t_peri.jd[lp] = self.obsdate.value[lp] - self.M.value[lp] \
                                 / np.sqrt(mu_bary.value / self.a.value[lp]**3)
            self.t_peri.jd[~lp] = self.obsdate.value[~lp] \
                                  + (2*np.pi - self.M.value[~lp]) \
                                  / np.sqrt(mu_bary.value / self.a.value[~lp]**3)
            self.t_peri = Time(self.t_peri, format='jd', scale='utc')
            SpaceRock.frame = 'barycentric'

        return self


    def to_helio(self):
        '''
        Method to convert barycentric coordinates to heliocentric coordinates.
        '''
        if SpaceRock.frame == 'barycentric':
            t = ts.tt(jd=self.obsdate.value)
            x_sun, y_sun, z_sun = sun.at(t).ecliptic_xyz().au * u.au
            vx_sun, vy_sun, vz_sun = sun.at(t).ecliptic_velocity().au_per_d * u.au / u.day
            # calculate the heliocentric xyz postion
            self.x -= x_sun
            self.y -= y_sun
            self.z -= z_sun
            self.vx -= vx_sun
            self.vy -= vy_sun
            self.vz -= vz_sun

            # calculate heliocentric keplerian elements
            self.xyz_to_kep(mu_helio)

            self.varpi = (self.arg + self.node).wrap_at(2 * np.pi * u.rad)
            lp = self.M < np.pi * u.rad
            self.t_peri.jd[lp] = self.obsdate.value[lp] - self.M.value[lp] \
                                 / np.sqrt(mu_bary.value / self.a.value[lp]**3)
            self.t_peri.jd[~lp] = self.obsdate.value[~lp] \
                                  + (2*np.pi - self.M.value[~lp]) \
                                  / np.sqrt(mu_bary.value / self.a.value[~lp]**3)
            self.t_peri = Time(self.t_peri, format='jd', scale='utc')
            SpaceRock.frame = 'heliocentric'

        return self


    def kep_to_xyz_pos(self, a, e, inc, arg, node, M):
        '''
        Just compute the xyz position of an object. Used for iterative equatorial
        calculation.
        '''
        # compute eccentric anomaly E
        E = self.calc_E(e.value, M.value) * u.rad

        # compute true anomaly ν
        ν = 2 * np.arctan2((1 + e)**0.5*np.sin(E/2), (1 - e)**0.5*np.cos(E/2))

        # compute the distance to the central body r
        r = a * (1 - e * np.cos(E))

        # obtain the position vector o
        o = [r * np.cos(ν), r * np.sin(ν), np.zeros(len(ν))]

        # Rotate o to the inertial frame
        x, y, z = euler_rotation(arg, inc, node, o) * u.au

        return x, y, z


    def plot_orbits(self):
        '''
        Plot the orbits of all your rocks.
        '''
        x = np.zeros([500, len(self.a)])
        y = np.zeros([500, len(self.a)])
        z = np.zeros([500, len(self.a)])
        for idx, M in enumerate(np.linspace(0, 2*np.pi, 500)):
            xx, yy, zz = self.kep_to_xyz_pos(self.a, self.e, self.inc.rad,
                                        self.arg.rad, self.node.rad,
                                        np.repeat(M, len(self.a)) * u.rad, False)
            x[idx] = xx
            y[idx] = yy
            z[idx] = zz

        fig = plt.figure(figsize=(12, 12))
        ax = fig.add_subplot(111, projection='3d')
        for idx in range(len(self.a)):
            ax.plot(x.T[idx], y.T[idx], z.T[idx], color=np.random.choice(['#FFCB05', '#00274C']))

        ax.axis('off')
        ax.view_init(90, 0)

        return fig, ax


    def get_dict(self):
        '''
        Create a dictionary of the object attributes. This method is used
        by the pandas_df and astropy_table methods.
        '''
        data = {}
        for attr in SpaceRock.attributes:
            data[attr] = getattr(self, attr)

        return data


    def pandas_df(self):
        '''
        Write the rocks to a pandas dataframe. Pandas can't handle astropy
        units (yet), so if you want to keep units intact you'll have to use
        an Astropy Table.
        '''
        return pd.DataFrame(self.get_dict())


    def astropy_table(self):
        '''
        Write the rocks to an astropy table. This can handle units, though
        it is generally less elegant than pandas.
        '''
        return Table(self.get_dict())


    def write_to_csv(self, path):
        '''
        Write the data to a csv.
        '''
        df = self.pandas_df()
        df.to_csv(path)
        return 'Data written to {}.csv.'.format(path)


    def plot_radec(self, color='black', alpha=0.5, zoom=False, galactic_plane=False, ecliptic_plane=True):
        '''
        Plot the right ascension and declination of each object on a Mollweide
        projection. (See https://en.wikipedia.org/wiki/Mollweide_projection)
        '''
        if zoom == True:

            fig = plt.figure(figsize=(12, 8))
            ax = fig.add_subplot(111, projection=ccrs.PlateCarree())

            xdata = self.ra.degree
            xdata[xdata > 180] -= 360

            xmin=np.min(xdata)
            xmax=np.max(xdata)
            ymin=np.min(self.dec.degree)
            ymax=np.max(self.dec.degree)

            ax.scatter(-xdata, self.dec.degree, color='black', alpha=0.5)

            xticks = np.linspace(-xmax, -xmin, 8)
            yticks = np.linspace(ymin, ymax, 8)

            ax.set_xticks(xticks)
            ax.set_yticks(yticks)
            xticklabels = [r'${:.2f}\degree$'.format(-value) for value in xticks]
            yticklabels = [r'${:.2f}\degree$'.format(value) for value in yticks]
            ax.set_xticklabels(xticklabels)
            ax.set_yticklabels(yticklabels)

            gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                              linewidth=1, color='gray', alpha=0.5, linestyle='-')
            gl.xlocator = mticker.FixedLocator(xticks)
            gl.ylocator = mticker.FixedLocator(yticks)

            gl.bottom_labels = False
            gl.top_labels = False
            gl.left_labels = False
            gl.right_labels = False

            xrange = xmax - xmin
            yrange = ymax - ymin

            try:
                ax.set_extent([-xmax - xrange * 0.05, -xmin + xrange * 0.05,
                               ymin - yrange * 0.05, ymax + yrange * 0.05], crs=ccrs.PlateCarree())
            except:
                ax.set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())

            if ecliptic_plane == True:
                def radec2project(ra, dec):
                    ra[ra>180] -= 360
                    return (ra, dec)

                plane_lon = np.linspace(-np.pi, np.pi, 1000)
                plane_lat = np.zeros(len(plane_lon))
                ecl_plane = np.zeros([len(plane_lat), 2])

                for i in range(len(plane_lat)):
                    ecl_plane[i] = ephem.Equatorial(ephem.Ecliptic(plane_lon[i], plane_lat[i])).get()

                x, y = radec2project(np.degrees(ecl_plane.T[0]), np.degrees(ecl_plane.T[1]))
                ax.plot(x[1:999], -y[1:999], 'r-', zorder=2)

            if galactic_plane == True:

                galactic = SkyCoord(l=np.linspace(0, 360, 1000)*u.degree, b=np.zeros(1000)*u.degree, frame='galactic')
                galxdata = galactic.icrs.ra.degree
                galxdata[galxdata > 180] -= 360
                galydata = galactic.icrs.dec.degree
                order = galxdata.argsort()
                ax.plot(-galxdata[order][1:999], -galydata[order][1:999], 'g-', zorder=3)

        else:

            fig = plt.figure(figsize=(12, 8))
            ax = fig.add_subplot(111, projection='mollweide')
            ax.grid(True)

            xdata = self.ra
            xdata[xdata.value > np.pi] -= 2*np.pi * u.rad
            ax.set_xticklabels([r'$150\degree$', r'$120\degree$', r'$90\degree$' ,
                                r'$60\degree$' , r'$30\degree$' , r'$0\degree$'  ,
                                r'$330\degree$', r'$300\degree$', r'$270\degree$',
                                r'$240\degree$', r'$210\degree$'])

            ax.scatter(-xdata, self.dec, color=color, alpha=alpha)

            if ecliptic_plane == True:
                def radec2project(ra, dec):
                    ra[ra>180] -= 360
                    return (ra, dec)

                plane_lon = np.linspace(-np.pi, np.pi, 1000)
                plane_lat = np.zeros(len(plane_lon))
                ecl_plane = np.zeros([len(plane_lat), 2])

                for i in range(len(plane_lat)):
                    ecl_plane[i] = ephem.Equatorial(ephem.Ecliptic(plane_lon[i], plane_lat[i])).get()

                x, y = radec2project(np.degrees(ecl_plane.T[0]), np.degrees(ecl_plane.T[1]))
                ax.plot(np.radians(x[1:999]), -np.radians(y[1:999]), 'r-', zorder=2)

            if galactic_plane == True:

                galactic = SkyCoord(l=np.linspace(0, 360, 1000)*u.degree, b=np.zeros(1000)*u.degree, frame='galactic')
                galxdata = galactic.icrs.ra.rad
                galxdata[galxdata > np.pi] -= 2*np.pi
                galydata = galactic.icrs.dec.rad
                order = galxdata.argsort()
                ax.plot(-galxdata[order][1:999], galydata[order][1:999], 'g-', zorder=3)

        return fig, ax
