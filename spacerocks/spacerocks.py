###############################################################################
# SpaceRocks, version 0.6.0
#
# Author: Kevin Napier kjnapier@umich.edu
################################################################################

import sys
import os
import random

import healpy as hp

import rebound

import numpy as np

from scipy.optimize import newton

from skyfield.api import Topos, Loader

from astropy import units as u
from astropy.table import Table
from astropy.coordinates import Angle
from astropy.constants import c

import matplotlib.pyplot as plt
from numba import jit
import pandas as pd

from .linalg3d import dot, norm, cross, euler_rotation

# Read in the observatory codes file and rehash as a dataframe.
observatories = pd.read_csv('/Users/kjnapier/research/spacerocks/data/observatories.csv')

# Load in planets for ephemeride calculation.
load = Loader('./Skyfield-Data', expire=False, verbose=False)
ts = load.timescale()
planets = load('de423.bsp')
sun = planets['sun']
jupiter = planets['jupiter barycenter']
saturn = planets['saturn barycenter']
uranus = planets['uranus barycenter']
neptune = planets['neptune barycenter']
pluto = planets['pluto barycenter']

Mercury_mass = 1.6601367952719304E-07
Venus_mass = 2.4478383396645447E-06
Earth_mass = 3.0404326462685257E-06
Mars_mass = 3.2271514450538743E-07
Sun_mass = 1 + Mercury_mass + Venus_mass + Earth_mass + Mars_mass
Jupiter_mass = 9.547919384243222E-04
Saturn_mass = 2.858859806661029E-04
Uranus_mass = 4.3662440433515637E-05
Neptune_mass = 5.151389020535497E-05
Pluto_mass = 7.361781606089469e-09

# standard gravitational parameter, GM. This value comes directly from Horizons.
mu_bary = 0.00029630927492415936 * u.radian**2 * u.au**3 / u.day**2
mu_helio = 0.00029591220828559093 * u.radian**2 * u.au**3 / u.day**2

# obliquity of the Earth. From Horizons.
epsilon = Angle(84381.448, u.arcsec).radian

# speed of light
c = c.to(u.au / u.year)

# G parameter for magnitude estimation
G = 0.15

class SpaceRock:

    def __init__(self, input_coordinates='keplerian', input_frame='barycentric',
                 input_angles='degrees', precise=False, NSIDE=None, obscode=None,
                 uncertainties=None, *args, **kwargs):

        self.frame = input_frame
        self.tau = kwargs.get('tau') * u.day

        if obscode is not None:
            self.obscode = str(obscode).zfill(3)
        else:
            self.obscode = None

        if self.frame == 'barycentric':
            mu = mu_bary
        elif self.frame == 'heliocentric':
            mu = mu_helio

        self.precise = precise

        if (self.precise is not True) and (self.precise is not False):
            raise ValueError('The parameter precise must be set to either True or False.')

        # Case-insensitive keyword arguments.
        kwargs = {key.lower(): data for key, data in kwargs.items()}
        keywords = ['a', 'e', 'inc', 'node', 'arg', 'M',
                    'x', 'y', 'z', 'vx', 'vy', 'vz',
                    'tau', 'epoch', 'h', 'name']

        if not all(key in keywords for key in [*kwargs]):
            raise ValueError('Keywords are limited to a, e, inc, node,\
                              arg, M, x, y, z, vx, vy, vz, tau, epoch,\
                              H, name')

        input_coordinates = input_coordinates.lower()
        input_frame = input_frame.lower()
        input_angles = input_angles.lower()

        # scalar input -> arrays
        for idx, key in enumerate([*kwargs]):
            if np.isscalar(kwargs.get(key)):
                kwargs[key] = np.array([kwargs.get(key)])

        if NSIDE is not None:
            if np.isscalar(NSIDE):
                NSIDE = np.array([NSIDE])
            self.NSIDE = NSIDE
        else:
            self.NSIDE = None

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

            if (kwargs.get('epoch') is None) and (kwargs.get('M') is not None):
                self.M = Angle(kwargs.get('M'), angle_unit) * u.rad
                self.epoch = self.tau - self.M / np.sqrt(mu / self.a**3)

            elif (kwargs.get('M') is None) and (kwargs.get('epoch') is not None):
                self.epoch = kwargs.get('epoch') * u.day # time at perihelion
                self.M = np.sqrt(mu / self.a**3) * (self.tau - self.epoch)

            # this looks redundant but it allows for broadcasring.
            self.tau = self.epoch +  self.M / np.sqrt(mu / self.a**3)
            self.kep_to_xyz(mu)

            if self.frame == 'barycentric':
                self.xyz_to_equa()

            elif self.frame == 'heliocentric':
                self.to_bary()
                self.xyz_to_equa()
                self.to_helio()


        elif input_coordinates == 'cartesian':

            self.x = kwargs.get('x') * u.au
            self.y = kwargs.get('y') * u.au
            self.z = kwargs.get('z') * u.au
            self.vx = kwargs.get('vx') * (u.au / u.day)
            self.vy = kwargs.get('vy') * (u.au / u.day)
            self.vz = kwargs.get('vz') * (u.au / u.day)

            self.xyz_to_kep(mu)
            self.epoch = self.tau - self.M / np.sqrt(mu / self.a**3)

            # this looks redundant but it allows for broadcasring.
            self.tau = self.epoch + self.M / np.sqrt(mu / self.a**3)

            if self.frame == 'barycentric':
                self.xyz_to_equa()

            elif self.frame == 'heliocentric':
                self.to_bary()
                self.xyz_to_equa()
                self.to_helio()

        self.varpi = (self.arg + self.node).wrap_at(2 * np.pi * u.rad)

        if kwargs.get('name') is not None:
            self.name = kwargs.get('name')
        else:
            # produces random, non-repeting integers between 0 and 1e10 - 1
            self.name = ['{:010}'.format(value) for value in random.sample(range(int(1e10)), len(self.a))]

        self.H = kwargs.get('h')
        if self.H is not None:
            self.mag = self.estimate_mag()
        if NSIDE is not None:
            for value in self.NSIDE:
                setattr(SpaceRock,
                        'HPIX_{}'.format(value),
                        self.radec_to_hpix(value))

    def calc_E(self, e, M):
        '''
        This method employs Newton's method to solve Kepler's Equation.
        '''
        f = lambda E, M, e: E - e * np.sin(E) - M
        E0 = M
        E = newton(f, E0, args=(M, e))
        return E

    @jit
    def calc_E_fast(self):
        E = self.M
        for kk in range(100):
            E = self.M + self.e * np.sin(E) * u.rad
        return E

    def xyz_to_equa(self):
        '''
        Transform from barycentric Cartesian coordinates to equatorial
        coordinates. If you need very precise values, the method will use
        a topocentric correction to the Earth's position.
        See https://en.wikipedia.org/wiki/Horizontal_coordinate_system.
        This results in a significant slowdown to the code.
        '''

        t = ts.tt(jd=self.tau.value)
        earth = planets['earth']

        # Only used for the topocentric calculation.
        if self.precise == True:
            if self.obscode is not None:
                print('doing topos')
                obs = observatories[observatories.obscode == self.obscode]
                earth += Topos(latitude_degrees=obs.lat.values,
                               longitude_degrees=obs.lon.values,
                               elevation_m=obs.elevation.values) # topocentric calculation

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
            if idx < max(range(10)):
                x0, y0, z0 = self.kep_to_xyz_pos(self.a, self.e, self.inc,
                                                 self.arg, self.node, M, self.precise)

        # Cartesian to spherical coordinate
        self.delta = norm([x, y, z])
        self.ltt = self.delta / c
        self.dec = Angle(np.arcsin(z / norm([x, y, z])), u.rad)
        self.ra = Angle(np.arctan2(y, x), u.rad).wrap_at(2 * np.pi * u.rad)
        self.phase_angle = Angle(np.arccos(-(earth_dis**2 - self.r**2 - self.delta**2)/(2 * self.r* self.delta)), u.rad)
        self.elong = Angle(np.arccos(-(self.r**2 - self.delta**2 - earth_dis**2)/(2 * self.delta * earth_dis)), u.rad)


        return self


    def radec_to_hpix(self, NSIDE):
        '''
        Convert (ra, dec) into healpix
        '''
        return hp.pixelfunc.ang2pix(NSIDE, np.pi/2 - self.dec.radian, self.ra.radian, nest=True)

    def estimate_mag(self):
        '''
        Estimate the apparent magnitude of a TNO
        '''
        # H + 5*log10(delta) + 5*log10(r) - 2.5*log10((1-G)*phi1 + G*phi2)
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
        anomaly, E. If you need very precise coordinates, the method uses Newton's
        root-finding method from scipy.optimize. This is very precise, but it is
        not vectorizable or compilable, so it is rather slow
        (> factor of 5 slowdown). Otherwise, the .ethod uses a fixed-point
        iteration method. See https://en.wikipedia.org/wiki/Kepler%27s_equation
        '''
        # compute eccentric anomaly E
        if self.precise == True:
            E = np.array(list(map(self.calc_E, self.e.value, self.M.value))) * u.rad
        else:
            E = self.calc_E_fast()

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
        if self.frame == 'heliocentric':
            t = ts.tt(jd=self.tau.value) #37 leap seconds
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
            self.epoch[self.M.value >= np.pi] = self.tau + (2*np.pi * u.rad - self.M) / np.sqrt(mu_bary / self.a**3)
            self.epoch[self.M.value < np.pi] = self.tau - self.M / np.sqrt(mu_bary / self.a**3)
            self.frame = 'barycentric'

    def to_helio(self):
        '''
        Method to convert barycentric coordinates to heliocentric coordinates.
        '''
        if self.frame == 'barycentric':

            t = ts.tt(jd=self.tau.value)
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
            self.epoch[self.M.value >= np.pi] = self.tau + (2*np.pi * u.rad - self.M) / np.sqrt(mu_helio / self.a**3)
            self.epoch[self.M.value < np.pi] = self.tau - self.M / np.sqrt(mu_helio / self.a**3)
            self.frame = 'heliocentric'

        return self

    def kep_to_xyz_pos(self, a, e, inc, arg, node, M, precision):
        # compute eccentric anomaly E
        if precision == True:
            E = np.array(list(map(self.calc_E, e.value, M.value))) * u.rad
        else:
            E = self.calc_E_fast()

        # compute true anomaly v
        ν = 2 * np.arctan2((1 + e)**0.5*np.sin(E/2.), (1 - e)**0.5*np.cos(E/2.))

        # compute the distance to the central body r
        r = a * (1 - e * np.cos(E))

        # obtain the position o and velocity ov vector
        o = [r * np.cos(ν), r * np.sin(ν), np.zeros(len(ν))]

        # Rotate o to the inertial frame
        x, y, z = euler_rotation(arg, inc, node, o) * u.au

        return x, y, z

    def get_dict(self):
        '''
        Create a dictionary of the object attributes. This method is used
        by the pandas_df and astropy_table methods.
        '''
        data = {'name':self.name,
                'a':self.a,
                'e':self.e,
                'inc':self.inc,
                'arg':self.arg,
                'node':self.node,
                'varpi':self.varpi,
                'epoch':self.epoch,
                'M':self.M,
                'tau':self.tau,
                'x':self.x,
                'y':self.y,
                'z':self.z,
                'vx':self.vx,
                'vy':self.vy,
                'vz':self.vz,
                'ra':self.ra,
                'dec':self.dec,
                'delta':self.delta,
                'ltt':self.ltt,
                'phase_angle':self.phase_angle,
                'elong':self.elong,
                'r':self.r}

        if self.NSIDE is not None:
            for value in self.NSIDE:
                data['HPIX_{}'.format(value)] = getattr(SpaceRock, 'HPIX_{}'.format(value))

        if self.H is not None:
            data['H'] = self.H
            data['mag'] = self.mag

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

    def plot_radec(self, color='black', alpha=0.5):
        '''
        Plot the right ascension and declination of each object on a Mollweide
        projection. (See https://en.wikipedia.org/wiki/Mollweide_projection)
        '''
        fig = plt.figure(figsize=(12, 8))
        ax = fig.add_subplot(111, projection="mollweide")

        xdata = self.ra
        xdata[xdata.value > np.pi] -= 2*np.pi * u.rad

        plt.scatter(xdata, self.dec, color=color, alpha=alpha)
        return


    def set_simulation(self, startdate):

        t = ts.tt(jd=startdate)

        active_bodies = [sun, jupiter, saturn, uranus, neptune, pluto]
        x, y, z = np.array([body.at(t).ecliptic_xyz().au for body in active_bodies]).T
        vx, vy, vz = np.array([body.at(t).ecliptic_velocity().au_per_d for body in active_bodies]).T
        names = ['Sun', 'Jupiter', 'Saturn', 'Uranus', 'Neptune', 'Pluto']

        # create a dataframe of the massive bodies in the solar system
        ss = pd.DataFrame()
        ss['x'] = x
        ss['y'] = y
        ss['z'] = z
        ss['vx'] = vx
        ss['vy'] = vy
        ss['vz'] = vz
        ss['mass'] = [Sun_mass, Jupiter_mass, Saturn_mass, Uranus_mass, Neptune_mass, Pluto_mass]
        ss['a'] = 1 / (2 / norm([ss.x, ss.y, ss.z]) - norm([ss.vx, ss.vy, ss.vz])**2 / mu_bary.value)
        ss['hill_radius'] = ss.a * pow(ss.mass/(3.*Sun_mass),1./3.)
        ss['name'] = names

        sim = rebound.Simulation()
        sim.units = ('day', 'AU', 'Msun')

        for p in ss.itertuples():
            sim.add(x=p.x, y=p.y, z=p.z,
                    vx=p.vx, vy=p.vy, vz=p.vz,
                    m=p.mass, hash=p.name, r=p.hill_radius)

        sim.move_to_com()

        sim.N_active = len(ss)
        sim.testparticle_type = 0

        sim.integrator = 'mercurius'
        sim.dt = 1 # one day

        sim.ri_ias15.min_dt = sim.dt / 1440 # one minute
        sim.ri_mercurius.hillfac = 3
        #sim.ri_whfast.safe_mode = 0
        #sim.ri_whfast.corrector = 11

        return sim

    def propagate(self, enddates):
        '''
        Integrate the bodies to the desired date. The logic could be cleaner
        but it should work.
        '''

        if np.isscalar(enddates):
            enddates = np.array([enddates])

        Nx = len(enddates)
        Ny = len(self.x)
        x_values = np.zeros([Nx, Ny])
        y_values = np.zeros([Nx, Ny])
        z_values = np.zeros([Nx, Ny])
        vx_values = np.zeros([Nx, Ny])
        vy_values = np.zeros([Nx, Ny])
        vz_values = np.zeros([Nx, Ny])
        tau_values = np.zeros([Nx, Ny])
        name_values = np.tile(self.name, Nx)
        if self.H is not None:
            H_values = np.tile(self.H, Nx)

        in_frame = self.frame

        # We need to (or should) work in barycentric coordinates in rebound
        if in_frame == 'heliocentric':
            self.to_bary()

        # Rehash as a dataframe for easy access
        df = self.pandas_df()

        # Integrate all particles to the same tau
        pickup_times = df.tau.values
        sim = self.set_simulation(np.min(pickup_times))
        sim.t = np.min(df.tau.values)

        for time in np.sort(np.unique(pickup_times)):
            ps = df[df.tau.values == time]
            for p in ps.itertuples():
                sim.add(x=p.x, y=p.y, z=p.z,
                vx=p.vx, vy=p.vy, vz=p.vz,
                hash=p.name)
                sim.integrate(time, exact_finish_time=1)

        for ii, time in enumerate(np.sort(enddates)):
            sim.integrate(time, exact_finish_time=1)
            for jj, name in enumerate(self.name):
                x_values[ii, jj] = sim.particles[name].x
                y_values[ii, jj] = sim.particles[name].y
                z_values[ii, jj] = sim.particles[name].z
                vx_values[ii, jj] = sim.particles[name].vx
                vy_values[ii, jj] = sim.particles[name].vy
                vz_values[ii, jj] = sim.particles[name].vz
                tau_values[ii, jj] = sim.t

        self.x = x_values.flatten() * u.au
        self.y = y_values.flatten() * u.au
        self.z = z_values.flatten() * u.au
        self.vx = vx_values.flatten() * (u.au / u.day)
        self.vy = vy_values.flatten() * (u.au / u.day)
        self.vz = vz_values.flatten() * (u.au / u.day)
        self.name = name_values.flatten()
        self.tau = tau_values.flatten() * u.day

        self.xyz_to_kep(mu_bary)

        self.epoch = np.zeros(Nx * Ny)
        lp = self.M < np.pi * u.rad
        self.epoch[lp] = self.tau[lp] - self.M[lp] / np.sqrt(mu_bary / self.a[lp]**3)
        self.epoch[~lp] = self.tau[~lp] + (2*np.pi * u.rad - self.M[~lp]) / np.sqrt(mu_bary / self.a[~lp]**3)

        self.xyz_to_equa()
        self.varpi = (self.arg + self.node).wrap_at(2 * np.pi * u.rad)

        if self.H is not None:
            self.H = H_values
            self.mag = self.estimate_mag()

        # calculate new hpix
        if self.NSIDE is not None:
            for value in self.NSIDE:
                setattr(SpaceRock,
                'HPIX_{}'.format(value),
                self.radec_to_hpix(value))

        # be polite and return orbital parameters in the input frame.
        if in_frame == 'heliocentric':
            self.to_helio()

        out_df = self.pandas_df()

        return df, out_df
