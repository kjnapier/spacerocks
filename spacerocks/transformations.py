import warnings

import numpy as np
from astropy import units as u


from .linalg3d import *
from .constants import *

from astropy.time import Time
from skyfield.api import Topos, Loader
# Load in planets for ephemeride calculation.
load = Loader('./Skyfield-Data', expire=False, verbose=False)
ts = load.timescale()
planets = load('de423.bsp')
sun = planets['sun']


class Transformations:

    def kep_to_xyz(self, mu):
        '''
        Transform from Keplerian to cartesian coordinates. There is no analytic
        solution to solve Kepler's Equation M = E - eSin[e] for the eccentric
        anomaly (E), but I use a close approximation.
        '''
        # compute eccentric anomaly E
        E = self.calc_E(self.e.value, self.M.value) * u.rad

        # compute true anomaly v
        ν = 2 * np.arctan2((1 + self.e)**0.5 * np.sin(E/2), (1 - self.e)**0.5 * np.cos(E/2))

        # compute the distance to the central body r
        r = self.a * (1 - self.e * np.cos(E))

        # obtain the position o and velocity ov vector
        # the /u.rad is to fix units...
        o = [r * np.cos(ν), r * np.sin(ν), np.zeros(len(ν))]
        ov = [(mu * self.a)**0.5 / r * (-np.sin(E))/ u.rad,
              (mu * self.a)**0.5 / r * ((1-self.e**2)**0.5 * np.cos(E))/ u.rad,
              np.zeros(len(ν))/ u.rad]

        # Rotate o and ov to the inertial frame
        self.x, self.y, self.z = euler_rotation(self.arg, self.inc, self.node, o) #* u.au
        self.vx, self.vy, self.vz = euler_rotation(self.arg, self.inc, self.node, ov) #* u.au / u.day
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
                     [hx, hy, hz])) / mu  - [self.x, self.y, self.z] * u.au/self.r
        self.e = norm([ex, ey, ez])

        # compute vector n
        ### hacky units
        nx, ny, nz = -hy.value, hx.value, np.zeros(len(hz))
        n = norm([nx, ny, nz])

        # compute true anomaly ν, the angle between e and r
        # true anomaly is undefined for circular orbits. Handle this more nicely.
        with np.errstate(invalid='ignore'):
            ν = np.arccos(dot([ex, ey, ez], [self.x, self.y, self.z]) / (self.e*self.r))
            ν[rrdot < 0] = 2 * np.pi * u.rad - ν[rrdot < 0]

        # compute inclination
        self.inc = Angle(np.arccos(hz/h), u.rad)

        # compute eccentric anomaly E
        E = 2 * np.arctan2(np.sqrt(1-self.e) * np.sin(ν/2), np.sqrt(1+self.e) * np.cos(ν/2))

        # compute ascending node
        with np.errstate(invalid='ignore'):
            node = np.arccos(nx/n)
            node[ny < 0] = 2 * np.pi - node[ny < 0]
            self.node = Angle(node, u.rad)

        # compute argument of periapsis, the angle between e and n
        with np.errstate(invalid='ignore'):
            arg = np.arccos(dot([nx, ny, nz], [ex, ey, ez]) / (n*self.e))
            arg[ez < 0] = 2 * np.pi * u.rad - arg[ez < 0]
            self.arg = Angle(arg, u.rad)

        # compute mean anomaly
        with np.errstate(invalid='ignore'):
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
        if self.__class__.frame == 'heliocentric':
            #t = ts.tai(jd=self.epoch.value + 37/86400)
            #t = ts.tdb(jd=self.epoch.tdb.jd)
            t = ts.tt(jd=self.epoch.tt.jd)
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
            with np.errstate(invalid='ignore'):
                lp = self.M < np.pi * u.rad
            self.t_peri.jd[lp] = self.epoch.value[lp] - self.M.value[lp] \
                                 / np.sqrt(mu_bary.value / self.a.value[lp]**3)
            self.t_peri.jd[~lp] = self.epoch.value[~lp] \
                                  + (2*np.pi - self.M.value[~lp]) \
                                  / np.sqrt(mu_bary.value / self.a.value[~lp]**3)
            self.t_peri = Time(self.t_peri, format='jd', scale='utc')
            self.__class__.frame = 'barycentric'

        return self


    def to_helio(self):
        '''
        Method to convert barycentric coordinates to heliocentric coordinates.
        '''
        if self.__class__.frame == 'barycentric':
            #t = ts.tai(jd=self.epoch.value + 37/86400)
            #t = ts.tdb(jd=self.epoch.tdb.jd)
            t = ts.tt(jd=self.epoch.tt.jd)
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
            with np.errstate(invalid='ignore'):
                lp = self.M < np.pi * u.rad
            self.t_peri.jd[lp] = self.epoch.value[lp] - self.M.value[lp] \
                                 / np.sqrt(mu_bary.value / self.a.value[lp]**3)
            self.t_peri.jd[~lp] = self.epoch.value[~lp] \
                                  + (2*np.pi - self.M.value[~lp]) \
                                  / np.sqrt(mu_bary.value / self.a.value[~lp]**3)
            self.t_peri = Time(self.t_peri, format='jd', scale='utc')
            self.__class__.frame = 'heliocentric'

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
        x, y, z = euler_rotation(arg, inc, node, o) #* u.au

        return x, y, z


    def calc_E(self, e, M):
        '''
        This method solves Kepler's Equation with the algorithm desscribed here:
        https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19950021346.pdf
        '''
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
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
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            E[E < 0] += 2 * np.pi
        return E

    def calc_t_peri(self):

        with np.errstate(invalid='ignore'):
            lp = self.M < np.pi * u.rad

        t_peri = np.zeros(len(self))
        t_peri[lp] = self.epoch.jd[lp] - self.M.value[lp] \
                          / np.sqrt(mu_bary.value / self.a.value[lp]**3)
        t_peri[~lp] = self.epoch.jd[~lp] \
                           + (2*np.pi - self.M.value[~lp]) \
                           / np.sqrt(mu_bary.value / self.a.value[~lp]**3)

        t_peri = Time(t_peri, format='jd', scale='utc')

        return t_peri
