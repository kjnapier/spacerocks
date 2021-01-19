import warnings
from healpy.pixelfunc import ang2pix

import numpy as np
from astropy import units as u
import pandas as pd

from .linalg3d import *
from .constants import *

from astropy.time import Time
from astropy.coordinates import SkyCoord
from skyfield.api import Topos, Loader
# Load in planets for ephemeride calculation.
load = Loader('./Skyfield-Data', expire=False, verbose=False)
ts = load.timescale()
planets = load('de423.bsp')
sun = planets['sun']
earth = planets['earth']

class Transformations:

    #def equa_to_ecl(self, x, y, z):

    #    R = np.array([[1, 0, 0],
    #                  [0, np.cos(epsilon), np.sin(epsilon)],
    #                  [0, -np.sin(epsilon), np.cos(epsilon)]])

    #    return R @ np.array([x, y, z])

    #def ecl_to_tel(self, x, y, z):

    #    #l =
    #    #b =

    #    R = np.array([[-np.sin(l), -np.cos(l) * np.sin(b), np.cos(l) * np.cos(b)],
    #                  [np.cos(l), -np.sin(l) * np.sin(b), np.sin(l) * np.cos(b)],
    #                  [0, np.cos(b), np.sin(b)]])

    #    return R @ np.array([x, y, z])

    #def xyz_to_abg(self):

    #    earth = planets['earth']

    #    if self.__class__.obscode != 500:
    #        earth += Topos(latitude_degrees=self.__class__.obslat,
    #                        longitude_degrees=self.__class__.obslon,
    #                        elevation_m=self.__class__.obselev) # topocentric calculation

    #    t = ts.tt(jd=self.epoch.tt.jd)
    #    x_earth, y_earth, z_earth = earth.at(t).position.au * u.au # earth ICRS position
    #    vx_earth, vy_earth, vz_earth = earth.at(t).velocity.au_per_d * u.au / u.day # earth ICRS position

    #    x, y, z = self.equa_to_ecl(self.x.value, self.y.value, self.z.value) * u.au
    #    vx, vy, vz = self.equa_to_ecl(self.vx.value, self.vy.value, self.vz.value) * u.au / u.day

    #    #x0 = self.x - x_earth
    #    #y0 = self.y - y_earth
    #    #z0 = self.z - z_earth
    #    #vx0 = self.vx - vx_earth
    #    #vy0 = self.vy - vy_earth
    #    #vz0 = self.vz - vz_earth

    #    x0 = x - x_earth
    #    y0 = y - y_earth
    #    z0 = z - z_earth
    #    vx0 = vx - vx_earth
    #    vy0 = vy - vy_earth
    #    vz0 = vz - vz_earth

    #    self.alpha = x0 / z0
    #    self.beta = y0 / z0
    #    self.gamma = 1 / z0
    #    self.dalpha = vx0 / z0
    #    self.dbeta = vy0 / z0
    #    self.dgamma = vz0 / z0

    #    return self

    #def abg_to_xyz(self):

    #    earth = planets['earth']

    #    if self.__class__.obscode != 500:
    #        earth += Topos(latitude_degrees=self.__class__.obslat,
    #                        longitude_degrees=self.__class__.obslon,
    #                        elevation_m=self.__class__.obselev) # topocentric calculation

    #    t = ts.tt(jd=self.epoch.tt.jd)
    #    x_earth, y_earth, z_earth = earth.at(t).position.au * u.au # earth ICRS position
    #    vx_earth, vy_earth, vz_earth = earth.at(t).velocity.au_per_d * u.au / u.day # earth ICRS position

    #    self.x = self.alpha / self.gamma + x_earth
    #    self.y = self.beta / self.gamma + y_earth
    #    self.z = 1 / self.gamma + z_earth
    #    self.vx = self.alpha_dot / self.gamma
    #    self.vy = self.beta_dot / self.gamma
    #    self.vz = self.gamma_dot / self.gamma

    #    return self

    def xyz_to_equa(self, rocks):
        '''
        Transform from barycentric Cartesian coordinates to equatorial
        coordinates. If you need very precise values, the method will use
        a topocentric correction to the Earth's position.
        See https://en.wikipedia.org/wiki/Horizontal_coordinate_system.
        This results in a significant slowdown to the code.
        '''
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

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

            earth_dis = norm([x_earth, y_earth, z_earth])

            x0, y0, z0 = rocks.x, rocks.y, rocks.z
            vx0, vy0, vz0 = rocks.vx, rocks.vy, rocks.vz
            for idx in range(5):
                # transfer ecliptic to ICRS and shift to Geocentric (topocentric)
                x = x0 - x_earth
                y = y0 * np.cos(epsilon) - z0 * np.sin(epsilon) - y_earth
                z = y0 * np.sin(epsilon) + z0 * np.cos(epsilon) - z_earth
                vx = vx0 - vx_earth
                vy = vy0 * np.cos(epsilon) - vz0 * np.sin(epsilon) - vy_earth
                vz = vy0 * np.sin(epsilon) + vz0 * np.cos(epsilon) - vz_earth
                delta = norm([x, y, z])
                ltt = delta / c
                M = rocks.M - ltt * (mu_bary / rocks.a**3)**0.5
                if idx < 4:
                    x0, y0, z0, vx0, vy0, vz0 = self.kep_to_xyz_temp(rocks.a, rocks.e, rocks.inc,
                                                                     rocks.arg, rocks.node, M)

            # Cartesian to spherical coordinate
            self.delta = norm([x, y, z])
            self.ltt = self.delta / c
            self.dec = Angle(np.arcsin(z / self.delta), u.rad)
            self.ra = Angle(np.arctan2(y, x), u.rad).wrap_at(2 * np.pi * u.rad)
            #self.phase_angle = Angle(np.arccos(-(earth_dis**2 - rocks.r**2 - self.delta**2)/(2 * rocks.r * self.delta)), u.rad)
            self.elong = Angle(np.arccos(-(rocks.r**2 - self.delta**2 - earth_dis**2)/(2 * self.delta * earth_dis)), u.rad)
            self.skycoord = SkyCoord(self.ra, self.dec, frame='icrs')
            self.dec_rate = (-z * (x * vx + y * vy) + (norm([x, y, np.zeros_like(x)])**2) * vz) \
                    / (norm([x, y, np.zeros_like(x)]) * self.delta**2) * u.rad
            self.ra_rate = - (y * vx - x * vy) / norm([x, y, np.zeros_like(z)])**2 * np.cos(self.dec) * u.rad

        return self

    def estimate_mag(self, rocks):
        '''
        Estimate the apparent magnitude of a TNO
        https://iopscience.iop.org/article/10.3847/1538-3881/ab18a9/pdf for light curves
        '''
        q = (rocks.r**2 + self.delta**2 - 1 * u.au**2)/(2 * rocks.r * self.delta)

        ## pyephem
        beta = np.zeros(len(q))

        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            beta = np.arccos(q)

        beta[np.where(q <= -1)[0]] = np.pi * u.rad
        beta[np.where(q >= 1)[0]] = 0 * u.rad

        Psi_1 = np.exp(-3.33 * np.tan(beta/2)**0.63)
        Psi_2 = np.exp(-1.87 * np.tan(beta/2)**1.22)
        mag = rocks.H + 5 * np.log10(rocks.r * self.delta / u.au**2)


        not_zero = np.where((Psi_1 != 0) | (Psi_2 != 0))[0]
        mag[not_zero] -= 2.5 * np.log10((1 - rocks.G[not_zero]) * Psi_1[not_zero] + rocks.G[not_zero] * Psi_2[not_zero])


        try:
            mag += rocks.delta_H * np.sin((rocks.epoch.jd - rocks.t0.value - self.ltt.value) * 2 * np.pi / rocks.rotation_period.value + rocks.phi0)
        except:
            pass

        return mag


    def radec_to_hpix(self, NSIDE):
        '''
        Convert (ra, dec) into healpix.
        '''
        return hp.pixelfunc.ang2pix(NSIDE, np.pi/2 - self.dec.radian, self.ra.radian, nest=True)


    def kep_to_xyz(self, mu):
        '''
        Transform from Keplerian to cartesian coordinates. There is no analytic
        solution to Kepler's Equation M = E - eSin[e] for the eccentric
        anomaly (E), but I use a close approximation.
        '''
        # compute eccentric anomaly E
        E = self.calc_E(self.e.value, self.M.value) * u.rad

        # compute true anomaly v
        ν = 2 * np.arctan2((1 + self.e)**0.5 * np.sin(E/2), (1 - self.e)**0.5 * np.cos(E/2))
        self.true = ν

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

        self.true = ν

        # compute inclination
        self.inc = Angle(np.arccos(hz/h), u.rad)

        # compute eccentric anomaly E
        E = 2 * np.arctan2(np.sqrt(1-self.e) * np.sin(ν/2), np.sqrt(1+self.e) * np.cos(ν/2))

        # compute ascending node
        node = np.arccos(nx/n) * u.rad
        node[ny < 0] = 2 * np.pi * u.rad - node[ny < 0]
        self.node = Angle(node, u.rad)

        # compute argument of periapsis, the angle between e and n
        arg = np.arccos(dot([nx, ny, nz], [ex, ey, ez]) / (n * self.e))
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
            self.t_peri.jd[lp] = self.epoch.jd[lp] - self.M.value[lp] \
                                 / np.sqrt(mu_bary.value / self.a.value[lp]**3)
            self.t_peri.jd[~lp] = self.epoch.jd[~lp] \
                                  + (2*np.pi - self.M.value[~lp]) \
                                  / np.sqrt(mu_bary.value / self.a.value[~lp]**3)

            #self.t_peri.jd[lp] = self.epoch.value[lp] - self.M.value[lp] \
            #                     / np.sqrt(mu_bary.value / self.a.value[lp]**3)
            #self.t_peri.jd[~lp] = self.epoch.value[~lp] \
            #                      + (2*np.pi - self.M.value[~lp]) \
            #                      / np.sqrt(mu_bary.value / self.a.value[~lp]**3)
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
            self.t_peri.jd[lp] = self.epoch.jd[lp] - self.M.value[lp] \
                                 / np.sqrt(mu_bary.value / self.a.value[lp]**3)
            self.t_peri.jd[~lp] = self.epoch.jd[~lp] \
                                  + (2*np.pi - self.M.value[~lp]) \
                                  / np.sqrt(mu_bary.value / self.a.value[~lp]**3)
            #self.t_peri.jd[lp] = self.epoch.value[lp] - self.M.value[lp] \
            #                     / np.sqrt(mu_bary.value / self.a.value[lp]**3)
            #self.t_peri.jd[~lp] = self.epoch.value[~lp] \
            #                      + (2*np.pi - self.M.value[~lp]) \
            #                      / np.sqrt(mu_bary.value / self.a.value[~lp]**3)
            self.t_peri = Time(self.t_peri, format='jd', scale='utc')
            self.__class__.frame = 'heliocentric'

        return self


    def kep_to_xyz_temp(self, a, e, inc, arg, node, M):
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
        ov = [(mu_bary * a)**0.5 / r * (-np.sin(E))/ u.rad,
              (mu_bary * a)**0.5 / r * ((1 - e**2)**0.5 * np.cos(E))/ u.rad,
              np.zeros(len(ν))/ u.rad]

        # Rotate o to the inertial frame
        x, y, z = euler_rotation(arg, inc, node, o) #* u.au
        vx, vy, vz = euler_rotation(arg, inc, node, ov) #* u.au / u.day

        return x, y, z, vx, vy, vz


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
            lp = self.M.rad < np.pi# * u.rad

        t_peri = np.zeros_like(self.M.rad)
        t_peri[lp] = self.epoch.jd[lp] - self.M.rad[lp] \
                          / np.sqrt(mu_bary.value / self.a.value[lp]**3)
        t_peri[~lp] = self.epoch.jd[~lp] \
                           + (2*np.pi - self.M.rad[~lp]) \
                           / np.sqrt(mu_bary.value / self.a.value[~lp]**3)

        t_peri = Time(t_peri, format='jd', scale='utc')

        return t_peri
