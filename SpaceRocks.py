import healpy as hp

import numpy as np
from scipy.optimize import newton
from skyfield.api import Topos, Loader
from astropy import units as u
from astropy.coordinates import Angle
from linalg3d import dot, norm, cross, euler_rotation
from astropy.constants import c
import matplotlib.pyplot as plt

# to do: implement obscode
# to do: implement plotting
# to do: implement polite kwarg rejection
# to do: implement input heliocentric

class SpaceRock:

    # Set precision for eccentric anomaly and radec calculations
    precise = False

    # Load in planets for ephemeride calculation.
    load = Loader('./Skyfield-Data', expire=False, verbose=False)
    planets = load('de423.bsp')
    sun = planets['sun']
    earth = planets['earth']
    earth += Topos('30.169 S', '70.804 W', elevation_m=2200)
    ts = load.timescale()

    # standard gravitational parameter, GM, M is the mass of sun + all planets
    μ_bary = (4 * np.pi**2) * 1.0013418393275118 * u.radian**2 * u.au**3 / u.year**2
    μ_helio = (4 * np.pi**2) * u.radian**2 * u.au**3 / u.year**2

    # obliquity of the Earth
    ε = Angle(23.43929111, u.degree).radian

    # speed of light
    c = c.to(u.au / u.year)

    def __init__(self, input_coordinates='keplerian', input_frame='barycentric',
                 input_angles='degrees', NSIDE=128, *args, **kwargs):

        # a, e, inc, node, arg, M0, x, y, z, vx, vy, vz, tau, epoch,

        # Case-insensitive keyword arguments.
        kwargs = {key.lower(): data for key, data in kwargs.items()}
        input_coordinates = input_coordinates.lower()
        input_frame = input_frame.lower()
        input_angles = input_angles.lower()

        # scalar input -> arrays
        #for idx, key in enumerate([*kwargs]):
        #    if np.isscalar(kwargs.get(key)):
        #        kwargs[key] = np.array([kwargs.get(key)])
        if np.isscalar(NSIDE):
            NSIDE = np.array([NSIDE])


        if input_angles == 'degrees':
            angle_unit = u.degree
        elif input_angles == 'radians':
            angle_unit = u.radian
        else:
            raise ValueError('The input_angles argument must be a string \
                              that reads either degrees or radians.')

        if input_coordinates == 'keplerian':

            self.frame = input_frame
            self.hash = kwargs.get('hash')
            self.a = kwargs.get('a') * u.au
            self.e = kwargs.get('e') * u.dimensionless_unscaled
            self.inc = Angle(kwargs.get('inc'), angle_unit).to(u.rad)
            self.node = Angle(kwargs.get('node'), angle_unit).to(u.rad)
            self.arg = Angle(kwargs.get('arg'), angle_unit).to(u.rad)

            self.epoch = kwargs.get('epoch') * u.day # time at perihelion
            self.tau = kwargs.get('tau') * u.day # obs_date
            self.M = np.sqrt(μ / self.a**3) * (self.tau - self.epoch)

            self.x, self.y, self.z, self.vx, self.vy, self.vz = self.kep_to_xyz(μ_bary)

            # If the input frame is heliocentric, convert to bary to calculate
            if input_frame == 'heliocentric':
                self.to_bary()
                self.frame = 'barycentric'

                self.r = norm([self.x, self.y, self.z])
                self.ra, self.dec, self.delta, self.ltt, self.phase_angle, self.elong = self.xyz_to_equa()

                for value in hpix:
                    setattr(blah, 'HPIX_{}'.format(value), self.func(value))

                self.to_helio()
                self.r = norm([self.x, self.y, self.z])
                self.frame = 'heliocentric'

            self.r = norm([self.x, self.y, self.z])
            self.ra, self.dec, self.delta, self.ltt, self.phase_angle, self.elong = self.xyz_to_equa()



        elif input_coordinates == 'cartesian':

            self.frame = input_frame
            self.hash = kwargs.get('hash')
            self.x = kwargs.get('x') * u.au
            self.y = kwargs.get('y') * u.au
            self.z = kwargs.get('z') * u.au
            self.vx = kwargs.get('vx') * (u.au / u.year)
            self.vy = kwargs.get('vy') * (u.au / u.year)
            self.vz = kwargs.get('vz') * (u.au / u.year)
            self.r = norm([self.x, self.y, self.z])
            self.a, self.e, self.inc, self.arg, self.node, self.M = self.xyz_to_kep(μ_bary)
            self.ra, self.dec, self.delta, self.r, self.ltt, self.phase_angle, self.elong = self.xyz_to_equa()
            for value in hpix:
                setattr(blah, 'HPIX_{}'.format(value), self.func(value))

        self.varpi = (self.arg + self.node).wrap_at(2 * np.pi * u.rad)

    if precise == True:

        def cal_E(self):
            # compute eccentric anomaly E
            f = lambda calc_E, M, e: calc_E - e * np.sin(calc_E) * u.rad - M
            E0 = self.M
            calc_E = newton(f, E0, args=(self.M, self.e))
            return calc_E

        def xyz_to_equa(self):
            t = ts.tt(jd=self.epoch.value)
            x_earth, y_earth, z_earth = earth.at(t).position.au * u.au # earth ICRS position
            earth_dis = norm([x_earth, y_earth, z_earth])
            x0, y0, z0 = self.x, self.y, self.z
            for idx in range(3):
                # transfer ecliptic to ICRS and shift to Geocentric (topocentric)
                x = x0 - x_earth
                y = y0 * np.cos(ε) - z0 * np.sin(ε) - y_earth
                z = y0 * np.sin(ε) + z0 * np.cos(ε) - z_earth
                delta = norm([x, y, z])
                ltt = delta / c
                M = self.M - ltt * (μ_bary / self.a**3)**0.5
                if idx < max(range(3)):
                    x0, y0, z0, _, _, _ = self.kep_to_xyz(μ_bary)
            # Cartesian to spherical coordinate
            dec = Angle(np.arcsin(z / norm([x, y, z])), u.rad)
            ra = Angle(np.arctan2(y, x), u.rad).wrap_at(2 * np.pi * u.rad)
            phase_angle = Angle(np.arccos(-(earth_dis**2 - self.r**2 - delta**2)/(2 * self.r* delta)), u.rad)
            elong = Angle(np.arccos(-(self.r**2 - delta**2 - earth_dis**2)/(2 * delta * earth_dis)), u.rad)
            return ra, dec, delta, ltt, phase_angle, elong

    elif precise == False:

        def cal_E(self):
            calc_E = self.M
            for kk in range(100):
                calc_E = self.M + self.e * np.sin(calc_E) * u.rad
            return calc_E

        def xyz_to_equa(self):
            # transfer ecliptic to ICRS and shift to Geocentric (topocentric)
            dx = self.x - x_earth
            dy = self.y * np.cos(ε) - self.z * np.sin(ε) - y_earth
            dz = self.y * np.sin(ε) + self.z * np.cos(ε) - z_earth
            delta = np.sqrt(dx**2 + dy**2+ dz**2)

            # Cartesian to spherical coordinate
            dec = np.arcsin(dz/delta)
            ra = np.arctan2(dy, dx) % (2 * pi)
            mag2 = H + 2.5 * log10(self.r**2 * delta**2)

            c = (self.r**2 + delta**2 - 1)/(2 * self.r * delta)

            beta = np.zeros(len(c))
            beta[np.where(c <= 1)[0]] = np.pi
            beta[np.where(c >= 1)[0]] = 0

            tb2 = np.tan(beta/2.0);

            psi_t = tb2**0.63
            Psi_1 = np.exp(-3.33 * psi_t)

            psi_t = tb2**1.22
            Psi_2 = np.exp(-1.87 * psi_t)
            mag = H + 5 * np.log10(r * delta)
            not_zero = np.where((Psi_1 != 0) | (Psi_2 != 0))[0]
            mag[not_zero] -= 2.5 * np.log10((1 - G) * Psi_1[not_zero] + G * Psi_2[not_zero])

            return ra, dec, delta, r, mag, mag2

    else:
        raise ValueError('The parameter precise must be set to either True or False.')

    def radec_to_hpix(self, NSIDE):
        # Convert (ra, dec) into healpix
        return hp.pixelfunc.ang2pix(NSIDE, np.pi/2 - self.dec, self.ra, nest)

    def kep_to_xyz(self, μ):
        # compute eccentric anomaly E
        E = self.cal_E()

        # compute true anomaly v
        ν = 2 * np.arctan2((1 + self.e)**0.5*np.sin(E/2.), (1 - self.e)**0.5*np.cos(E/2.))

        # compute the distance to the central body r
        r = self.a * (1 - self.e*np.cos(E))

        # obtain the position o and velocity ov vector
        o = [r * np.cos(ν), r * np.sin(ν), np.zeros(len(ν))]
        ov = [(μ * self.a)**0.5 / r * (-np.sin(E)),
              (μ * self.a)**0.5 / r * ((1-self.e**2)**0.5 * np.cos(E)),
              np.zeros(len(ν))]

        # Rotate o and ov to the inertial frame
        x, y, z = euler_rotation(self.arg, self.inc, self.node, o) * u.au
        vx, vy, vz = euler_rotation(self.arg, self.inc, self.node, ov) * u.au / u.year

        return x, y, z, vx, vy, vz


    def xyz_to_kep(self, μ):
        ## there are a lot of hacky units in here because astropy doesn't behave as expected.
        # compute the barycentric distance r
        r = norm([self.x, self.y, self.z])
        rrdot = dot([self.x, self.y, self.z], [self.vx, self.vy, self.vz])

        # compute the specific angular momentum h
        hx, hy, hz = cross([self.x, self.y, self.z], [self.vx, self.vy, self.vz])
        h = norm([hx, hy, hz])

        # compute eccentricity vector
        ### hacky units
        ex, ey, ez = u.au**3 / u.yr**2 * np.array(cross([self.vx, self.vy, self.vz], [hx, hy, hz])) *u.rad**2 / μ  - [self.x, self.y, self.z]*u.au/r
        e = norm([ex, ey, ez])

        # compute vector n
        ### hacky units
        nx, ny, nz = -hy.value, hx.value, np.zeros(len(hz))
        n = norm([nx, ny, nz])

        # compute true anomaly ν, the angle between e and r
        ν = np.arccos(dot([ex, ey, ez], [self.x, self.y, self.z]) / (e*r))
        ν[rrdot < 0] = 2 * np.pi * u.rad - ν[rrdot < 0]

        # compute inclination
        inc = np.arccos(hz/h)

        # compute eccentric anomaly E
        E = 2 * np.arctan2(np.sqrt(1-e) * np.sin(ν/2), np.sqrt(1+e) * np.cos(ν/2))

        # compute ascending node
        Ω = np.arccos(nx/n)
        Ω[ny < 0] = 2 * np.pi * u.rad - Ω[ny < 0]

        # compute argument of periapsis, the angle between e and n
        ω = np.arccos(dot([nx, ny, nz], [ex, ey, ez]) / (n*e))
        ω[ez < 0] = 2 * np.pi * u.rad - ω[ez < 0]

        # compute mean anomaly
        M = E - e * np.sin(E) * u.rad
        M[M < 0] += 2 * np.pi * u.rad

        # compute a
        a = 1 / (2 / r - norm([self.vx, self.vy, self.vz])**2 / μ * u.rad**2)

        return a, e, inc, ω, Ω, M

    def to_bary(self):
        if self.frame == 'heliocentric':
            t = ts.tt(jd=self.epoch.value) #37 leap seconds
            x_sun, y_sun, z_sun = sun.at(t).position.au * u.au
            vx_sun, vy_sun, vz_sun = (sun.at(t).velocity.au_per_d * u.au / u.day).to(u.au / u.year)
            # calculate the barycentric xyz postion
            self.x = self.x + x_sun
            self.y = self.y + y_sun * np.cos(ε) - z_sun * np.sin(-ε)
            self.z = self.z + y_sun * np.sin(-ε) + z_sun * np.cos(ε)
            self.vx = self.vx + vx_sun
            self.vy = self.vy + vy_sun * np.cos(ε) - vz_sun * np.sin(-ε)
            self.vz = self.vz + vy_sun * np.sin(-ε) + vz_sun * np.cos(ε)

            # calculate barycentric keplerian elements
            self.a, self.e, self.i, self.arg, self.node, self.M = self.xyz_to_kep(μ_bary)

            self.frame = 'barycentric'

        return


    def to_helio(self):

        if self.frame == 'barycentric':
            t = ts.tt(jd=self.epoch.value)
            x_sun, y_sun, z_sun = sun.at(t).position.au * u.au
            vx_sun, vy_sun, vz_sun = (sun.at(t).velocity.au_per_d * u.au / u.day).to(u.au / u.year)
            # calculate the heliocentric xyz postion
            self.x = self.x - x_sun
            self.y = self.y - y_sun * np.cos(ε) + z_sun * np.sin(-ε)
            self.z = self.z - y_sun * np.sin(-ε) - z_sun * np.cos(ε)
            self.vx = self.vx - vx_sun
            self.vy = self.vy - vy_sun * np.cos(ε) + vz_sun * np.sin(-ε)
            self.vz = self.vz - vy_sun * np.sin(-ε) - vz_sun * np.cos(ε)
            self.r = norm([self.x, self.y, self.z])
            # calculate heliocentric keplerian elements
            self.a, self.e, self.i, self.arg, self.node, self.M = self.xyz_to_kep(μ_helio)

            self.frame = 'heliocentric'

        return

    def plot_radec(self, color):
        fig = plt.figure(figsize=(12, 8))
        ax = fig.add_subplot(111, projection="mollweide")

        plt.scatter(self.ra, self.dec, color=color)
