import healpy as hp

import numpy as np
from scipy.optimize import newton
from skyfield.api import Topos, Loader
from astropy import units as u
from astropy.coordinates import Angle
from aeiderivs import aeiderivs
from linalg3d import dot, norm, cross, euler_rotation

load = Loader('./Skyfield-Data', expire=False, verbose=False)
planets = load('de423.bsp')
sun = planets['sun']
earth = planets['earth']
earth += Topos('30.169 S', '70.804 W', elevation_m=2200)
ts = load.timescale()

# set global constants
# standard gravitational parameter, GM, M is the mass of sun + all planets
μ_bary = (4 * np.pi**2) * 1.0013418393275118 * u.radian**2 * u.au**3 / u.year**2
μ_helio = (4 * np.pi**2) * u.radian**2 * u.au**3 / u.year**2

# obliquity of the Earth
ε = Angle(23.43929111, u.degree).radian

class SpaceRock:

    def __init__(self, coordinates='keplerian', frame='barycentric', precise=False,
                 input_angles='degrees', NSIDE=128, *args, **kwargs):

        if frame.lower() == 'barycentric':
            μ = μ_bary
        elif frame.lower() == 'heliocentric':
            μ = μ_helio
        else:
            raise ValueError('The frame must bea string that reads either \
                              barycentric or heliocentric.')

        # a, e, inc, node, omega, M0, x, y, z, vx, vy, vz, tau, epoch,

        # Case-insensitive keyword arguments.
        kwargs =  {key.lower(): data for key, data in kwargs.items()}

        # scalar input -> arrays
        for idx, key in enumerate([*kwargs]):
            if np.isscalar(kwargs.get(key)):
                kwargs[key] = np.array([kwargs.get(key)])

        if input_angles == 'degrees':
            angle_unit = u.degree
        elif input_angles == 'radians':
            angle_unit = u.radian
        else:
            raise ValueError('The input_angles argument must be a string \
                              that reads either degrees or radians.')

        if coordinates == 'keplerian':

            self.a = kwargs.get('a') * u.au
            self.e = kwargs.get('e') * u.dimensionless_unscaled
            self.inc = Angle(kwargs.get('inc'), angle_unit).to(u.rad)
            self.node = Angle(kwargs.get('node'), angle_unit).to(u.rad)
            self.omega = Angle(kwargs.get('omega'), angle_unit).to(u.rad)

            # self.M0 = Angle(kwargs.get('M0'), angle_unit).radian
            self.epoch = kwargs.get('epoch') * u.day # time at perihelion
            self.tau = kwargs.get('tau') * u.day # obs_date

            self.M = np.sqrt(μ / self.a**3) * (self.tau - self.epoch)
            # self.M = np.sqrt(μ / self.a**3) * (self.tau - self.epoch) + self.M0
            self.x, self.y, self.z, self.vx, self.vy, self.vz = self.kep_to_xyz()
            #self.ra, self.dec, self.delta, self.r, self.ltt, self.phase_angle, self.elong = self.xyz_to_equa()

        elif coordinates == 'cartesian':

            self.x = kwargs.get('x') * u.au
            self.y = kwargs.get('y') * u.au
            self.z = kwargs.get('z') * u.au
            self.vx = kwargs.get('vx') * (u.au / u.year)
            self.vy = kwargs.get('vy') * (u.au / u.year)
            self.vz = kwargs.get('vz') * (u.au / u.year)
            self.r = norm([self.x, self.y, self.z])
            self.a, self.e, self.inc, self.omega, self.node, self.M = self.xyz_to_kep()
            #self.ra, self.dec, self.delta, self.r, self.ltt, self.phase_angle, self.elong = self.xyz_to_equa()

        elif coordinates == 'equatorial':

            self.ra = Angle(kwargs.get('ra'), angle_unit).to(u.rad)
            self.dec = Angle(kwargs.get('dec'), angle_unit).to(u.rad)
            # self.orbitalelements = Orbitfit()

        else:
            raise ValueError('The specified coordinate system must be keplerian, \
                              cartesian, or equatorial.')
            sys.exit(1)

        self.varpi = (self.omega + self.node).wrap_at(2 * np.pi * u.rad)
        #self.ra, self.dec, self.delta, self.r, self.ltt, self.phase_angle, self.elong = self.xyz_to_equa()

    if precise == False:

        def cal_E(self):
            calc_E = self.M
            for kk in range(100):
                calc_E = self.M + self.e * np.sin(calc_E) * u.rad
            return calc_E

        def calcradec(self):
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

        def cal_E(self):
            # compute eccentric anomaly E
            f = lambda calc_E, M, e: calc_E - e * np.sin(calc_E) * u.rad - M
            E0 = self.M
            calc_E = newton(f, E0, args=(self.M, self.e))
            return calc_E


    # Convert (ra, dec) into healpix
    def radec_to_hpix(ra, dec, NSIDE=64, nest=True):
        return hp.pixelfunc.ang2pix(NSIDE, np.pi/2 - dec, ra, nest)

    def kep_to_xyz(self):
        # compute eccentric anomaly E
        E = self.cal_E()

        # compute true anomaly v
        v = 2 * np.arctan2((1 + self.e)**0.5*np.sin(E/2.), (1 - self.e)**0.5*np.cos(E/2.))

        # compute the distance to the central body r
        r = self.a * (1 - self.e*np.cos(E))

        # obtain the position o and velocity ov vector
        o = [r * np.cos(v), r * np.sin(v), 0]
        ov = [(μ * self.a)**0.5 / r * (-np.sin(E)),
              (μ * self.a)**0.5 / r * ((1-self.e**2)**0.5 * np.cos(E)),
              0]

        # Rotate o and ov to the inertial frame
        x, y, z = euler_rotation(self.omega, self.inc, self.node, [ox, oy, oz])
        vx, vy, vz = euler_rotation(self.omega, self.inc, self.node, [ovx, ovy, ovz])

        return x, y, z, vx/u.rad, vy/u.rad, vz/u.rad


    def xyz_to_kep(self):
        # compute the barycentric distance r
        r = norm([self.x, self.y, self.z])
        rrdot = dot([self.x, self.y, self.z], [self.vx, self.vy, self.vz])

        # compute the specific angular momentum h
        hx, hy, hz = cross([self.x, self.y, self.z], [self.vx, self.vy, self.vz])
        h = norm([hx, hy, hz])

        # compute eccentricity vector
        ex, ey, ez = cross([vx, vy, vz], [hx, hy, hz]) * u.rad**2 / μ - [x, y, z]/r
        e = norm([ex, ey, ez])

        # compute vector n
        nx, ny, nz = -hy, hx, 0
        n = norm([nx, ny, nz])

        # compute true anomaly ν, the angle between e and r

        ν = np.arccos(dot([ex, ey, ez], [self.x, self.y, self.z]) / (e*r))
        ν[rrdot < 0] = 2 * np.pi - ν[rrdot < 0]

        # compute inclination
        inc = np.arccos(hz/h)

        # compute eccentric anomaly E
        E = 2 * np.arctan2(np.sqrt(1-e) * np.sin(ν/2), np.sqrt(1+e) * np.cos(ν/2))

        # compute ascending node
        Ω = np.arccos(nx/n)
        Ω[ny < 0] = 2 * np.pi - Ω[ny < 0]

        # compute argument of periapsis, the angle between e and n

        ω = np.arccos(dot([nx, ny, nz], [ex, ey, ez]) / (n*e))
        ω[ez < 0] = 2 * np.pi - ω[ez < 0]

        # compute mean anomaly
        M = E - e * np.sin(E)
        M[M < 0] += 2 * np.pi

        # compute a
        a = 1 / (2 / r - norm([self.vx, self.vy, self.vz])**2 / μ)

        return a, e, inc, ω, Ω, M

    def helio_to_bary(self, a0, e0, i0, arg0, node0, M0, epoch0):
        # This is heliocentric xyz postion
        X0, Y0, Z0, VX0, VY0, VZ0 =  self.kep_to_xyz(a0, e0, i0, arg0, node0, M0, self.u_helio)
        # extract barycentric postion of Sun
        sun = planets['Sun']
        ts = load.timescale()
        t = ts.tt(jd=epoch0 ) #37 leap seconds
        x_sun, y_sun, z_sun = sun.at(t).position.au
        vx_sun, vy_sun, vz_sun = sun.at(t).velocity.au_per_d
        # now we have barycentric xyz postion
        X = X0 + x_sun
        Y = Y0 + y_sun * np.cos(-ε) - z_sun * np.sin(-ε)
        Z = Z0 + y_sun * np.sin(-ε) + z_sun * np.cos(-ε)
        VX = VX0 + vx_sun
        VY = VY0 + vy_sun * np.cos(-ε) - vz_sun * np.sin(-ε)
        VZ = VZ0 + vy_sun * np.sin(-ε) + vz_sun * np.cos(-ε)
        # transfer back to keplerian elements
        a, e, i, arg, node, M = self.xyz_to_kep(X, Y, Z, VX, VY, VZ, self.u_bary)
        return a, e, i, arg, node, M

    def bary_to_helio(self, a0, e0, i0, arg0, node0, M0, epoch0):
        # This is barycentric xyz postion
        X0, Y0, Z0, VX0, VY0, VZ0 =  self.kep_to_xyz(a0, e0, i0, arg0, node0, M0, self.u_bary)
        # extract barycentric postion of Sun
        sun = planets['Sun']
        ts = load.timescale()
        t = ts.tt(jd=epoch0) #37 leap seconds
        x_sun, y_sun, z_sun = sun.at(t).position.au
        vx_sun, vy_sun, vz_sun = sun.at(t).velocity.au_per_d
        # now we have barycentric xyz postion
        X = X0 - x_sun
        Y = Y0 - y_sun * np.cos(-self.epsilon) + z_sun * np.sin(-self.epsilon)
        Z = Z0 - y_sun * np.sin(-self.epsilon) - z_sun * np.cos(-self.epsilon)
        VX = VX0 - vx_sun
        VY = VY0 - vy_sun * np.cos(-self.epsilon) + vz_sun * np.sin(-self.epsilon)
        VZ = VZ0 - vy_sun * np.sin(-self.epsilon) - vz_sun * np.cos(-self.epsilon)
        # transfer back to keplerian elements
        a, e, i, arg, node, M = self.xyz_to_kep(X, Y, Z, VX, VY, VZ, self.u_helio)
        return a, e, i, arg, node, M

    def xyz_to_equa(self, X0, Y0, Z0, epoch):
        c = 173.1446323547978
        r = (X0**2 + Y0**2 + Z0**2)**0.5
        earth = planets['earth']
        earth = earth + Topos('30.169 S', '70.804 W', elevation_m=2200) #turn off the topocentric calculation should run much faster
        ts = load.timescale()
        t = ts.tt(jd=epoch)
        x_earth, y_earth, z_earth = earth.at(t).position.au # earth ICRS position
        earth_dis = (x_earth**2 + y_earth**2 + z_earth**2)**0.5
        for i in range(3):
            # transfer ecliptic to ICRS and shift to Geocentric (topocentric)
            X = X0 - x_earth
            Y = Y0 * np.cos(self.epsilon) - Z0 * np.sin(self.epsilon) - y_earth
            Z = Y0 * np.sin(self.epsilon) + Z0 * np.cos(self.epsilon) - z_earth
            delta = (X**2 + Y**2+ Z**2)**0.5
            ltt = delta / c
            M = (self.u_bary/self.a**3)**0.5 * (-ltt) + self.M
            X0, Y0, Z0, VX0, VY0, VZ0 = self.kep_to_xyz(self.a, self.e, self.i, self.arg, self.node, M, self.u_bary)
        # Cartesian to spherical coordinate
        dec = np.arcsin(Z/(X**2+Y**2+Z**2)**0.5)
        ra = np.arctan2(Y, X) % (2*np.pi)
        phase_angle = np.arccos(-(earth_dis**2-r**2-delta**2)/(2*r*delta))
        elong = np.arccos(-(r**2-delta**2-earth_dis**2)/(2*delta*earth_dis))
        return ra, dec, delta, r, ltt, phase_angle, elong
