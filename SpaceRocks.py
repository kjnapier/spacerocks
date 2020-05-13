import healpy as hp

import numpy as np
from scipy.optimize import newton
from skyfield.api import Topos, Loader
from astropy import units as u
from astropy.coordinates import Angle
#from aeiderivs import aeiderivs

load = Loader('./Skyfield-Data', expire=False, verbose=False)
planets = load('de423.bsp')
earth = planets['earth']
earth += Topos('30.169 S', '70.804 W', elevation_m=2200)
ts = load.timescale()

# standard gravitational parameter, GM, M is the mass of sun + all planets
μ = (4 * np.pi**2) * 1.0013418393275118 * u.radian**2 * u.au**3 / u.year**2

# obliquity of the Earth
ε = Angle(23.43929111, u.degree).radian
precise=False
class SpaceRock:

    def __init__(self, coordinates='keplerian', precise=False, input_angles='degrees', *args, **kwargs):

        # a, e, inc, node, omega, M0, x, y, z, vx, vy, vz, tau, epoch,

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

        # self.M0 = Angle(kwargs.get('M0'), angle_unit).radian
        self.epoch = kwargs.get('epoch') * u.day # time at perihelion
        self.tau = kwargs.get('tau') * u.day # obs_date

        if coordinates == 'keplerian':

            self.a = kwargs.get('a') * u.au
            self.e = kwargs.get('e') * u.dimensionless_unscaled
            self.inc = Angle(kwargs.get('inc'), angle_unit).to(u.rad)
            self.node = Angle(kwargs.get('node'), angle_unit).to(u.rad)
            self.omega = Angle(kwargs.get('omega'), angle_unit).to(u.rad)

            self.M = np.sqrt(μ / self.a**3) * (self.tau - self.epoch)
            # self.M = np.sqrt(μ / self.a**3) * (self.tau - self.epoch) + self.M0
            self.X, self.Y, self.Z, self.VX, self.VY, self.VZ = self.kep_to_xyz()

        elif coordinates == 'cartesian':
            self.X = kwargs.get('X') * u.au
            self.Y = kwargs.get('Y') * u.au
            self.Z = kwargs.get('Z') * u.au
            self.VX = kwargs.get('VX') * u.au / u.year
            self.VY = kwargs.get('VY') * u.au / u.year
            self.VZ = kwargs.get('VZ') * u.au / u.year
            self.r = np.sqrt(self.X**2 + self.Y**2 + self.Z**2)
            self.a, self.e, self.inc, self.omega, self.node, self.M = self.xyz_to_kep()

        else:
            self.ra = Angle(kwargs.get('ra'), angle_unit).to(u.rad)
            self.dec = Angle(kwargs.get('dec'), angle_unit).to(u.rad)

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
            dX = self.X - x_earth
            dY = self.Y * np.cos(ε) - self.Z * np.sin(ε) - y_earth
            dZ = self.Y * np.sin(ε) + self.Z * np.cos(ε) - z_earth
            delta = np.sqrt(dX**2 + dY**2+ dZ**2)

            # Cartesian to spherical coordinate
            dec = np.arcsin(dZ/delta)
            ra = np.arctan2(dY, dX) % (2 * pi)
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

    def Rotation(self):

        R1 = np.array([[np.cos(self.omega), -np.sin(self.omega), 0],
                         [np.sin(self.omega), np.cos(self.omega), 0],
                         [0, 0, 1]])

        R2 = np.array([[1, 0, 0],
                         [0, np.cos(self.inc), -np.sin(self.inc)],
                         [0, np.sin(self.inc), np.cos(self.inc)]])

        R3 = np.array([[np.cos(self.node), -np.sin(self.node), 0],
                         [np.sin(self.node), np.cos(self.node), 0],
                         [0, 0, 1]])

        return np.linalg.multi_dot([R3, R2, R1])

    def kep_to_xyz(self):
        # compute eccentric anomaly E
        E = self.cal_E()
        # compute true anomaly v
        v = 2 * np.arctan2((1 + self.e)**0.5*np.sin(E/2.), (1 - self.e)**0.5*np.cos(E/2.))
        # compute the distance to the central body r
        r = self.a * (1 - self.e*np.cos(E))
        # obtain the position o and velocity ov vector
        ox = r * np.cos(v)
        oy = r * np.sin(v)
        oz = 0
        ovx = (μ * self.a)**0.5 / r * (-np.sin(E))
        ovy = (μ * self.a)**0.5 / r * ((1-self.e**2)**0.5 * np.cos(E))
        ovz = 0
        # Transform o and ov to the inertial frame
        X = ox * (np.cos(self.omega)*np.cos(self.node) - np.sin(self.omega)*np.sin(self.node)*np.cos(self.inc)) - oy * (np.sin(self.omega)*np.cos(self.node) + np.cos(self.omega)*np.sin(self.node)*np.cos(self.inc))
        Y = ox * (np.cos(self.omega)*np.sin(self.node) + np.sin(self.omega)*np.cos(self.node)*np.cos(self.inc)) + oy * (np.cos(self.omega)*np.cos(self.node)*np.cos(self.inc) - np.sin(self.omega)*np.sin(self.node))
        Z = ox * (np.sin(self.omega)*np.sin(self.inc)) + oy * (np.cos(self.omega)*np.sin(self.inc))
        VX = ovx * (np.cos(self.omega)*np.cos(self.node) - np.sin(self.omega)*np.sin(self.node)*np.cos(self.inc)) - ovy * (np.sin(self.omega)*np.cos(self.node) + np.cos(self.omega)*np.sin(self.node)*np.cos(self.inc))
        VY = ovx * (np.cos(self.omega)*np.sin(self.node) + np.sin(self.omega)*np.cos(self.node)*np.cos(self.inc)) + ovy * (np.cos(self.omega)*np.cos(self.node)*np.cos(self.inc) - np.sin(self.omega)*np.sin(self.node)) 
        VZ = ovx * (np.sin(self.omega)*np.sin(self.inc)) + ovy * (np.cos(self.omega)*np.sin(self.inc))
        return X, Y, Z, VX/u.rad, VY/u.rad, VZ/u.rad

    #def kep_to_xyz(self).
    #    # compute eccentric anomaly E
    #    #E = np.array(list(map(self.cal_E, self.e, self.M)))
    #    E = self.cal_E()

    #    # compute true anomaly v
    #    ν = 2 * np.arctan2(np.sqrt(1 + self.e) * np.sin(E/2), np.sqrt(1 - self.e) * np.cos(E/2))

    #    # compute the distance to the central body r
    #    r = self.a * (1 - self.e * np.cos(E))

    #    # obtain the position o and velocity ov vector
    #    ox, oy, oz = r * np.array([np.cos(ν), np.sin(ν), 0])
    #    ovx, ovy, ovz = np.sqrt(μ * self.a) / r * np.array([-np.sin(E), np.sqrt(1 - self.e**2) * np.cos(E), 0])

    #    X, Y, Z = np.dot([ox, oy, oz], self.Rotation())
    #    VX, VY, VZ = np.dot([ovx, ovy, ovz], self.Rotation())

    #    return X, Y, Z, VX, VY, VZ


    def xyz_to_kep(self):
        rrdot = self.X * self.VX + self.Y * self.VY + self.Z * self.VZ
        # compute the specific angular momentum h
        hx = self.Y * self.VZ - self.Z * self.VY
        hy = self.Z * self.VX - self.X * self.VZ
        hz = self.X * self.VY - self.Y * self.VX
        h = (hx**2 + hy**2 + hz**2)**0.5
        # compute eccentricity vector
        ex = (self.VY * hz - self.VZ * hy)/μ - self.X/self.r
        ey = (self.VZ * hx - self.VX * hz)/μ - self.Y/self.r
        ez = (self.VX * hy - self.VY * hx)/μ - self.Z/self.r
        e = (ex**2 + ey**2 + ez**2)**0.5
        # compute vector n
        nx = -hy
        ny = hx
        nz = 0
        n = np.sqrt(nx**2 + ny**2)
        # compute true anomaly v, the angle between e and r
        ν = np.arccos((ex * self.X + ey * self.Y + ez * self.Z) / (e*self.r))
        ν[rrdot < 0] = 2 * np.pi - v[rrdot < 0]
        # compute inclination
        i = np.arccos(hz / h)
        # compute eccentric anomaly E
        E = 2 * np.sqrt(np.arctan2((1 - e)) * np.sin(ν / 2), np.sqrt(1 + e) * np.cos(ν / 2))
        # compute ascending node
        Ω = np.arccos(nx / n)
        Ω[ny < 0] = 2 * np.pi - Ω[ny < 0]
        # compute argument of periapsis, the angle between e and n
        ω = np.arccos((nx * ex + ny * ey + nz *ez) / (n * e))
        ω[ez < 0] = 2 * np.pi - ω[ez < 0]
        # compute mean anomaly
        M = E - e * np.sin(E)
        M[M < 0] += 2 * np.pi
        # compute a
        a = 1/(2/self.r - (self.VX**2 + self.VY**2 + self.VZ**2)/μ)
        return a, e, Angle(i, u.radian), Angle(ω, u.radian), \
               Angle(Ω, u.radian), Angle(M, u.radian)



class Date:

    def __init__(self, day, month, year, UT):

        #for idx, key in enumerate([args]):
        #    if np.isscalar(kwargs.get(key)):
        #        kwargs[key] = np.array([kwargs.get(key)])

        self.day = day
        self.month = month
        self.year = year
        self.UT = UT # Universal Time

    def JD(self):

        m = np.zeros(len(self.month))
        y = np.zeros(len(self.year))
        B = np.zeros(len(self.year))

        aa = self.month > 2
        y[aa] = self.year[aa]
        m[aa] = self.month[aa]
        y[~aa] = self.year[~aa] - 1
        m[~aa] = self.month[~aa] + 12

        bb = self.year > 1582
        B[bb] = y[bb] // 400 - y[bb] // 100

        cc = self.year < 1582
        B[cc] = -2

        dd = self.year == 1582
        ee = self.month < 10
        B[dd & ee] = -2

        ff = self.month > 10
        B[dd & ff] = y[dd & ff] // 400 - y[dd & ff] // 100

        gg = self.month == 10
        hh = self.day <= 4
        B[dd & gg & hh] = -2

        ii = self.day >= 15
        B[dd & gg & ii] = y[dd & gg & ii] // 400 - y[dd & gg & ii] // 100

        return np.floor(365.25 * y) + np.floor(30.6001 * (m + 1)) + \
               1720996.5 + B + self.day + self.UT / 24


    def MJD(self):
        return self.JD() - 2400000.5
