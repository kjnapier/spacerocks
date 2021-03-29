from numpy import pi, sqrt, cbrt, sin, cos, tan, exp, log10, zeros_like, arccos, arctan2, array
import pandas as pd

from astropy import units as u

from .linalg3d import norm, dot, cross, euler_rotation
from astropy.time import Time

class OrbitFuncs:


    '''Private Methods'''

    def __calc_E(self, e, M):
        '''
        Calculate eccentric anomaly from eccentricity and mean anomaly

        This method solves Kepler's Equation with the algorithm desscribed here:
        https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19950021346.pdf
        '''

        M[M > pi] -= 2 * pi
        α = (3 * pi**2 + 1.6 * (pi**2 - pi * abs(M))/(1 + e))/(pi**2 - 6)
        d = 3 * (1 - e) + α * e
        q = 2 * α * d * (1 - e) - M**2
        r = 3 * α * d * (d - 1 + e) * M + M**3
        w = (abs(r) + sqrt(q**3 + r**2))**(2/3)
        E1 = (2 * r * w / (w**2 + w*q + q**2) + M)/d
        f2 = e * sin(E1)
        f3 = e * cos(E1)
        f0 = E1 - f2 - M
        f1 = 1 - f3
        δ3 = -f0 / (f1 - f0 * f2 / (2 * f1))
        δ4 = -f0 / (f1 + f2 * δ3 / 2 + δ3**2 * f3/6)
        δ5 = -f0 / (f1 + δ4*f2/2 + δ4**2*f3/6 - δ4**3*f2/24)
        E = E1 + δ5
        E = E % (2 * pi)

        return E


    def __calc_rrdot(self, x, y, z, vx, vy, vz):
        '''Calculate the dot product of the position and velocity vectors'''
        return dot([x, y, z], [vx, vy, vz])


    def __calc_true_anomaly_from_kep(self, e, E):
        '''Calculate true anomaly from keplerian elements'''
        return 2 * arctan2(sqrt(1 + e) * sin(E/2), sqrt(1 - e) * cos(E/2))


    def __calc_true_anomaly_from_xyz(self, x, y, z, evec, rrdot):
        '''Calculate true anomaly'''
        true_anomaly = arccos(dot(evec, [x, y, z]) / (norm(evec) * norm([x, y, z])))
        true_anomaly[rrdot < 0] = 2 * pi - true_anomaly[rrdot < 0]
        return true_anomaly


    def __calc_r_from_kep(self, a, e, E):
        '''Compure the distance of the rock from the center of the coordinate system'''
        return a * (1 - e * cos(E))

    def __calc_r_from_xyz(self, x, y, z):
        '''Compure the distance of the rock from the center of the coordinate system'''
        return norm([x, y, z])


    def __calc_ovec(self, r, true_anomaly):
        '''Compute the position of the object in its orbital plane.'''
        return r * array([cos(true_anomaly), sin(true_anomaly), zeros_like(true_anomaly)])


    def __calc_vovec(self, a, e, E, r, mu):
        '''Compute the velocity of the object in its orbital plane.'''
        return (mu * a**2 / r) * array([-sin(E), sqrt(1 - e**2) * cos(E), zeros_like(E)])


    def __calc_xyz(self, arg, inc, node, ovec):
        '''Compute state vector position from Keplerian elements.'''
        return euler_rotation(arg, inc, node, ovec)


    def __calc_vxyz(self, arg, inc, node, vovec):
        '''Compute state vector velocity from Keplerian elements.'''
        return euler_rotation(arg, inc, node, vovec)


    def __calc_hvec(self, x, y, z, vx, vy, vz):
        '''Compute an object's angular momentum vector from its state vector.'''
        return cross([x, y, z], [vx, vy, vz])


    def __calc_evec(self, x, y, z, vx, vy, vz, hvec, mu):
        '''Compute an object's eccentricity vector from its state and angular momentum vectors.'''
        return array(cross([vx, vy, vz], hvec)) / mu  - array([x, y, z]) / norm([x, y, z])


    def __calc_nvec(self, hvec):
        '''Compute an object's node vector from its angular momentum vector.'''
        return array([-hvec[1], hvec[0], zeros_like(hvec[2])])


    def __calc_inc(self, hvec):
        '''Compute orbital inclination'''
        return arccos(hvec[2] / norm(hvec))


    def __calc_arg(self, nvec, evec):
        '''Compute argument of pericenter'''
        e = norm(evec)
        n = norm(nvec)

        arg = arccos(dot(nvec, evec) / (n * e))
        arg[evec[2] < 0] = 2 * pi - arg[evec[2] < 0]
        arg[e == 0] = 0
        return arg


    def __calc_node(self, nvec, inc):
        '''Compute longitude of ascending node'''
        # compute ascending node
        node = arccos(nvec[0] / norm(nvec))
        node[inc == 0] = 0
        node[nvec[1] < 0] = 2 * pi - node[nvec[1] < 0]
        return node


    def __calc_M(self, e, E):
        '''Calculate mean anomaly from eccentric anomaly'''
        M = E - e * sin(E)
        return M % (2 * pi)


    def __calc_e(self, evec):
        return norm(evec)


    def __calc_a(self, x, y, z, vx, vy, vz, mu):
        return 1 / (2 / norm([x, y, z]) - norm([vx, vy, vz])**2 / mu)


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
            self.xyz_to_kep()
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
            self.xyz_to_kep()
            self.__class__.frame = 'heliocentric'

        return self


    def kep_to_xyz(self):

        '''Transform from Keplerian to Cartesian coordinates.'''

        self.E = self.__calc_E(self.e, self.M.rad)

        self.true_anomaly = self.__calc_true_anomaly_from_kep(self.e, self.E)

        self.r = self.__calc_r_from_kep(self.a, self.e, self.E)

        self.ovec = self.__calc_ovec(self.r, self.true_anomaly)

        self.vovec = self.__calc_vovec(self.a, self.e, self.E, self.r, self.mu)

        self.x, self.y, self.z = self.__calc_xyz(self.arg, self.inc, self.node, self.ovec)

        self.vx, self.vy, self.vz = self.__calc_vxyz(self.arg, self.inc, self.node, self.vovec)


    def xyz_to_kep(self):

        '''Transform from Cartesian to Keplerian coordinates.'''

        self.r = self.__calc_r_from_xyz(self.x, self.y, self.z)

        rrdot = self.__calc_rrdot(self.x, self.y, self.z, self.vx, self.vy, self.vz)

        self.hvec = self.__calc_hvec(self.x, self.y, self.z, self.vx, self.vy, self.vz)

        self.evec = self.__calc_evec(self.x, self.y, self.z, self.vx, self.vy, self.vz, self.hvec, self.mu)

        self.nvec = self.__calc_nvec(self.hvec)

        self.true_anomaly = self.__calc_true_anomaly_from_xyz(self.x, self.y, self.z, self.evec, self.rrdot)

        self.inc = self.__calc_inc(self.hvec)

        self.arg = self.__calc_arg(self.nvec, self.evec)

        self.node = self.__calc_node(self.nvec, self.inc)

        self.e = self.__calc_e(self.evec)

        self.a = self.__calc_a(self.x, self.y, self.z, self.vx, self.vy, self.vz, self.mu)


    '''Keplerian Elements'''

    @property
    def a(self):
        return self._a

    @a.setter
    def a(self, value):
        self._a = value

    @a.deleter
    def a(self):
        del self._a


    @property
    def e(self):
        return self._e

    @e.setter
    def e(self, value):
        self._e = value

    @e.deleter
    def e(self):
        del self._e


    @property
    def inc(self):
        return self._inc

    @inc.setter
    def inc(self, value):
        self._inc = value

    @inc.deleter
    def inc(self):
        del self._inc


    @property
    def node(self):
        if not hasattr(self, '_node'):
            self.node = (self.varpi - self.arg) % (2 * pi * u.rad)
        return self._node

    @node.setter
    def node(self, value):
        self._node = value

    @node.deleter
    def node(self):
        del self._node


    @property
    def arg(self):
        if not hasattr(self, '_arg'):
            self.arg = (self.varpi - self.node) % (2 * pi * u.rad)
        return self._arg

    @arg.setter
    def arg(self, value):
        self._arg = value

    @arg.deleter
    def arg(self):
        del self._arg


    @property
    def varpi(self):
        if not hasattr(self, '_varpi'):
            self.varpi = (self.node + self.arg) % (2 * pi * u.rad)
        return self._varpi

    @varpi.setter
    def varpi(self, value):
        self._varpi = value

    @varpi.deleter
    def varpi(self):
        del self._varpi

    @property
    def M(self):
        if not hasattr(self, '_M'):
            self.M = self.E * self.e * sin(self.E)
        return self._M

    @M.setter
    def M(self, value):
        self._M = value

    @M.deleter
    def M(self):
        del self._M

    @property
    def E(self):
        return self._E

    @E.setter
    def E(self, value):
        self._E = value

    @E.deleter
    def E(self):
        del self._E

    @property
    def true_anomaly(self):
        return self._true_anomaly

    @true_anomaly.setter
    def true_anomaly(self, value):
        self._true_anomaly = value

    @true_anomaly.deleter
    def true_anomaly(self):
        del self._true_anomaly

    @property
    def true_longitude(self):
        if not hasattr(self, '_true_longitude'):
            self.true_longitude = (self.true_anomaly + self.varpi) % (2 * pi * u.rad)
        return self._true_longitude

    @true_longitude.setter
    def true_longitude(self, value):
        self._true_longitude = value

    @true_longitude.deleter
    def true_longitude(self):
        del self._true_longitude

    @property
    def mean_longitude(self):
        if not hasattr(self, '_true_longitude'):
            self.mean_longitude = (self.M + self.varpi) % (2 * pi * u.rad)
        return self._mean_longitude

    @mean_longitude.setter
    def mean_longitude(self, value):
        self._mean_longitude = value

    @mean_longitude.deleter
    def mean_longitude(self):
        del self._mean_longitude


    @property
    def t_peri(self):
        if not hasattr(self, '_t_peri'):
            self.t_peri = Time(self.epoch.jd * u.day - self.M / self.n, format='jd', scale='utc')
        return self._t_peri

    @t_peri.setter
    def t_peri(self, value):
        self._t_peri = value

    @t_peri.deleter
    def t_peri(self):
        del self._t_peri


    @property
    def b(self):
        if not hasattr(self, '_b'):
            self.b = self.a * sqrt(1 - self.e**2)
        return self._b

    @b.setter
    def b(self, value):
        self._b = value

    @b.deleter
    def b(self):
        del self._b


    @property
    def p(self):
        if not hasattr(self, '_p'):
            self.p = self.a * (1 - self.e**2)
        return self._p

    @p.setter
    def p(self, value):
        self._p = value

    @p.deleter
    def p(self):
        del self._p


    @property
    def q(self):
        if not hasattr(self, '_q'):
            self.q = self.a * (1 - self.e)
        return self._q

    @q.setter
    def q(self, value):
        self._q = value

    @q.deleter
    def q(self):
        del self._q


    @property
    def Q(self):
        if not hasattr(self, '_Q'):
            self.Q = self.a * (1 + self.e)
        return self._Q

    @Q.setter
    def Q(self, value):
        self._Q = value

    @Q.deleter
    def Q(self):
        del self._Q


    @property
    def n(self):
        if not hasattr(self, '_n'):
            self.n = sqrt(self.mu / self.a**3)
        return self._n

    @n.setter
    def n(self, value):
        self._n = value

    @n.deleter
    def n(self):
        del self._n


    @property
    def mass(self):
        if not hasattr(self, '_mass'):
            self.mass = 0
        return self._mass

    @mass.setter
    def mass(self, value):
        self._mass = value


    @property
    def radius(self):
        if not hasattr(self, '_radius'):
            self.radius = 1e-16
        return self._radius

    @radius.setter
    def radius(self, value):
        self._radius = value


    @property
    def density(self):
        if not hasattr(self, '_density'):
            self.density = 1e-16
        return self._density

    @density.setter
    def density(self, value):
        self._density = value


    @property
    def hill_radius(self):
        if not hasattr(self, '_hill_radius'):
            self.hill_radius = self.a * (1 - self.e) * cbrt(self.m / 3)
        return self._hill_radius

    @hill_radius.setter
    def hill_radius(self, value):
        self._hill_radius = value

    @hill_radius.deleter
    def hill_radius(self):
        del self._hill_radius


    @property
    def TisserandJ(self):
        if not hasattr(self, '_TisserandJ'):
            aJ = 5.2 * u.au
            self.TisserandJ = aJ / self.a * 2 * cos(self.inc) * sqrt(self.p / aJ)
        return self._TisserandJ

    @TisserandJ.setter
    def TisserandJ(self, value):
        self._TisserandJ = value

    @TisserandJ.deleter
    def TisserandJ(self):
        del self._TisserandJ
