import numpy as np
from astropy import units as u
from astropy.coordinates import Angle, Distance
from astropy.time import Time
from .constants import mu_bary, frames

from .vector import Vector
from .cbindings import calc_kep_from_xyz, calc_vovec_from_kep, calc_M_from_E, calc_E_from_f, calc_E_from_M, calc_f_from_E
from .spice import SpiceBody

class KeplerOrbit:

    def change_origin(self, new_origin: str):
        '''
        Method to convert coordinates to a new origin.
        Modifies the SpaceRock object inplace.

        new_origin: string of the spiceid of the new origin.
        '''
        if self.origin != new_origin:

            current_origin = SpiceBody(spiceid=self.origin)
            new = SpiceBody(spiceid=new_origin)
            co = current_origin.at(self.epoch)
            no = new.at(self.epoch)

            x = self.x + co.x - no.x
            y = self.y + co.y - no.y
            z = self.z + co.z - no.z
            vx = self.vx + co.vx - no.vx
            vy = self.vy + co.vy - no.vy
            vz = self.vz + co.vz - no.vz

            self.origin = new_origin

            if new_origin == 'ssb':
                self.mu = mu_bary
            else:
                self.mu = new.mu

            # clear the keplerian variables because they need to be recomputed
            self.clear_kep()

            self.position = Vector(x, y, z)
            self.velocity = Vector(vx, vy, vz)

        return self
        

    def change_frame(self, new_frame: str):
        '''
        Method to convert coordinates to a new frame.
        Modifies the SpaceRock object inplace.

        new_frame: string of the spiceid of the new frame.
        '''

        if self.frame != new_frame:

            rotation_matrix = frames[new_frame.upper()] @ np.linalg.inv(frames[self.frame.upper()])

            distance_units = self.x._unit
            velo_units = self.vx._unit

            x, y, z = rotation_matrix.dot([self.x, self.y, self.z])
            vx, vy, vz = rotation_matrix.dot([self.vx, self.vy, self.vz]) * velo_units

            self.frame = new_frame
    
            # clear the keplerian variables because they need to be recomputed
            self.clear_kep()

            self.x = Distance(x, distance_units, allow_negative=True)
            self.y = Distance(y, distance_units, allow_negative=True)
            self.z = Distance(z, distance_units, allow_negative=True)
            self.vx = vx
            self.vy = vy
            self.vz = vz

    
            # self.position = Vector(Distance(x, distance_units, allow_negative=True), 
            #                        Distance(y, distance_units, allow_negative=True), 
            #                        Distance(z, distance_units, allow_negative=True))
            # self.velocity = Vector(vx, vy, vz)

        return self


    def to_bary(self):
        '''
        Method to convert heliocentric coordinates to barycentric coordinates.
        '''
        self.change_origin('ssb')
        

    def to_helio(self):
        '''
        Method to convert barycentric coordinates to heliocentric coordinates.
        '''
        self.change_origin('sun')
       

    def clear_kep(self):

        to_delete = ['_a', 
                     '_e', 
                     '_inc', 
                     '_arg', 
                     '_node', 
                     '_varpi', 
                     '_M', 
                     '_E',
                     '_f', 
                     '_true_longitude', 
                     '_mean_longitude',
                     '_q', 
                     '_t_peri', 
                     '_b', 
                     '_p', 
                     '_n', 
                     '_Q', 
                     '_hill_radius', 
                     '_r',
                     '_ovec', 
                     '_vovec', 
                     '_position', 
                     '_velocity', 
                     '_rrdot',
                     '_x', 
                     '_y', 
                     '_z', 
                     '_vx', 
                     '_vy', 
                     '_vz']

        for attr in to_delete:
            self.__dict__.pop(attr, None)

    def kep_from_xyz(self):
        self.a, self.e, self.inc, self.arg, self.node, self.f = calc_kep_from_xyz(self.mu.value, self.x.au, self.y.au, self.z.au, self.vx.value, self.vy.value, self.vz.value)

    ''' Vector Quantities '''

    @property
    def ovec(self):
        if not hasattr(self, '_ovec'):
            self.ovec = Vector(self.r * np.cos(self.f), self.r * np.sin(self.f), np.zeros_like(self.f))
        return self._ovec

    @ovec.setter
    def ovec(self, value):
        self._ovec = value

    @ovec.deleter
    def ovec(self):
        del self._ovec

    @property
    def vovec(self):
        if not hasattr(self, '_vovec'):
            vx, vy, vz = calc_vovec_from_kep(self.mu.value, self.a.au, self.e, self.r.au, self.E.rad)
            self.vovec = Vector(vx, vy, vz)
        return self._vovec

    @vovec.setter
    def vovec(self, value):
        self._vovec = value

    @vovec.deleter
    def vovec(self):
        del self._vovec

    @property
    def position(self):
        if not hasattr(self, '_position'):
            if hasattr(self, '_x') and hasattr(self, '_y') and hasattr(self, '_z'):
                self.position = Vector(self.x, self.y, self.z)
            else:
                self.position = self.ovec.euler_rotation(self.arg, self.inc, self.node)
        return self._position

    @position.setter
    def position(self, value):
        self._position = value

    @position.deleter
    def position(self):
        del self._position

    @property
    def velocity(self):
        if not hasattr(self, '_velocity'):
            if hasattr(self, '_vx') and hasattr(self, '_vy') and hasattr(self, '_vz'):
                self.velocity = Vector(self.vx, self.vy, self.vz)
            else:
                self.velocity = self.vovec.euler_rotation(self.arg, self.inc, self.node)
        return self._velocity

    @velocity.setter
    def velocity(self, value):
        self._velocity = value

    @velocity.deleter
    def velocity(self):
        del self._velocity

    ''' Cartesian Elements '''

    @property
    def x(self):
        if not hasattr(self, '_x'):
            self.x = self.position.x
        return self._x

    @x.setter
    def x(self, value):
        self._x = value

    @x.deleter
    def x(self):
        del self._x

    @property
    def y(self):
        if not hasattr(self, '_y'):
            self.y = self.position.y
        return self._y

    @y.setter
    def y(self, value):
        self._y = value

    @y.deleter
    def y(self):
        del self._y

    @property
    def z(self):
        if not hasattr(self, '_z'):
            self.z = self.position.z
        return self._z

    @z.setter
    def z(self, value):
        self._z = value

    @z.deleter
    def z(self):
        del self._z

    @property
    def vx(self):
        if not hasattr(self, '_vx'):
            self.vx = self.velocity.x
        return self._vx

    @vx.setter
    def vx(self, value):
        self._vx = value

    @vx.deleter
    def vx(self):
        del self._vx

    @property
    def vy(self):
        if not hasattr(self, '_vy'):
            self.vy = self.velocity.y
        return self._vy

    @vy.setter
    def vy(self, value):
        self._vy = value

    @vy.deleter
    def vy(self):
        del self._vy

    @property
    def vz(self):
        if not hasattr(self, '_vz'):
            self.vz = self.velocity.z
        return self._vz

    @vz.setter
    def vz(self, value):
        self._vz = value

    @vz.deleter
    def vz(self):
        del self._vz

    ''' Required Keplerian Elements '''

    @property
    def a(self):
        if not hasattr(self, '_a'):
            if hasattr(self, '_position') and hasattr(self, '_velocity'):
                self.kep_from_xyz()
            elif hasattr(self, '_e') and hasattr(self, '_q'):
                self.a = self.q / (1 - self.e)
            elif hasattr(self, '_e') and hasattr(self, '_b'):
                self.a = self.b / \
                    np.sqrt(abs(1 - self.e**2)) * (-1 * (self.e > 1))
            elif hasattr(self, '_Q') and hasattr(self, '_q'):
                self.a = self.q / (1 - self.e)
            elif hasattr(self, '_v_inf'):
                self.a = Distance(-self.mu / self.v_inf**2 /
                                  u.rad**2, u.au, allow_negative=True)
        return self._a

    @a.setter
    def a(self, value):
        self._a = value

    @a.deleter
    def a(self):
        del self._a

    @property
    def e(self):
        if not hasattr(self, '_e'):

            if hasattr(self, '_position') and hasattr(self, '_velocity'):
                self.kep_from_xyz()

            elif hasattr(self, '_a') and hasattr(self, '_q'):
                self.e = 1 - self.q.au / self.a.au

            elif hasattr(self, '_Q') and hasattr(self, '_q'):
                self.e = (1 - self.q.au / self.Q.au) / (1 + self.q.au / self.Q.au)

            elif hasattr(self, '_b') and hasattr(self, '_v_inf'):
                self.e = np.sqrt(1 + self.b.au**2 / self.a.au**2)

        return self._e

    @e.setter
    def e(self, value):
        self._e = value

    @e.deleter
    def e(self):
        del self._e

    @property
    def inc(self):
        if not hasattr(self, '_inc'):
            self.kep_from_xyz()
        return self._inc

    @inc.setter
    def inc(self, value):
        self._inc = value

    @inc.deleter
    def inc(self):
        del self._inc

    ''' Keplerian choice-of-three '''

    @property
    def node(self):
        if not hasattr(self, '_node'):
            if hasattr(self, '_varpi') and hasattr(self, '_arg'):
                node = np.zeros(len(self))
                prograde = self.inc.deg < 90
                retrograde = ~prograde
                node[prograde] = self.varpi[prograde].rad - self.arg[prograde].rad
                node[retrograde] = self.varpi[retrograde].rad + self.arg[retrograde].rad
                self.node = Angle(node % (2 * np.pi), u.rad)

            elif hasattr(self, '_position') and hasattr(self, '_velocity'):
                self.kep_from_xyz()
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

            if hasattr(self, '_varpi') and hasattr(self, 'node'):
                arg = np.zeros(len(self))
                prograde = self.inc.deg < 90
                retrograde = ~prograde
                arg[prograde] = self.varpi[prograde].rad - self.node[prograde].rad
                arg[retrograde] = self.node[retrograde].rad - self.varpi[retrograde].rad

                #self.arg = Angle((self.varpi.rad - self.node.rad) % (2 * np.pi), u.rad)
                self.arg = Angle(arg % (2 * np.pi), u.rad)

            elif hasattr(self, '_position') and hasattr(self, '_velocity'):
                self.kep_from_xyz()

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
            #self.varpi = Angle((self.node.rad + self.arg.rad) %
            #                   (2 * np.pi), u.rad)

            varpi = np.zeros(len(self))
            prograde = self.inc.deg < 90
            retrograde = ~prograde
            varpi[prograde] = self.node[prograde].rad + self.arg[prograde].rad
            varpi[retrograde] = self.node[retrograde].rad - self.arg[retrograde].rad
            self.varpi = Angle(varpi % (2 * np.pi), u.rad)

        return self._varpi

    @varpi.setter
    def varpi(self, value):
        self._varpi = value

    @varpi.deleter
    def varpi(self):
        del self._varpi

    '''

    Keplerian choice-of-six

    There are several options for describing the location of a rock on its orbit:
    - mean anomaly
    - true anomaly
    - eccentric anomaly
    - mean longitude
    - true longitude
    - time of pericenter

    As long as you know the longitude of pericenter, you can get any of these
    quantities from any of the others.

    '''

    @property
    def M(self):
        if not hasattr(self, '_M'):

            if hasattr(self, '_mean_longitude'):
                M = (self.mean_longitude.rad - self.varpi.rad) % (2 * np.pi)
                self.M = Angle(M, u.rad)

            elif hasattr(self, '_f') or hasattr(self, '_true_longitude') or hasattr(self, '_E') or (hasattr(self, '_position') and hasattr(self, '_velocity')):
                self.M = calc_M_from_E(self.e, self.E.rad)

            elif hasattr(self, '_t_peri'):
                M = self.n * (self.epoch - self.t_peri)
                M = M % (2 * np.pi * u.rad)
                self.M = Angle(M, u.rad)

        return self._M

    @M.setter
    def M(self, value):
        self._M = value

    @M.deleter
    def M(self):
        del self._M

    @property
    def E(self):
        if not hasattr(self, '_E'):

            if hasattr(self, '_f') or hasattr(self, '_true_longitude') or (hasattr(self, '_position') and hasattr(self, '_velocity')):
                self.E = calc_E_from_f(self.e, self.f.rad)

            elif hasattr(self, '_M') or hasattr(self, '_mean_longitude') or hasattr(self, '_t_peri'):
                self.E = calc_E_from_M(self.e, self.M.rad)
        return self._E

    @E.setter
    def E(self, value):
        self._E = value

    @E.deleter
    def E(self):
        del self._E

    @property
    def f(self):
        if not hasattr(self, '_f'):

            if hasattr(self, '_true_longitude'):
                self.f = Angle(self.true_longitude.rad - self.varpi.rad, u.rad)

            elif hasattr(self, '_E') or hasattr(self, '_M') or hasattr(self, '_mean_longitude') or hasattr(self, '_t_peri'):
                self.f = calc_f_from_E(self.e, self.E.rad)

            elif hasattr(self, '_position') and hasattr(self, '_velocity'):
                self.kep_from_xyz()

        return self._f

    @f.setter
    def f(self, value):
        self._f = value

    @f.deleter
    def f(self):
        del self._f

    @property
    def true_longitude(self):
        if not hasattr(self, '_true_longitude') or (hasattr(self, '_position') and hasattr(self, '_velocity')):
            self.true_longitude = Angle((self.f.rad + self.varpi.rad) % (2 * np.pi), u.rad)
        return self._true_longitude

    @true_longitude.setter
    def true_longitude(self, value):
        self._true_longitude = value

    @true_longitude.deleter
    def true_longitude(self):
        del self._true_longitude

    @property
    def mean_longitude(self):
        if not hasattr(self, '_mean_longitude'):
            self.mean_longitude = Angle(
                (self.M.rad + self.varpi.rad) % (2 * np.pi), u.rad)
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

    ''' Physical Properties '''

    @property
    def mass(self):
        if not hasattr(self, '_mass'):
            if hasattr(self, '_radius') and hasattr(self, '_density'):
                self.mass = self.density * (4/3) * np.pi * self.radius**3
            else:
                raise ValueError('Insufficient information for calculating mass.')
        return self._mass

    @mass.setter
    def mass(self, value):
       self._mass = value

    @property
    def radius(self):
        if not hasattr(self, '_radius'):
            if hasattr(self, '_mass') and hasattr(self, '_density'):
               self.radius = Distance(np.cbrt((3 / (4 * np.pi * self.density)) * self.mass).to(u.km))
            else:
                raise ValueError('Insufficient information for calculating radius.')
        return self._radius

    @radius.setter
    def radius(self, value):
       self._radius = value

    @property
    def density(self):
        if not hasattr(self, '_density'):
            if hasattr(self, '_mass') and hasattr(self, '_radius'):
                self.density = (3 * self.mass) / (4 * np.pi * self.radius**3)
            else:
                raise ValueError('Insufficient information for calculating density.')
        return self._density

    @density.setter
    def density(self, value):
       self._density = value

    ''' Derived Quantities '''

    @property
    def b(self):
        if not hasattr(self, '_b'):
            self.b = abs(self.a) * np.sqrt(abs(1 - self.e**2))
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
            self.p = abs(self.a) * abs(1 - self.e**2)
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
    def v_inf(self):
        if not hasattr(self, '_v_inf'):
            self.v_inf = np.sqrt(- self.mu / self.a) / u.rad
        return self._v_inf

    @v_inf.setter
    def v_inf(self, value):
        self._v_inf = value

    @v_inf.deleter
    def v_inf(self):
        del self._v_inf

    @property
    def n(self):
        if not hasattr(self, '_n'):
            self.n = np.sqrt(self.mu / abs(self.a)**3)
        return self._n

    @n.setter
    def n(self, value):
        self._n = value

    @n.deleter
    def n(self):
        del self._n

    @property
    def rrdot(self):
        '''Calculate the dot product of the position and velocity vectors'''
        return self.position.dot(self.velocity)

    @property
    def r(self):
        '''Compute the distance of the rock from the center of the coordinate system'''
        if not hasattr(self, '_r'):
            if hasattr(self, '_position'):
                self.r = Distance(self.position.norm)
            else:
                rs = np.zeros(len(self))
                rs[self.e < 1] = self.a[self.e < 1] * \
                    (1 - self.e[self.e < 1] * np.cos(self.E[self.e < 1].rad))
                rs[self.e >= 1] = self.p[self.e >= 1] / \
                    (1 + self.e[self.e >= 1] *
                     np.cos(self.f[self.e >= 1]))
                self.r = Distance(rs, u.au)
        return self._r

    @r.setter
    def r(self, value):
        self._r = value

    @r.deleter
    def r(self):
        del self._r