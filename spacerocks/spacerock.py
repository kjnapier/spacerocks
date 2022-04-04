import uuid

from astropy import units as u
from astropy.coordinates import Angle, Distance
from astropy.time import Time

import datetime

import numpy as np
import pandas as pd

import rebound
import asdf
import copy

from .constants import mu_bary, c, epsilon
from .keplerorbit import KeplerOrbit
from .convenience import Convenience
from .units import Units
from .vector import Vector
from .ephemerides import Ephemerides 
from .observer import Observer
from .cbindings import kepM_to_xyz, correct_for_ltt
from .spice import SpiceBody
import os
import pkg_resources
import spiceypy as spice


SPICE_PATH = pkg_resources.resource_filename('spacerocks', 'data/spice')
spice.furnsh(os.path.join(SPICE_PATH, 'latest_leapseconds.tls'))
spice.furnsh(os.path.join(SPICE_PATH, 'de440s.bsp'))
#spice.furnsh(os.path.join(SPICE_PATH, 'hst.bsp'))
#spice.furnsh(os.path.join(SPICE_PATH, 'nh.bsp'))

sun = SpiceBody(spiceid='Sun')
earth = SpiceBody(spiceid='Earth')

DATA_PATH = pkg_resources.resource_filename('spacerocks', 'data/observatories.csv')
observatories = pd.read_csv(DATA_PATH)


class SpaceRock(KeplerOrbit, Convenience):
    '''
    SpaceRock objects provide an interface to work with Keplerian orbits. 
    '''
    def __init__(self, origin='ssb', frame='eclipJ2000', units=Units(), *args, **kwargs):

        coords = self.detect_coords(kwargs)
        #frame = frame.lower()
        origin = origin.lower()

        # input -> arrays
        for key in [*kwargs]:
            k = kwargs.get(key)
            if not isinstance(k, np.ndarray):
                if isinstance(k, list):
                    if key == 'name':
                        kwargs[key] = np.array(k, dtype=object)
                    else:
                        kwargs[key] = np.array(k)
                else:
                    kwargs[key] = np.atleast_1d(k)

        self.frame = frame

        self.origin = origin

        if self.origin == 'ssb':
            self.mu = mu_bary
        else:
            o = SpiceBody(spiceid=self.origin)
            self.mu = o.mu

        if coords == 'kep':

            if kwargs.get('a') is not None:
                self.a = Distance(kwargs.get('a'), units.distance, allow_negative=True)

            if kwargs.get('b') is not None:
                self.b = Distance(kwargs.get('b'), units.distance, allow_negative=True)

            if kwargs.get('v_inf') is not None:
                self.v_inf = (kwargs.get('v_inf') * units.speed).to(u.au / u.day)

            if kwargs.get('e') is not None:
                self.e = kwargs.get('e')
                if np.any(self.e < 0):
                    raise ValueError('Eccentricity must be positive')
                if hasattr(self, '_a'):
                    if np.any((self.a.au < 0) * (self.e < 1)):
                        raise ValueError('Orbital elements mismatch. a must be positive for e < 1.') 
                    if np.any((self.a.au > 0) * (self.e > 1)):
                        raise ValueError('Orbital elements mismatch. a must be negative for e > 1.') 

            self.inc = Angle(kwargs.get('inc'), units.angle)

            if kwargs.get('q') is not None:
                self.q = Distance(kwargs.get('q'), units.distance)
                if np.any(self.q.au < 0):
                    raise ValueError('pericenter distance must be positive')

            if kwargs.get('Q') is not None:
                self.Q = Distance(kwargs.get('Q'), units.distance)
                if np.any(self.Q.au < 0):
                    raise ValueError('apocenter distance must be positive')
                if np.any(self.Q < self.q):
                    raise ValueError('Apocenter distance must be >= pericenter distance.')

            if kwargs.get('node') is not None:
                self.node = Angle(kwargs.get('node'), units.angle)

            if kwargs.get('arg') is not None:
                self.arg = Angle(kwargs.get('arg'), units.angle)
                if np.any(self.arg.deg[self.e < 1e-8] != 0):
                    self.arg[self.e < 1e-8] = 0
                    #raise ValueError('Cannot have e = 0 with arg != 0')

            if kwargs.get('varpi') is not None:
                self.varpi = Angle(kwargs.get('varpi'), units.angle)

            if kwargs.get('t_peri') is not None:
                if units.timeformat is None:
                    self.t_peri = self.detect_timescale(kwargs.get('t_peri'), units.timescale)
                else:
                    self.t_peri = Time(kwargs.get('t_peri'), format=units.timeformat, scale=units.timescale)

            if kwargs.get('M') is not None:
                self.M = Angle(kwargs.get('M'), units.angle)

            if kwargs.get('E') is not None:
                self.E = Angle(kwargs.get('E'), units.angle)

            if kwargs.get('f') is not None:
                self.f = Angle(kwargs.get('f'), units.angle)

            if kwargs.get('true_longitude') is not None:
                self.true_longitude = Angle(kwargs.get('true_longitude'), units.angle)

            if kwargs.get('mean_longitude') is not None:
                self.mean_longitude = Angle(kwargs.get('mean_longitude'), units.angle)

            if kwargs.get('epoch') is not None:
                if units.timeformat is None:
                    self.epoch = self.detect_timescale(kwargs.get('epoch'), units.timescale)
                else:
                    self.epoch = Time(kwargs.get('epoch'), format=units.timeformat, scale=units.timescale)
            else:
                self.epoch = Time(np.zeros(len(self.inc)), format='jd', scale='utc')

            if kwargs.get('name') is not None:
                self.name = np.atleast_1d(kwargs.get('name'))
            else:
                # produces random, non-repeting integers between 0 and 1e10 - 1
                #self.name = np.array(['{:010}'.format(value) for value in random.sample(range(int(1e10)), len(self.inc))])
                self.name = np.array([uuid.uuid4().hex for _ in range(len(self.inc))])


        elif coords == 'xyz':

            x = Distance(kwargs.get('x'), units.distance, allow_negative=True).to(u.au)
            y = Distance(kwargs.get('y'), units.distance, allow_negative=True).to(u.au)
            z = Distance(kwargs.get('z'), units.distance, allow_negative=True).to(u.au)
            vx = (kwargs.get('vx') * units.speed).to(u.au / u.day)
            vy = (kwargs.get('vy') * units.speed).to(u.au / u.day)
            vz = (kwargs.get('vz') * units.speed).to(u.au / u.day)

            self.position = Vector(x, y, z)
            self.velocity = Vector(vx, vy, vz)

            if kwargs.get('epoch') is not None:

                if units.timeformat is None:
                    self.epoch = self.detect_timescale(kwargs.get('epoch'), units.timescale)
                else:
                    self.epoch = Time(kwargs.get('epoch'), format=units.timeformat, scale=units.timescale)
            else:
                self.epoch = Time(np.zeros(len(self.x)), format='jd', scale='utc')

            if kwargs.get('name') is not None:
                self.name = np.atleast_1d(kwargs.get('name'))
            else:
                # produces random, non-repeting integers between 0 and 1e10 - 1
                #self.name = np.array(['{:010}'.format(value) for value in random.sample(range(int(1e10)), len(self.x))])
                self.name = np.array([uuid.uuid4().hex for _ in range(len(self.x))])

        
        if kwargs.get('H') is not None:
            
            if units.rotation_curves == False:
                self.H = kwargs.get('H')

                if kwargs.get('G') is not None:
                    self.G = kwargs.get('G')
                else:
                    self.G = np.repeat(0.15, len(self))
            else:
                curves = kwargs.get('H')
                curve_funcs = []
                for curve in curves:
                    if callable(curve):
                        curve_funcs.append(curve)
                    else:
                        curve_funcs.append(lambda _, x=curve: x)
                self.H_func = np.array(curve_funcs)

                if kwargs.get('G') is not None:
                    self.G = kwargs.get('G')
                else:
                    self.G = np.repeat(0.15, len(self))

            

        elif kwargs.get('mag') is not None:
            curves = kwargs.get('mag')
            curve_funcs = []
            for curve in curves:
                if callable(curve):
                    curve_funcs.append(curve)
                else:
                    curve_funcs.append(lambda _, x=curve: x)
            self.mag_func = np.array(curve_funcs)

            if kwargs.get('G') is not None:
                self.G = kwargs.get('G')
            else:
                self.G = np.repeat(0.15, len(self))

        if kwargs.get('radius') is not None:
            self.radius = Distance(kwargs.get('radius'), units.size, allow_negative=False)

        if kwargs.get('mass') is not None:
            self.mass = kwargs.get('mass') * units.mass

        if kwargs.get('density') is not None:
            self.density = kwargs.get('density') * units.density


    @classmethod
    def from_horizons(cls, name):

        from urllib.parse import urlencode
        from urllib.request import urlopen
       

        def quote(text):
            return "'{}'".format(text)

        start = datetime.date.today()
        end = start + datetime.timedelta(days=1)

        xs = []
        ys = []
        zs = []
        vxs = []
        vys = []
        vzs = []
        Hs = []
        Gs = []
        epochs = []
        names = []

        for n in np.atleast_1d(name):
        
            params = {
                      'format':      "text",
                      'COMMAND':     quote(n),
                      'START_TIME':  quote(start),
                      'STOP_TIME':   quote(end),
                      'MAKE_EPHEM':  quote("YES"),
                      'EPHEM_TYPE':  quote("VECTORS"),
                      'CENTER':      quote("@sun"),
                      'REF_PLANE':   quote("ecliptic"),
                      'STEP_SIZE':   quote("1"),
                      'REF_SYSTEM':  quote("J2000"),
                      'VEC_CORR':    quote("NONE"),
                      'OUT_UNITS':   quote("KM-S"),
                      'CSV_FORMAT':  quote("NO"),
                      'VEC_DELTA_T': quote("NO"),
                      'VEC_TABLE':   quote("3"),
                      'VEC_LABELS':  quote("NO")
                     }
    
            try:
                url = "https://ssd.jpl.nasa.gov/api/horizons.api?" + urlencode(params)
                with urlopen(url) as f:
                    body = f.read().decode()

                lines = body.split("$$SOE")[-1].split("\n")
                x, y, z = [float(i) for i in lines[2].split()]
                vx, vy, vz = [float(i) for i in lines[3].split()]
        
                try:
                    pps = body.split('physical parameters')[-1].split('\n')[2].split()
                    H = float(pps[1])
                    G = float(pps[3])
                except Exception as E:
                    print(f'Magnitude information not available for {name}. Is this a comet?')
                    H = None
                    G = None
                    pass
    
                xs.append(x)
                ys.append(y)
                zs.append(z)
                vxs.append(vx)
                vys.append(vy)
                vzs.append(vz)
                Hs.append(H)
                Gs.append(G)
                epochs.append(start.ctime())
                names.append(n)

            except Exception as E:
                raise ValueError(f'Body {n} not found.')
            
        units = Units()
        units.distance = u.km
        units.speed = u.km/u.s
        units.timescale = 'tdb'

        return cls(x=xs, y=ys, z=zs, vx=vxs, vy=vys, vz=vzs, name=names,
                   origin='sun', frame='eclipJ2000', units=units, H=Hs, G=Gs, epoch=epochs)
        
    @classmethod
    def from_file(cls, name):
        with asdf.open(name) as f:
            N = len(f['x'])
            if len(f['name']) == 1:
                names = np.repeat(f['name'], N)
            else:
                names = f['name']
            units = Units()
            units.timescale = 'tdb'
            return cls(x=f['x'], 
                       y=f['y'], 
                       z=f['z'], 
                       vx=f['vx'], 
                       vy=f['vy'], 
                       vz=f['vz'], 
                       epoch=f['epoch'], 
                       name=names, 
                       origin=f['origin'], 
                       frame=f['frame'], 
                       units=units)

    @property
    def H(self):
        if hasattr(self, 'H_func'):
            return np.array([func(epoch) for epoch, func in zip(self.epoch.jd, self.H_func)])
        elif hasattr(self, '_H'):
            return self._H

    @H.setter
    def H(self, value):
        self._H = value

    @property
    def mag(self):
        if hasattr(self, 'mag_func'):
            return np.array([func(epoch) for epoch, func in zip(self.epoch.jd, self.mag_func)])

    def propagate(self, epochs, model, units=Units(), gr=False):
        '''
        Numerically integrate all bodies to the desired date.
        This routine synchronizes the epochs.
        '''

        epochs = self.detect_timescale(np.atleast_1d(epochs), units.timescale)
        origin = copy.copy(self.origin)
        frame = copy.copy(self.frame)

        # We need to integrate in barycentric coordinates
        self.to_bary()

        # Put the simulation into ecliptic coordinates
        self.change_frame('eclipJ2000')

        # Integrate all particles to the same obsdate
        pickup_times = self.epoch.tdb.jd
        sim, planet_names = self.set_simulation(np.min(pickup_times), model=model, gr=gr)
        sim.t = np.min(pickup_times)

        # need to ensure these are computed
        self.x
        self.y
        self.z
        self.vx
        self.vy
        self.vz

        sim.move_to_com()

        for time in np.sort(np.unique(pickup_times)):
            ps = self[self.epoch.tdb.jd == time]
            for x, y, z, vx, vy, vz, name in zip(ps.x.value, ps.y.value, ps.z.value, ps.vx.value, ps.vy.value, ps.vz.value, ps.name):
                sim.add(x=x, y=y, z=z, vx=vx, vy=vy, vz=vz, hash=name)
                sim.integrate(time, exact_finish_time=1)

        Nx = len(epochs)
        x_values = np.zeros((Nx, sim.N))
        y_values = np.zeros((Nx, sim.N))
        z_values = np.zeros((Nx, sim.N))
        vx_values = np.zeros((Nx, sim.N))
        vy_values = np.zeros((Nx, sim.N))
        vz_values = np.zeros((Nx, sim.N))
        name_values = np.zeros((Nx, sim.N), dtype='object')
        obsdate_values = np.zeros((Nx, sim.N))

        p_names = planet_names + [name for name in self.name]

        for ii, time in enumerate(np.sort(epochs.tdb.jd)):
            sim.integrate(time, exact_finish_time=1)
            a = np.zeros((sim.N, 3), dtype='float64')
            b = np.zeros((sim.N, 3), dtype='float64')
            sim.serialize_particle_data(xyz=a, vxvyvz=b)

            x, y, z = a.T
            vx, vy, vz = b.T

            x_values[ii] = x
            y_values[ii] = y
            z_values[ii] = z
            vx_values[ii] = vx
            vy_values[ii] = vy
            vz_values[ii] = vz
            name_values[ii] = p_names
            obsdate_values[ii] = np.repeat(sim.t, sim.N)

        Nactive = sim.N_active
        x = x_values[:, Nactive:].flatten()
        y = y_values[:, Nactive:].flatten()
        z = z_values[:, Nactive:].flatten()
        vx = vx_values[:, Nactive:].flatten()
        vy = vy_values[:, Nactive:].flatten()
        vz = vz_values[:, Nactive:].flatten()
        name = name_values[:, Nactive:].flatten()
        epoch = obsdate_values[:, Nactive:].flatten()

        px = x_values[:, :Nactive].flatten()
        py = y_values[:, :Nactive].flatten()
        pz = z_values[:, :Nactive].flatten()
        pvx = vx_values[:, :Nactive].flatten()
        pvy = vy_values[:, :Nactive].flatten()
        pvz = vz_values[:, :Nactive].flatten()
        pname = name_values[:, :Nactive].flatten()
        pepoch = obsdate_values[:, :Nactive].flatten()

        units = Units()
        units.timescale = 'tdb'
        rocks = self.__class__(x=x, y=y, z=z, vx=vx, vy=vy, vz=vz, name=name, epoch=epoch, origin='ssb', frame='eclipJ2000', units=units)
        planets = self.__class__(x=px, y=py, z=pz, vx=pvx, vy=pvy, vz=pvz, name=pname, epoch=pepoch, origin='ssb', frame='eclipJ2000', units=units)

        # be polite and return orbital parameters using the input origin.
        if origin != 'ssb':
            self.change_origin(origin)
            rocks.change_origin(origin)
            planets.change_origin(origin)

        if frame != 'eclipJ2000':
            self.change_frame(frame)
            rocks.change_frame(frame)
            planets.change_frame(frame)

        if hasattr(self, 'G'):
            rocks.G = np.tile(self.G, Nx)

        if hasattr(self, 'H_func'):
            rocks.H_func = np.tile(self.H_func, Nx)
        elif hasattr(self, '_H'):
            rocks.H = np.tile(self.H, Nx)

        if hasattr(self, 'mag_func'):
            rocks.mag_func = np.tile(self.mag_func, Nx)

        return rocks, planets, sim

    def calc_H(self, obscode):
        return self.__calc_H_from_mag(obscode=obscode)

    def __calc_H_from_mag(self, obscode):
        obs = self.observe(obscode=obscode)

        e = earth.at(self.epoch)
        s = sun.at(self.epoch)

        x_helio = self.x - s.x
        y_helio = self.y - s.y
        z_helio = self.z - s.z

        r_helio = Distance(np.sqrt(x_helio*x_helio + y_helio*y_helio + z_helio*z_helio).value, u.au)
        earth_dist = ((e.x-s.x)**2 + (e.y-s.y)**2 + (e.z-s.z)**2)**0.5

        q = (r_helio.au**2 + obs.delta.au**2 - earth_dist.value**2)/(2 * r_helio.au * obs.delta.au)

        beta = np.arccos(q)
        beta[np.where(q <= -1)[0]] = np.pi * u.rad
        beta[np.where(q >= 1)[0]] = 0 * u.rad

        Psi_1 = np.exp(-3.332 * np.tan(beta/2)**0.631)
        Psi_2 = np.exp(-1.862 * np.tan(beta/2)**1.218)

        H = self.mag - 5 * np.log10(r_helio.au * obs.delta.au)

        not_zero = np.where((Psi_1 != 0) | (Psi_2 != 0))[0]
        H[not_zero] += 2.5 * np.log10((1 - self.G[not_zero]) * Psi_1[not_zero] + self.G[not_zero] * Psi_2[not_zero])

        return H


    def observe(self, **kwargs) -> Ephemerides:

        '''
        Calculate the ephemerides of a SpaceRock object from a observer's location.

        All MPC observatory codes are supported, as well as the Geocenter.

        Currently-supported spacecraft are:
        - The Hubble Space Telescope
        - The New Horizons Spacecraft

        Ephemerides can also be computed as observed from the Sun and the 
        solar system barycenter (ssb).

        The James Webb Space Telescope will be supported as soon as it launches
        and NASA provides the necessary spk files.
        '''

        if kwargs.get('obscode') is not None:
            x, y, z, vx, vy, vz = self.xyz_to_tel(obscode=kwargs.get('obscode'))
        elif kwargs.get('spiceid') is not None:
            x, y, z, vx, vy, vz = self.xyz_to_tel(spiceid=kwargs.get('spiceid'))
        else:
            raise ValueError('Must pass either an obscode or spiceid.')

        if not (hasattr(self, 'H_func') or hasattr(self, '_H')):
            return Ephemerides(x=x, y=y, z=z, vx=vx, vy=vy, vz=vz, epoch=self.epoch, name=self.name)
        else:
            if hasattr(self, 'H_func'):
                return Ephemerides(x=x, y=y, z=z, vx=vx, vy=vy, vz=vz, epoch=self.epoch, name=self.name, H_func=self.H_func, G=self.G)
            else:
                return Ephemerides(x=x, y=y, z=z, vx=vx, vy=vy, vz=vz, epoch=self.epoch, name=self.name, H=self.H, G=self.G)


    def xyz_to_tel(self, **kwargs):
        '''
        Transform from barycentric ecliptic Cartesian coordinates to
        telescope-centric coordinates.

        Routine corrects iteratively for light travel time.
        '''

        if kwargs.get('obscode') is not None:
            observer = Observer(obscode=kwargs.get('obscode'), epoch=self.epoch.utc.jd)
        elif kwargs.get('spiceid') is not None:
            observer = Observer(spiceid=kwargs.get('spiceid'), epoch=self.epoch.utc.jd)
        else:
            raise ValueError('Must pass either an obscode or spiceid.')

        in_origin = copy.copy(self.origin)
        self.to_bary()

        in_frame = copy.copy(self.frame)
        self.change_frame('eclipJ2000') 

        dx, dy, dz, dvx, dvy, dvz = correct_for_ltt(self, observer)

        # Be polite
        if in_origin != self.origin:
            self.change_origin(in_origin)

        if in_frame != self.frame:
            self.change_frame(in_frame)

        # Transform to the equatorial frame
        yrot = dy * np.cos(epsilon) - dz * np.sin(epsilon)
        zrot = dy * np.sin(epsilon) + dz * np.cos(epsilon)
        vyrot = dvy * np.cos(epsilon) - dvz * np.sin(epsilon)
        vzrot = dvy * np.sin(epsilon) + dvz * np.cos(epsilon)

        dy = yrot
        dz = zrot
        dvy = vyrot
        dvz = vzrot

        return dx, dy, dz, dvx, dvy, dvz

    def set_simulation(self, startdate, model, gr):

        sun = SpiceBody(spiceid='Sun')
        mercury = SpiceBody(spiceid='Mercury Barycenter')
        venus = SpiceBody(spiceid='Venus Barycenter')
        earth = SpiceBody(spiceid='Earth')
        moon = SpiceBody(spiceid='Moon')
        mars = SpiceBody(spiceid='Mars Barycenter')
        jupiter = SpiceBody(spiceid='Jupiter Barycenter')
        saturn = SpiceBody(spiceid='Saturn Barycenter')
        uranus = SpiceBody(spiceid='Uranus Barycenter')
        neptune = SpiceBody(spiceid='Neptune Barycenter')
        pluto = SpiceBody(spiceid='Pluto Barycenter')

        M_sun = sun.mass.value
        M_mercury = mercury.mass.value
        M_venus = venus.mass.value
        M_earth = earth.mass.value
        M_moon = moon.mass.value
        M_mars = mars.mass.value
        M_jupiter = jupiter.mass.value
        M_saturn = saturn.mass.value
        M_uranus = uranus.mass.value
        M_neptune = neptune.mass.value
        M_pluto = pluto.mass.value

        if model == 0:
            active_bodies = [sun]
            names = ['Sun']
            M_sun += M_mercury + M_venus + M_earth + M_mars \
                     + M_jupiter + M_saturn + M_uranus + M_neptune
            masses = [M_sun]

        elif model == 1:
            active_bodies = [sun, jupiter, saturn, uranus, neptune]
            names = ['Sun', 'Jupiter', 'Saturn', 'Uranus', 'Neptune']
            M_sun += M_mercury + M_venus + M_earth + M_mars
            masses = [M_sun, M_jupiter, M_saturn, M_uranus, M_neptune]


        elif model == 2:
            active_bodies = [sun, mercury, venus, earth, moon, mars, jupiter, saturn, uranus, neptune, pluto]
            names = ['Sun', 'Mercury', 'Venus', 'Earth', 'Moon', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune', 'Pluto']
            masses = [M_sun, M_mercury, M_venus, M_earth, M_moon, M_mars, M_jupiter, M_saturn, M_uranus, M_neptune, M_pluto]

        elif model == 3:

            spice.furnsh(os.path.join(SPICE_PATH, 'sb441-n16s.bsp'))
            # spice.furnsh(os.path.join(SPICE_PATH, 'asteroids', '2000001.bsp'))
            # spice.furnsh(os.path.join(SPICE_PATH, 'asteroids', '2000002.bsp'))
            # spice.furnsh(os.path.join(SPICE_PATH, 'asteroids', '2000003.bsp'))
            # spice.furnsh(os.path.join(SPICE_PATH, 'asteroids', '2000004.bsp'))
            # spice.furnsh(os.path.join(SPICE_PATH, 'asteroids', '2000007.bsp'))
            # spice.furnsh(os.path.join(SPICE_PATH, 'asteroids', '2000010.bsp'))
            # spice.furnsh(os.path.join(SPICE_PATH, 'asteroids', '2000015.bsp'))
            # spice.furnsh(os.path.join(SPICE_PATH, 'asteroids', '2000016.bsp'))
            # spice.furnsh(os.path.join(SPICE_PATH, 'asteroids', '2000031.bsp'))
            # spice.furnsh(os.path.join(SPICE_PATH, 'asteroids', '2000052.bsp'))
            # spice.furnsh(os.path.join(SPICE_PATH, 'asteroids', '2000065.bsp'))
            # spice.furnsh(os.path.join(SPICE_PATH, 'asteroids', '2000087.bsp'))
            # spice.furnsh(os.path.join(SPICE_PATH, 'asteroids', '2000088.bsp'))
            # spice.furnsh(os.path.join(SPICE_PATH, 'asteroids', '2000107.bsp'))
            # spice.furnsh(os.path.join(SPICE_PATH, 'asteroids', '2000511.bsp'))
            # spice.furnsh(os.path.join(SPICE_PATH, 'asteroids', '2000704.bsp'))

            ceres = SpiceBody(spiceid='Ceres')
            vesta = SpiceBody(spiceid='Vesta')
            pallas = SpiceBody(spiceid='Pallas')
            juno = SpiceBody(spiceid='2000003')
            iris = SpiceBody(spiceid='2000007')
            hygiea = SpiceBody(spiceid='2000010')
            eunomia = SpiceBody(spiceid='2000015')
            psyche = SpiceBody(spiceid='2000016')
            euphrosyne = SpiceBody(spiceid='2000031')
            europa = SpiceBody(spiceid='2000052')
            cybele = SpiceBody(spiceid='2000065')
            sylvia = SpiceBody(spiceid='2000087')
            thisbe = SpiceBody(spiceid='2000088')
            camilla = SpiceBody(spiceid='2000107')
            davida = SpiceBody(spiceid='2000511')
            interamnia = SpiceBody(spiceid='2000704')
            
            M_vesta = vesta.mass.value
            M_ceres = ceres.mass.value
            M_pallas = pallas.mass.value
            M_juno = juno.mass.value
            M_iris = iris.mass.value
            M_hygiea = hygiea.mass.value
            M_eunomia = eunomia.mass.value
            M_psyche = psyche.mass.value
            M_europa = europa.mass.value
            M_cybele = cybele.mass.value
            M_sylvia = sylvia.mass.value
            M_thisbe = thisbe.mass.value
            M_davida = davida.mass.value
            M_interamnia = interamnia.mass.value
            M_camilla = camilla.mass.value
            M_euphrosyne = euphrosyne.mass.value

            active_bodies = [sun, mercury, venus, earth, moon, mars, jupiter, saturn, uranus, neptune, pluto, 
                             ceres, vesta, pallas, interamnia, juno, camilla, iris, hygiea, eunomia, psyche, 
                             euphrosyne, europa, cybele, sylvia, thisbe, davida]
            names = ['Sun', 'Mercury', 'Venus', 'Earth', 'Moon', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune', 'Pluto', 
                     'Ceres', 'Vesta', 'Pallas', 'Interamnia', 'Juno', 'Camilla', 'Iris', 'Hygiea', 'Eunomia', 
                     'Psyche', 'Euphrosyne', 'Europa', 'Cybele', 'Sylvia', 'Thisbe', 'Davida']
            masses = [M_sun, M_mercury, M_venus, M_earth, M_moon, M_mars, M_jupiter, M_saturn, M_uranus, M_neptune, M_pluto, 
                      M_ceres, M_vesta, M_pallas, M_interamnia, M_juno, M_camilla, M_iris, M_hygiea, M_eunomia, M_psyche, 
                             M_euphrosyne, M_europa, M_cybele, M_sylvia, M_thisbe, M_davida]


        else:
            raise ValueError('Model not recognized. Check the documentation.')

        startdate = Time(startdate, scale='tdb', format='jd')

        bodies = [body.at(startdate) for body in active_bodies]

        # create a dataframe of the massive bodies in the solar system
        ss = pd.DataFrame()
        ss['x'] = [body.x.au[0] for body in bodies]
        ss['y'] = [body.y.au[0] for body in bodies]
        ss['z'] = [body.z.au[0] for body in bodies]
        ss['vx'] = [body.vx.value[0] for body in bodies]
        ss['vy'] = [body.vy.value[0] for body in bodies]
        ss['vz'] = [body.vz.value[0] for body in bodies]
        ss['mass'] = masses
        ss['name'] = names

        sim = rebound.Simulation()
        sim.units = ('day', 'AU', 'Msun')

        for p in ss.itertuples():
            sim.add(x=p.x, y=p.y, z=p.z,
                    vx=p.vx, vy=p.vy, vz=p.vz,
                    m=p.mass, hash=p.name)

        sim.N_active = len(ss)

        sim.testparticle_type = 0
        sim.integrator = 'ias15'
        #sim.ri_ias15.epsilon = 1e-12
        

        if gr == True:
            import reboundx
            rebx = reboundx.Extras(sim)
            gr_force = rebx.load_force('gr_full')
            gr_force.params["c"] = c.to(u.au / u.day).value
            rebx.add_force(gr_force)
        
        return sim, names

    def analytic_propagate(self, epoch, units=Units()):
        '''
        analytically propagate all rocks to a common epoch
        '''
        epoch = self.detect_timescale(np.atleast_1d(epoch), units.timescale)
        dt = (epoch.utc.jd - self.epoch.utc.jd) * u.day
        dM = self.n * dt
        return self.__class__(a=self.a.au, 
                              e=self.e, 
                              inc=self.inc.deg, 
                              node=self.node.deg, 
                              arg=self.arg.deg, 
                              M=(self.M + dM).deg, 
                              name=self.name, 
                              epoch=np.repeat(epoch.utc.jd, len(self)), 
                              origin=self.origin, 
                              frame=self.frame)
        

    def orbits(self, N=1000):

        M = Angle(np.linspace(0, 2*np.pi, N), u.rad)

        xs = []
        ys = []
        zs = []

        for r in self:
            x, y, z, _, _, _ = kepM_to_xyz(np.repeat(r.a, N),
                                           np.repeat(r.e, N),
                                           np.repeat(r.inc.rad, N),
                                           np.repeat(r.arg.rad, N),
                                           np.repeat(r.node.rad, N),
                                           M.rad)
            xs.append(x)
            ys.append(y)
            zs.append(z)

        return xs, ys, zs

    def write_to(self, path, compression='zlib'):
        uniquenames = np.unique(self.name)
        if len(uniquenames) == 1:
            name = uniquenames
        else:
            name = self.name.tolist()
        tree = {
            'frame': self.frame,
            'origin': self.origin,
            'mu': self.mu,
            'epoch': self.epoch.tdb.jd,
            'name': name,
            'x': self.x.au,
            'y': self.y.au,
            'z': self.z.au,
            'vx': self.vx.to(u.au/u.day).value,
            'vy': self.vy.to(u.au/u.day).value,
            'vz': self.vz.to(u.au/u.day).value
        }

        # Create the ASDF file object from our data tree
        af = asdf.AsdfFile(tree)
        af.write_to(path, all_array_compression=compression)