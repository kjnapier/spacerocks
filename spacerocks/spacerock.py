################################################################################
# SpaceRocks, version 1.1.0
#
# Author: Kevin Napier kjnapier@umich.edu
################################################################################
import random

from astropy import units as u
from astropy.coordinates import Angle, Distance
from astropy.time import Time

from numpy import sqrt, array
import numpy as np
import pandas as pd

import rebound

#import reboundx
#from reboundx import constants

from .constants import *
from .keplerorbit import KeplerOrbit
from .convenience import Convenience
from .units import Units
from .vector import Vector
from .ephemerides import Ephemerides
from .observer import Observer
from .cbindings import kep_to_xyz_temp
from .spice import *
import os
import pkg_resources
import spiceypy as spice


SPICE_PATH = pkg_resources.resource_filename('spacerocks', 'data/spice')
spice.furnsh(os.path.join(SPICE_PATH, 'latest_leapseconds.tls'))
spice.furnsh(os.path.join(SPICE_PATH, 'de440s.bsp'))
spice.furnsh(os.path.join(SPICE_PATH, 'hst.bsp'))
spice.furnsh(os.path.join(SPICE_PATH, 'nh.bsp'))

sun = SpiceBody(spiceid='Sun')
earth = SpiceBody(spiceid='Earth')

DATA_PATH = pkg_resources.resource_filename('spacerocks', 'data/observatories.csv')
observatories = pd.read_csv(DATA_PATH)


class SpaceRock(KeplerOrbit, Convenience):

    '''
    SpaceRock objects provide an interface to work with Keplerian orbits. 
    '''

    def __init__(self, origin='ssb', frame='ECLIPJ2000', units=Units(), *args, **kwargs):

        coords = self.detect_coords(kwargs)
        frame = frame.lower()
        origin = origin.lower()

        # input -> arrays
        for key in [*kwargs]:
            kwargs[key] = np.atleast_1d(kwargs.get(key))

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
                        raise ValueError('Orbital elements mismatch. a must be negative for e < 1.') 

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
                    raise ValueError('Cannot have e = 0 with arg != 0')

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
                self.name = array(['{:010}'.format(value) for value in random.sample(range(int(1e10)), len(self.inc))])


        elif coords == 'xyz':

            x = Distance(kwargs.get('x'), units.distance, allow_negative=True)
            y = Distance(kwargs.get('y'), units.distance, allow_negative=True)
            z = Distance(kwargs.get('z'), units.distance, allow_negative=True)
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
                self.name = array(['{:010}'.format(value) for value in random.sample(range(int(1e10)), len(self.x))])

        if kwargs.get('G') is not None:
            self.G = kwargs.get('G')
        else:
            self.G = np.repeat(0.15, len(self))

        if kwargs.get('H') is not None:
            curves = kwargs.get('H')
            curve_funcs = []
            for curve in curves:
                if callable(curve):
                    curve_funcs.append(curve)
                else:
                    curve_funcs.append(lambda x: curve)
            self.H_func = np.array(curve_funcs)

        if kwargs.get('mag') is not None:
            curves = kwargs.get('mag')
            curve_funcs = []
            for curve in curves:
                if callable(curve):
                    curve_funcs.append(curve)
                else:
                    curve_funcs.append(lambda x: curve)
            self.mag_func = np.array(curve_funcs)

        if kwargs.get('radius') is not None:
            self.radius = Distance(kwargs.get('radius'), units.size, allow_negative=False)

        if kwargs.get('mass') is not None:
            self.mass = kwargs.get('mass') * units.mass

        if kwargs.get('density') is not None:
            self.density = kwargs.get('density') * units.density
        

    @property
    def H(self):
        return np.array([func(epoch) for epoch, func in zip(self.epoch.jd, self.H_func)])

    def propagate(self, epochs, model, units=Units()):
        '''
        Numerically integrate all bodies to the desired date.
        This routine synchronizes the epochs.
        '''

        epochs = self.detect_timescale(np.atleast_1d(epochs), units.timescale)
        origin = self.origin

        # We need to integrate in barycentric coordinates
        self.to_bary()

        # Integrate all particles to the same obsdate
        pickup_times = self.epoch.tdb.jd
        sim, planet_names = self.set_simulation(np.min(pickup_times), model=model)
        sim.t = np.min(pickup_times)

        # need to ensure these are computed
        self.x
        self.y
        self.z
        self.vx
        self.vy
        self.vz

        

        for time in np.sort(np.unique(pickup_times)):
            ps = self[self.epoch.tdb.jd == time]
            for x, y, z, vx, vy, vz, name in zip(ps.x.value, ps.y.value, ps.z.value, ps.vx.value, ps.vy.value, ps.vz.value, ps.name):
                sim.add(x=x, y=y, z=z, vx=vx, vy=vy, vz=vz, hash=name)
                sim.integrate(time, exact_finish_time=1)

        sim.move_to_com()

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
        rocks = self.__class__(x=x, y=y, z=z, vx=vx, vy=vy, vz=vz, name=name, epoch=epoch, origin='ssb', units=units)
        planets = self.__class__(x=px, y=py, z=pz, vx=pvx, vy=pvy, vz=pvz, name=pname, epoch=pepoch, origin='ssb', units=units)

        # be polite and return orbital parameters using the input origin.
        if origin == 'sun':
            rocks.to_helio()
            self.to_helio()

        if hasattr(self, 'G'):
            rocks.G = np.tile(self.G, Nx)

        if hasattr(self, 'H_func'):
            rocks.H_func = np.tile(self.H_func, Nx)

        if hasattr(self, 'mag_func'):
            rocks.mag_func = np.tile(self.mag_func, Nx)

        return rocks, planets, sim


    def __calc_H_from_mag(self, obscode):
        obs = self.observe(obscode=obscode)

        e = earth.at(self.epoch)
        s = sun.at(self.epoch)

       
        earth_dist = ((e.x-s.x)**2 + (e.y-s.y)**2 + (e.z-s.z)**2)**0.5

        q = (obs.r_helio.au**2 + obs.delta.au**2 - earth_dist**2)/(2 * obs.r_helio.au * obs.delta.au)

        beta = np.arccos(q)
        beta[np.where(q <= -1)[0]] = np.pi * u.rad
        beta[np.where(q >= 1)[0]] = 0 * u.rad

        Psi_1 = np.exp(-3.332 * np.tan(beta/2)**0.631)
        Psi_2 = np.exp(-1.862 * np.tan(beta/2)**1.218)


        H = self.mag - 5 * np.log10(obs.r_helio.au * obs.delta.au)

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


        if not hasattr(self, 'H_func'):
            return Ephemerides(x=x, y=y, z=z, vx=vx, vy=vy, vz=vz, epoch=self.epoch, name=self.name)
        else:
            return Ephemerides(x=x, y=y, z=z, vx=vx, vy=vy, vz=vz, epoch=self.epoch, name=self.name, H=self.H_func, G=self.G)

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

        in_origin = self.origin
        self.to_bary()

        x0, y0, z0 = self.x, self.y, self.z
        vx0, vy0, vz0 = self.vx, self.vy, self.vz

        ltt0 = 0

        a = self.a.au.astype(np.double)
        e = self.e.astype(np.double)
        inc = self.inc.rad.astype(np.double)
        arg = self.arg.rad.astype(np.double)
        node = self.node.rad.astype(np.double)

        for idx in range(10):

            dx = x0 - observer.x
            dy = y0 - observer.y
            dz = z0 - observer.z
            dvx = vx0 - observer.vx
            dvy = vy0 - observer.vy
            dvz = vz0 - observer.vz

            delta = sqrt(dx**2 + dy**2 + dz**2)

            ltt = delta / c
            dltt = (ltt - ltt0)

            if np.all(abs(max(dltt.value)) < 1e-12):
                break

            M = self.M - (ltt * self.n)
            x0, y0, z0, vx0, vy0, vz0 = kep_to_xyz_temp(a, e, inc, arg, node, M.rad.astype(np.double))

            ltt0 = ltt

        # Be polite
        if in_origin != self.origin:
            self.change_origin(in_origin)

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

    def set_simulation(self, startdate, model):

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
        ss['a'] = 1 / (2 / sqrt(ss.x**2 + ss.y**2 + ss.z**2) - (ss.vx**2 + ss.vy**2 + ss.vz**2) / mu_bary.value)
        ss['hill_radius'] = ss.a * pow(ss.mass / (3 * M_sun), 1/3)
        ss['name'] = names

        sim = rebound.Simulation()
        sim.units = ('day', 'AU', 'Msun')

        for p in ss.itertuples():
            sim.add(x=p.x, y=p.y, z=p.z,
                    vx=p.vx, vy=p.vy, vz=p.vz,
                    m=p.mass, hash=p.name)

        

        sim.N_active = len(ss)

        #if gr == True:
        #    rebx = reboundx.Extras(sim)
        #    gr = rebx.load_force('gr_full')
        #    gr.params["c"] = constants.C
        #    rebx.add_force(gr)

        sim.testparticle_type = 0
        sim.integrator = 'mercurius'
        

        return sim, names

    def orbits(self, N=1000):

        M = Angle(np.linspace(0, 2*np.pi, N), u.rad)

        xs = []
        ys = []
        zs = []

        for r in self:
            x, y, z, _, _, _ = kep_to_xyz_temp(np.repeat(r.a, N),
                                               np.repeat(r.e, N),
                                               np.repeat(r.inc.rad, N),
                                               np.repeat(r.arg.rad, N),
                                               np.repeat(r.node.rad, N),
                                               M.rad)
            xs.append(x)
            ys.append(y)
            zs.append(z)

        return xs, ys, zs

    def analytic_propagate(self, epoch: list, propagate_origin: str='sun'):
        '''
        propagate all bodies to the desired date using Keplerian orbit.
        '''
        in_origin = self.origin
        if propagate_origin != in_origin:
            if propagate_origin == 'sun':
                self.to_helio()
            else:
                self.to_bary()

        M = (self.n.value * (epoch - self.epoch.jd) + self.M.rad*180/np.pi)%360

        rocks = SpaceRock(a=self.a,
                          e=self.e,
                          inc=self.inc,
                          node=self.node,
                          arg=self.arg,
                          M=M,
                          name=self.name,
                          epoch=epoch,
                          origin=propagate_origin, 
                          frame=self.frame)

        # be polite and return orbital parameters in the input frame.
        if in_origin != self.origin:
            if in_origin == 'sun':
                self.to_helio()
            else:
                self.to_bary()

        if hasattr(self, 'G'):
            rocks.G = self.G

        if hasattr(self, 'mag'):
            rocks.mag = self.mag

        if hasattr(self, 'delta_H'):
            rocks.delta_H = self.delta_H
            rocks.rotation_period = self.rotation_period
            rocks.phi0 = self.phi0
            rocks.t0 = Time(self.t0.jd, format='jd')

            rocks.H = self.H + rocks.delta_H * np.sin(2 * np.pi * (rocks.epoch.jd - rocks.t0.jd) / rocks.rotation_period  - rocks.phi0)

        elif hasattr(self, 'H'):
            rocks.H = self.H

        return rocks