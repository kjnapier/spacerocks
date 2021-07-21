################################################################################
# SpaceRocks, version 1.0.7
#
# Author: Kevin Napier kjnapier@umich.edu
################################################################################
import random
import copy
import os

from astropy import units as u
from astropy.coordinates import Angle, Distance
from astropy.time import Time

from numpy import sin, cos, arctan2, sqrt, array, pi, zeros, sinh, arctanh, arccosh, cosh, tanh, zeros_like
import numpy as np
import pandas as pd

import rebound
import reboundx
from reboundx import constants
from scipy.optimize import newton

from .constants import *
from .orbitfuncs import OrbitFuncs
from .convenience import Convenience
from .units import Units
from .vector import Vector
from .ephemerides import Ephemerides

from skyfield.api import Topos, Loader
# Load in planets for ephemeride calculation.
load = Loader('./Skyfield-Data', expire=False, verbose=False)
ts = load.timescale()
planets = load('de440.bsp')
earth = planets['earth']
sun = planets['sun']

observatories = pd.read_csv(os.path.join(os.path.dirname(__file__),
                            'data',
                            'observatories.csv'))

from skyfield.api import wgs84
from skyfield.data import iers

ts = load.timescale()

class SpaceRock(OrbitFuncs, Convenience):

    def __init__(self, frame='barycentric', units=Units(), *args, **kwargs):

        coords = self.detect_coords(kwargs)
        frame = frame.lower()

        # input -> arrays
        for idx, key in enumerate([*kwargs]):
            kwargs[key] = np.atleast_1d(kwargs.get(key))


        self.frame = frame
        if self.frame == 'barycentric':
            self.mu = mu_bary
        elif self.frame == 'heliocentric':
            self.mu = mu_helio


        if coords == 'kep':

            if kwargs.get('a') is not None:
                self.a = Distance(kwargs.get('a'), units.distance, allow_negative=True)

            if kwargs.get('b') is not None:
                self.b = Distance(kwargs.get('b'), units.distance, allow_negative=True)

            if kwargs.get('v_inf') is not None:
                self.v_inf = (kwargs.get('v_inf') * units.speed).to(u.au / u.day)

            if kwargs.get('e') is not None:
                self.e = kwargs.get('e') #* u.dimensionless_unscaled

            self.inc = Angle(kwargs.get('inc'), units.angle)

            if kwargs.get('q') is not None:
                self.q = Distance(kwargs.get('q'), units.distance)

            if kwargs.get('Q') is not None:
                self.Q = Distance(kwargs.get('Q'), units.distance)

            if kwargs.get('node') is not None:
                self.node = Angle(kwargs.get('node'), units.angle)

            if kwargs.get('arg') is not None:
                self.arg = Angle(kwargs.get('arg'), units.angle)

            if kwargs.get('varpi') is not None:
                self.varpi = Angle(kwargs.get('varpi'), units.angle)

            if kwargs.get('t_peri') is not None:
                if units.timeformat is None:
                    self.t_peri = self.detect_timescale(kwargs.get('t_peri'), units.timescale)
                else:
                    self.t_peri = Time(kwargs.get('t_peri'), format=units.timeformat, scale=units.timescale)

            if kwargs.get('M') is not None:
                self.M = Angle(kwargs.get('M'), units.angle)
                #self.M[(self.e > 1) & (self.M.rad > np.pi)] -= 2 * np.pi * u.rad

            if kwargs.get('E') is not None:
                self.E = Angle(kwargs.get('E'), units.angle)

            if kwargs.get('true_anomaly') is not None:
                self.true_anomaly = Angle(kwargs.get('true_anomaly'), units.angle)

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


        if kwargs.get('H0') is not None:
            self.H0 = kwargs.get('H0')


        if kwargs.get('G') is not None:
            self.G = kwargs.get('G')
        else:
            self.G = np.repeat(0.15, len(self))


        if kwargs.get('mag') is not None:
            self.mag = kwargs.get('mag')


        if (kwargs.get('rotation_period') is not None) and (kwargs.get('delta_H') is not None) and (kwargs.get('phi0') is not None):
            self.rotation_period = kwargs.get('rotation_period')
            self.delta_H = kwargs.get('delta_H')
            self.phi0 = Angle(kwargs.get('phi0'), units.angle)
            self.t0 = Time(self.epoch.jd, format='jd', scale=units.timescale)


    def analytic_propagate(self, epoch, propagate_frame = 'heliocentric'):
        '''
        propagate all bodies to the desired date using Kieprian orbit.
        '''
        in_frame = self.frame
        if propagate_frame != in_frame:
            if propagate_frame == 'heliocentric':
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
                          frame=propagate_frame)

        # be polite and return orbital parameters in the input frame.
        if in_frame != self.frame:
            if in_frame == 'heliocentric':
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

    def propagate(self, epochs, model, units=Units(), gr=False):
        '''
        Integrate all bodies to the desired date. The logic could be cleaner
        but it works.
        '''

        epochs = self.detect_timescale(np.atleast_1d(epochs), units.timescale)

        Nx = len(epochs)
        Ny = len(self)
        x_values = zeros([Nx, Ny])
        y_values = zeros([Nx, Ny])
        z_values = zeros([Nx, Ny])
        vx_values = zeros([Nx, Ny])
        vy_values = zeros([Nx, Ny])
        vz_values = zeros([Nx, Ny])
        obsdate_values = np.zeros([Nx, Ny])
        name_values = np.tile(self.name, Nx)
        frame = self.frame

        # We need to integrate in barycentric coordinates
        self.to_bary()

        # Integrate all particles to the same obsdate
        pickup_times = self.epoch.tdb.jd
        sim = self.set_simulation(np.min(pickup_times), units=units, model=model, gr=gr)
        sim.t = np.min(pickup_times) #np.min(df.epoch)

        # need to ensure these are computed
        self.x
        self.y
        self.z
        self.vx
        self.vy
        self.vz

        for time in np.sort(np.unique(pickup_times)):
            ps = self[self.epoch.tdb.jd == time]
            for p in ps:
                sim.add(x=p.x.value, y=p.y.value, z=p.z.value,
                        vx=p.vx.value, vy=p.vy.value, vz=p.vz.value,
                        hash=p.name)

                sim.integrate(time, exact_finish_time=1)

        for ii, time in enumerate(np.sort(epochs.tdb.jd)):
            sim.integrate(time, exact_finish_time=1)
            for jj, name in enumerate(self.name):
                x_values[ii, jj] = sim.particles[name].x
                y_values[ii, jj] = sim.particles[name].y
                z_values[ii, jj] = sim.particles[name].z
                vx_values[ii, jj] = sim.particles[name].vx
                vy_values[ii, jj] = sim.particles[name].vy
                vz_values[ii, jj] = sim.particles[name].vz
                obsdate_values[ii, jj] = sim.t

        x = x_values.flatten()
        y = y_values.flatten()
        z = z_values.flatten()
        vx = vx_values.flatten()
        vy = vy_values.flatten()
        vz = vz_values.flatten()
        name = name_values.flatten()
        epoch = obsdate_values.flatten()

        units = Units()
        units.timescale = 'tdb'
        rocks = self.__class__(x=x, y=y, z=z, vx=vx, vy=vy, vz=vz, name=name, epoch=epoch, frame='barycentric', units=units)


        # be polite and return orbital parameters in the input frame.
        if frame == 'heliocentric':
            rocks.to_helio()
            self.to_helio()


        if hasattr(self, 'G'):
            rocks.G = np.tile(self.G, Nx)


        if hasattr(self, 'delta_H'):
            rocks.delta_H = np.tile(self.delta_H, Nx)
            rocks.rotation_period = np.tile(self.rotation_period, Nx)
            rocks.phi0 = np.tile(self.phi0, Nx)
            rocks.t0 = Time(np.tile(self.t0.jd, Nx), format='jd')
            rocks.H0 = np.tile(self.H0, Nx)


        elif hasattr(self, 'H0'):
            rocks.H0 = np.tile(self.H0, Nx)

        return rocks

    def calc_H(self, obscode):
        obs = self.observe(obscode=obscode)

        t = ts.tdb(jd=self.epoch.tdb.jd)
        e = earth.at(t)
        s = sun.at(t)

        sx, sy, sz = s.ecliptic_xyz().au
        ex, ey, ez = e.ecliptic_xyz().au
        earth_dist = ((ex-sx)**2 + (ey-sy)**2 + (ez-sz)**2)**0.5

        q = (obs.r_helio.au**2 + obs.delta.au**2 - earth_dist**2)/(2 * obs.r_helio.au * obs.delta.au)

        beta = np.arccos(q)
        beta[np.where(q <= -1)[0]] = np.pi * u.rad
        beta[np.where(q >= 1)[0]] = 0 * u.rad

        Psi_1 = np.exp(-3.332 * np.tan(beta/2)**0.631)
        Psi_2 = np.exp(-1.862 * np.tan(beta/2)**1.218)


        H = self.mag - 5 * np.log10(obs.r_helio.au * obs.delta.au)# / earth_dist**2)

        not_zero = np.where((Psi_1 != 0) | (Psi_2 != 0))[0]
        H[not_zero] += 2.5 * np.log10((1 - self.G[not_zero]) * Psi_1[not_zero] + self.G[not_zero] * Psi_2[not_zero])

        return H


    def observe(self, obscode):

        if self.frame == 'barycentric':
            self.to_helio()
            r_helio = self.r
            self.to_bary()
        else:
            r_helio = self.r

        x, y, z, vx, vy, vz = self.xyz_to_tel(obscode)

        if not hasattr(self, 'H0'):
            return Ephemerides(x=x, y=y, z=z, vx=vx, vy=vy, vz=vz, epoch=self.epoch, name=self.name, r_helio=r_helio)
        elif (not hasattr(self, 'delta_H')) and hasattr(self, 'H0'):
            return Ephemerides(x=x, y=y, z=z, vx=vx, vy=vy, vz=vz, epoch=self.epoch, name=self.name, r_helio=r_helio, H0=self.H0, G=self.G)
        else:
            return Ephemerides(x=x, y=y, z=z, vx=vx, vy=vy, vz=vz, epoch=self.epoch, name=self.name, r_helio=r_helio, H0=self.H0, G=self.G, delta_H=self.delta_H, t0=self.t0, rotation_period=self.rotation_period, phi0=self.phi0)

    def xyz_to_tel(self, obscode=500):
        '''
        Transform from barycentric ecliptic Cartesian coordinates to
        telescope-centric coordinates.

        Routine corrects iteratively for light travel time.
        '''

        if obscode != 500:
            obscode = str(obscode).zfill(3)
            obs = observatories[observatories.obscode == obscode]
            obslat = obs.lat.values
            obslon = obs.lon.values
            obselev = obs.elevation.values


        in_frame = self.frame
        #if in_frame == 'heliocentric':
        self.to_bary()


        alltimes = self.epoch.tdb.jd
        unique_times = np.unique(alltimes)
        t = ts.tdb(jd=unique_times)
        earth = planets['earth']



        # Only used for the topocentric calculation.
        if (obscode != 500) and (obscode != '500'):
            earth += wgs84.latlon(obslat, obslon, obselev)

        ee = earth.at(t)

        xx, yy, zz = ee.position.au * u.au # earth ICRS position
        vxx, vyy, vzz = ee.velocity.au_per_d * u.au / u.day # earth ICRS position


        exs = {t:x.value for t, x in zip(unique_times, xx)}
        eys = {t:y.value for t, y in zip(unique_times, yy)}
        ezs = {t:z.value for t, z in zip(unique_times, zz)}
        evxs = {t:vx.value for t, vx in zip(unique_times, vxx)}
        evys = {t:vy.value for t, vy in zip(unique_times, vyy)}
        evzs = {t:vz.value for t, vz in zip(unique_times, vzz)}

        x_earth = np.array([exs[t] for t in alltimes]) * u.au
        y_earth = np.array([eys[t] for t in alltimes]) * u.au
        z_earth = np.array([ezs[t] for t in alltimes]) * u.au

        vx_earth = np.array([evxs[t] for t in alltimes]) * u.au / u.day
        vy_earth = np.array([evys[t] for t in alltimes]) * u.au / u.day
        vz_earth = np.array([evzs[t] for t in alltimes]) * u.au / u.day


        x0, y0, z0 = self.x, self.y, self.z
        vx0, vy0, vz0 = self.vx, self.vy, self.vz

        ltt0 = 0

        for idx in range(10):

            dx = x0 - x_earth
            dy = y0 * np.cos(epsilon) - z0 * np.sin(epsilon) - y_earth
            dz = y0 * np.sin(epsilon) + z0 * np.cos(epsilon) - z_earth
            dvx = vx0 - vx_earth
            dvy = vy0 * np.cos(epsilon) - vz0 * np.sin(epsilon) - vy_earth
            dvz = vy0 * np.sin(epsilon) + vz0 * np.cos(epsilon) - vz_earth

            delta = sqrt(dx**2 + dy**2 + dz**2)

            ltt = delta / c
            dltt = (ltt - ltt0)

            if np.all(abs(max(dltt.value)) < 1e-12):
                break

            M = self.M - (ltt * self.n)
            x0, y0, z0, vx0, vy0, vz0 = self.kep_to_xyz_temp(self.a,
                                                             self.e,
                                                             self.inc,
                                                             self.arg,
                                                             self.node,
                                                             M)

            ltt0 = ltt

        # Be polite
        if in_frame == 'heliocentric':
            self.to_helio()

        return dx, dy, dz, dvx, dvy, dvz


    def set_simulation(self, startdate, units, model, gr=False):

        sun = planets['sun']
        mercury = planets['mercury barycenter']
        venus = planets['venus barycenter']
        earth = planets['earth']
        moon = planets['moon']
        mars = planets['mars barycenter']
        jupiter = planets['jupiter barycenter']
        saturn = planets['saturn barycenter']
        uranus = planets['uranus barycenter']
        neptune = planets['neptune barycenter']
        pluto = planets['pluto barycenter']

        M_sun = 1
        M_mercury = 1.6601367952719304e-7 # Mercury Barycenter
        M_venus = 2.4478383396645447e-6 # Venus Barycenter
        M_earth = 3.0034896161241036e-06
        M_moon = M_earth / 81.3005690769
        M_mars = 3.2271560375549977e-7 # Mars Barycenter

        M_jupiter = 9.547919384243222e-4
        M_saturn = 2.858859806661029e-4
        M_uranus = 4.3662440433515637e-5
        M_neptune = 5.151389020535497e-5
        M_pluto = 7.361781606089469e-9


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


        startdate = Time(startdate, scale=units.timescale, format='jd')
        t = ts.tdb(jd=startdate.tdb.jd)

        x, y, z = np.array([body.at(t).ecliptic_xyz().au for body in active_bodies]).T
        vx, vy, vz = np.array([body.at(t).ecliptic_velocity().au_per_d for body in active_bodies]).T

        # create a dataframe of the massive bodies in the solar system
        ss = pd.DataFrame()
        ss['x'] = x
        ss['y'] = y
        ss['z'] = z
        ss['vx'] = vx
        ss['vy'] = vy
        ss['vz'] = vz
        ss['mass'] = masses
        ss['a'] = 1 / (2 / sqrt(ss.x**2 + ss.y**2 + ss.z**2) - (ss.vx**2 + ss.vy**2 + ss.vz**2) / mu_bary.value)
        ss['hill_radius'] = ss.a * pow(ss.mass / (3 * M_sun), 1/3)
        ss['name'] = names


        sim = rebound.Simulation()
        sim.units = ('day', 'AU', 'Msun')

        for p in ss.itertuples():
            sim.add(x=p.x, y=p.y, z=p.z,
                    vx=p.vx, vy=p.vy, vz=p.vz,
                    m=p.mass, hash=p.name, r=p.hill_radius)

        sim.N_active = len(ss)

        if gr == True:
            rebx = reboundx.Extras(sim)
            gr = rebx.load_force('gr_full')
            gr.params["c"] = constants.C
            rebx.add_force(gr)

        sim.testparticle_type = 0
        sim.integrator = 'ias15'

        sim.move_to_com()

        return sim


    def orbits(self):

        M = Angle(np.linspace(0, 2*np.pi, 1000), u.rad)

        xs = []
        ys = []
        zs = []

        for r in self:
            x, y, z, _, _, _ = self.kep_to_xyz_temp(np.repeat(r.a, 1000),
                                                    np.repeat(r.e, 1000),
                                                    np.repeat(r.inc, 1000),
                                                    np.repeat(r.arg, 1000),
                                                    np.repeat(r.node, 1000),
                                                    M)
            xs.append(x)
            ys.append(y)
            zs.append(z)

        return xs, ys, zs


    def kep_to_xyz_temp(self, a, e, inc, arg, node, M):
        '''
        Just compute the xyz position of an object. Used for iterative
        equatorial calculation.
        '''

        E = np.zeros(len(M))

        MM = M[e < 1].rad
        ee = e[e < 1]

        MM[MM > pi] -= 2 * pi
        alpha = (3 * pi**2 + 1.6 * (pi**2 - pi * abs(MM))/(1 + ee))/(pi**2 - 6)
        d = 3 * (1 - ee) + alpha * ee
        q = 2 * alpha * d * (1 - ee) - MM**2
        r = 3 * alpha * d * (d - 1 + ee) * MM + MM**3
        w = (abs(r) + sqrt(q**3 + r**2))**(2/3)
        E1 = (2 * r * w / (w**2 + w*q + q**2) + MM)/d
        f2 = ee * sin(E1)
        f3 = ee * cos(E1)
        f0 = E1 - f2 - MM
        f1 = 1 - f3
        d3 = -f0 / (f1 - f0 * f2 / (2 * f1))
        d4 = -f0 / (f1 + f2 * d3 / 2 + d3**2 * f3 / 6)
        d5 = -f0 / (f1 + d4 * f2 / 2 + d4**2 * f3 / 6 - d4**3 * f2 / 24)
        E[e < 1] = (E1 + d5) % (2 * pi)

        if np.any(e >= 1):
            MMM = M[e >= 1].rad
            eee = e[e >= 1]
            f = lambda E, MMM, eee: eee * sinh(E) - E - MMM
            E0 = MMM
            E[e >= 1] = np.array([newton(f, E0[idx], args=(MMM[idx], eee[idx]), maxiter=10000) for idx in range(len(MMM))])

        # compute true anomaly ν
        #true_anomaly = 2 * np.arctan2((1 + e)**0.5*np.sin(E/2), (1 - e)**0.5*np.cos(E/2))

        true_anomaly = np.zeros(len(M))
        true_anomaly[e < 1] = 2 * arctan2(sqrt(1 + e[e < 1]) * sin(E[e < 1] / 2), sqrt(1 - e[e < 1]) * cos(E[e < 1] / 2))
        true_anomaly[e >= 1] = 2 * arctan2(sqrt(e[e >= 1] + 1) * tanh(E[e >= 1] / 2), sqrt(e[e >= 1] - 1))


        r = np.zeros(len(M))
        r[e < 1] = a[e < 1] * (1 - e[e < 1] * cos(E[e < 1]))
        r[e >= 1] = abs(a[e >= 1]) * abs(1 - e[e >= 1]**2) / (1 + e[e >= 1] * cos(true_anomaly[e >= 1]))

        # compute the distance to the central body r
        #r = a * (1 - e * np.cos(E))

        # obtain the position vector o
        o = Vector(r * cos(true_anomaly), r * sin(true_anomaly), np.zeros_like(true_anomaly))

        vx = np.zeros(len(M))
        vy = np.zeros(len(M))
        vz = np.zeros(len(M))
        xi_bound = sqrt(mu_bary / u.rad**2 * abs(a[e < 1])) / r[e < 1]
        xi_unbound = sqrt(mu_bary / u.rad**2 * abs(a[e >= 1])) / r[e >= 1]

        vx[e < 1] = - xi_bound * sin(E[e < 1])
        vx[e >= 1] = - xi_unbound * sinh(E[e >= 1])

        vy[e < 1] = xi_bound * sqrt(abs(1 - e[e < 1]**2)) * cos(E[e < 1])
        vy[e >= 1] = xi_unbound * sqrt(abs(1 - e[e >= 1]**2)) * np.cosh(E[e >= 1])

        vz[e < 1] = xi_bound * zeros_like(E[e < 1])
        vz[e >= 1] = xi_unbound * zeros_like(E[e >= 1])

        ov = Vector(vx, vy, vz)

        #ov = Vector((mu_bary * a)**0.5 / r * (-np.sin(E))/ u.rad, (mu_bary * a)**0.5 / r * ((1 - e**2)**0.5 * np.cos(E))/ u.rad, np.zeros(len(true_anomaly))/ u.rad)

        # Rotate o to the inertial frame
        position = o.euler_rotation(arg, inc, node) #* u.au
        velocity = ov.euler_rotation(arg, inc, node) #* u.au / u.day

        return position.x * u.au, position.y * u.au, position.z * u.au, velocity.x * u.au/u.d, velocity.y * u.au/u.d, velocity.z * u.au/u.d
