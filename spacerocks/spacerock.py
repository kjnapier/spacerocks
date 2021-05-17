################################################################################
# SpaceRocks, version 1.0.0
#
# Author: Kevin Napier kjnapier@umich.edu
################################################################################
import random

from astropy import units as u
from astropy.coordinates import Angle, Distance
from astropy.time import Time

from numpy import sin, cos, arctan2, sqrt, array, pi, zeros
import numpy as np
import pandas as pd

import rebound
import reboundx
from reboundx import constants

from .constants import *
from .orbitfuncs import OrbitFuncs
from .convenience import Convenience
from .units import Units
from .vector import Vector

from skyfield.api import Topos, Loader
# Load in planets for ephemeride calculation.
load = Loader('./Skyfield-Data', expire=False, verbose=False)
ts = load.timescale()
planets = load('de423.bsp')

class SpaceRock(OrbitFuncs, Convenience):

    def __init__(self, input_frame='barycentric', units=Units(), *args, **kwargs):

        coords = self.detect_coords(kwargs)
        input_frame = input_frame.lower()

        # input -> arrays
        for idx, key in enumerate([*kwargs]):
            if not hasattr(kwargs.get(key), '__len__'):
                kwargs[key] = array([kwargs.get(key)])
            else:
                kwargs[key] = array(kwargs.get(key))

        SpaceRock.frame = input_frame
        if SpaceRock.frame == 'barycentric':
            SpaceRock.mu = mu_bary
        elif SpaceRock.frame == 'heliocentric':
            SpaceRock.mu = mu_helio

        if units.timeformat is None:
            self.epoch = self.detect_timescale(kwargs.get('epoch'), units.timescale)
        else:
            self.epoch = Time(kwargs.get('epoch'), format=units.timeformat, scale=units.timescale)

        self.t0 = Time(self.epoch.jd, format='jd', scale=units.timescale)

        if kwargs.get('name') is not None:
            self.name = kwargs.get('name').astype(str)
        else:
            # produces random, non-repeting integers between 0 and 1e10 - 1
            self.name = array(['{:010}'.format(value) for value in random.sample(range(int(1e10)), len(self.epoch))])



        if coords == 'kep':

            self.a = Distance(kwargs.get('a'), units.distance)
            self.e = kwargs.get('e') #* u.dimensionless_unscaled
            self.inc = Angle(kwargs.get('inc'), units.angle)

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

            if kwargs.get('E') is not None:
                self.E = Angle(kwargs.get('E'), units.angle)

            if kwargs.get('true_anomaly') is not None:
                self.true_anomaly = Angle(kwargs.get('true_anomaly'), units.angle)

            if kwargs.get('true_longitude') is not None:
                self.true_longitude = Angle(kwargs.get('true_longitude'), units.angle)

            if kwargs.get('mean_longitude') is not None:
                self.mean_longitude = Angle(kwargs.get('mean_longitude'), units.angle)


        elif coords == 'xyz':

            x = Distance(kwargs.get('x'), units.distance, allow_negative=True)
            y = Distance(kwargs.get('y'), units.distance, allow_negative=True)
            z = Distance(kwargs.get('z'), units.distance, allow_negative=True)
            vx = (kwargs.get('vx') * units.speed).to(u.au / u.day)
            vy = (kwargs.get('vy') * units.speed).to(u.au / u.day)
            vz = (kwargs.get('vz') * units.speed).to(u.au / u.day)

            self.position = Vector(x, y, z)
            self.velocity = Vector(vx, vy, vz)


    def propagate(self, obsdates, model, add_pluto=False, gr=False):
        '''
        Integrate all bodies to the desired date. The logic could be cleaner
        but it works.
        '''
        Nx = len(obsdates)
        Ny = len(self)
        x_values = zeros([Nx, Ny])
        y_values = zeros([Nx, Ny])
        z_values = zeros([Nx, Ny])
        vx_values = zeros([Nx, Ny])
        vy_values = zeros([Nx, Ny])
        vz_values = zeros([Nx, Ny])
        obsdate_values = np.zeros([Nx, Ny])
        name_values = np.tile(self.name, Nx)
        in_frame = self.frame

        # We need to (or should) work in barycentric coordinates in rebound
        if in_frame == 'heliocentric':
            self.to_bary()

        # Integrate all particles to the same obsdate
        pickup_times = self.epoch.jd #df.epoch
        sim = self.set_simulation(np.min(pickup_times), model=model, add_pluto=add_pluto, gr=gr)
        sim.t = np.min(pickup_times) #np.min(df.epoch)

        for time in np.sort(np.unique(pickup_times)):
            ps = self[self.epoch.jd == time] #df[df.epoch == time]
            for p in ps:
                sim.add(x=p.x.value, y=p.y.value, z=p.z.value,
                        vx=p.vx.value, vy=p.vy.value, vz=p.vz.value,
                        hash=p.name)

                sim.integrate(time, exact_finish_time=1)

        for ii, time in enumerate(np.sort(obsdates)):
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

        rocks = SpaceRock(x=x, y=y, z=z, vx=vx, vy=vy, vz=vz, name=name, epoch=epoch, input_frame='barycentric')

        # be polite and return orbital parameters in the input frame.
        if in_frame == 'heliocentric':
            rocks.to_helio()

        #try:
        #    self.t0 = np.tile(self.t0, Nx)
        #    self.H = np.tile(self.H, Nx)
        #    self.G = np.tile(self.G, Nx)
        #    self.delta_H = np.tile(self.delta_H, Nx)
        #    self.rotation_period = np.tile(self.rotation_period, Nx)
        #    self.phi0 = np.tile(self.phi_0, Nx)
        #    #self.t0 = np.tile(self.t0, Nx)
        #except:
        #    pass

        '''
        Should make this return a dictionary?
        Maybe a sortby argument?
        '''

        return rocks

    def observe(self, rocks, obscode):
        pass

    def set_simulation(self, startdate, model, add_pluto=False, gr=False):

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

        M_mercury = 1.6601367952719304e-7 # Mercury Barycenter
        M_venus = 2.4478383396645447e-6 # Venus Barycenter
        M_earth = 3.0034896161241036e-06
        M_moon = M_earth / 81.3005690769
        M_mars = 3.2271560375549977e-7 # Mars Barycenter
        M_sun = 1
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
            active_bodies = [sun, earth, jupiter, saturn, uranus, neptune]
            names = ['Sun', 'Earth', 'Jupiter', 'Saturn', 'Uranus', 'Neptune']
            M_sun += M_mercury + M_venus + M_mars
            masses = [M_sun, M_earth, M_jupiter, M_saturn, M_uranus, M_neptune]


        elif model == 3:
            active_bodies = [sun, earth, mars, jupiter, saturn, uranus, neptune]
            names = ['Sun', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune']
            M_sun += M_mercury + M_venus
            masses = [M_sun, M_earth, M_mars, M_jupiter, M_saturn, M_uranus, M_neptune]


        elif model == 4:
            active_bodies = [sun, venus, earth, mars, jupiter, saturn, uranus, neptune]
            names = ['Sun', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune']
            M_sun += M_mercury
            masses = [M_sun, M_venus, M_earth, M_mars, M_jupiter, M_saturn, M_uranus, M_neptune]


        elif model == 5:
            active_bodies = [sun, mercury, venus, earth, mars, jupiter, saturn, uranus, neptune]
            names = ['Sun', 'Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune']
            masses = [M_sun, M_mercury, M_venus, M_earth, M_mars, M_jupiter, M_saturn, M_uranus, M_neptune]


        elif model == 6:
            active_bodies = [sun, mercury, venus, earth, moon, mars, jupiter, saturn, uranus, neptune]
            names = ['Sun', 'Mercury', 'Venus', 'Earth', 'Moon', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune']
            masses = [M_sun, M_mercury, M_venus, M_earth, M_moon, M_mars, M_jupiter, M_saturn, M_uranus, M_neptune]


        else:
            raise ValueError('Model not recognized. Check the documentation.')

        if add_pluto == True:
            pluto = planets['pluto barycenter']
            active_bodies.append(pluto)
            names.append('Pluto')
            masses.append(M_pluto)


        startdate = Time(startdate, scale='utc', format='jd')
        t = ts.tt(jd=startdate.tt.jd)

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


        #if Propagate.gr == True:
        #    bodies = sim.particles
        #    rebx = reboundx.Extras(sim)
        #    gr = rebx.load_force('gr')
        #    gr.params['c'] = constants.C
        #    rebx.add_force(gr)
        #    bodies['Sun'].params['gr_source'] = 1

        if gr == True:
            rebx = reboundx.Extras(sim)
            gr = rebx.load_force('gr_full')
            gr.params["c"] = constants.C
            rebx.add_force(gr)

        sim.testparticle_type = 0
        sim.integrator = 'ias15'
        #sim.ri_ias15.epsilon = 1e-12

        sim.move_to_com()

        return sim

    def orbits(self):

        M = np.linspace(0, 2*np.pi, 1000)

        xs = []
        ys = []
        zs = []

        for r in self:
            x, y, z, _, _, _ = self.kep_to_xyz_temp(r.a, r.e, r.inc, r.arg, r.node, M)
            xs.append(x)
            ys.append(y)
            zs.append(z)

        return xs, ys, zs


    def kep_to_xyz_temp(self, a, e, inc, arg, node, M):
        '''
        Just compute the xyz position of an object. Used for iterative equatorial
        calculation.
        '''
        # compute eccentric anomaly E
        M = array(M)
        M[M > pi] -= 2 * pi
        alpha = (3 * pi**2 + 1.6 * (pi**2 - pi * abs(M))/(1 + e))/(pi**2 - 6)
        d = 3 * (1 - e) + alpha * e
        q = 2 * alpha * d * (1 - e) - M**2
        r = 3 * alpha * d * (d - 1 + e) * M + M**3
        w = (abs(r) + sqrt(q**3 + r**2))**(2/3)
        E1 = (2 * r * w / (w**2 + w*q + q**2) + M)/d
        f2 = e * sin(E1)
        f3 = e * cos(E1)
        f0 = E1 - f2 - M
        f1 = 1 - f3
        d3 = -f0 / (f1 - f0 * f2 / (2 * f1))
        d4 = -f0 / (f1 + f2 * d3 / 2 + d3**2 * f3/6)
        d5 = -f0 / (f1 + d4*f2/2 + d4**2*f3/6 - d4**3*f2/24)
        E = E1 + d5
        E = E % (2 * pi)

        # compute true anomaly ν
        true_anomaly = 2 * np.arctan2((1 + e)**0.5*np.sin(E/2), (1 - e)**0.5*np.cos(E/2))

        # compute the distance to the central body r
        r = a * (1 - e * np.cos(E))

        # obtain the position vector o
        o = Vector(r * cos(true_anomaly), r * sin(true_anomaly), np.zeros_like(true_anomaly))
        ov = Vector((mu_bary * a)**0.5 / r * (-np.sin(E))/ u.rad, (mu_bary * a)**0.5 / r * ((1 - e**2)**0.5 * np.cos(E))/ u.rad, np.zeros(len(true_anomaly))/ u.rad)

        # Rotate o to the inertial frame
        position = o.euler_rotation(arg, inc, node) #* u.au
        velocity = ov.euler_rotation(arg, inc, node) #* u.au / u.day

        return position.x, position.y, position.z, velocity.x, velocity.y, velocity.z