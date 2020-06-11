from spacerocks import SpaceRock
import numpy as np
import pandas as pd
import rebound
import reboundx
from reboundx import constants

import warnings

from astropy import units as u
from astropy.table import Table
from astropy.coordinates import Angle
from astropy.time import Time
from astropy.coordinates import SkyCoord

from skyfield.api import Topos, Loader
# Load in planets for ephemeride calculation.
load = Loader('./Skyfield-Data', expire=False, verbose=False)
ts = load.timescale()
planets = load('de423.bsp')
#planets = load('de431t.bsp')

sun = planets['sun']

from .linalg3d import *
from .constants import *

class Propagate(SpaceRock):

    def __init__(self, rocks, obsdates, model=0, gr=False,
                 gr_fast=False, gr_full=False, add_pluto=False,
                 close_encounter=False, large_asteroids=False):

        #J2=False, J4=False,

        for key in list(rocks.__dict__.keys()):
            setattr(self, key, rocks.__dict__[key])

        self.frame = rocks.frame
        self.NSIDE = rocks.NSIDE
        Propagate.gr = gr
        Propagate.gr_fast = gr_fast
        Propagate.gr_full = gr_full
        Propagate.add_pluto = add_pluto
        Propagate.CloseEncounter = close_encounter
        #Propagate.J2 = J2
        #Propagate.J4 = J4
        Propagate.LargeAsteroids = large_asteroids

        if np.isscalar(obsdates):
            obsdates = np.array([obsdates])

        Propagate.model = model
        Propagate.obsdates = obsdates
        self.propagate()


    def propagate(self):
        '''
        Integrate the bodies to the desired date. The logic could be cleaner
        but it works.
        '''

        Nx = len(Propagate.obsdates)
        Ny = len(self.x)
        x_values = np.zeros([Nx, Ny])
        y_values = np.zeros([Nx, Ny])
        z_values = np.zeros([Nx, Ny])
        vx_values = np.zeros([Nx, Ny])
        vy_values = np.zeros([Nx, Ny])
        vz_values = np.zeros([Nx, Ny])
        obsdate_values = np.zeros([Nx, Ny])
        name_values = np.tile(self.name, Nx)
        if self.H is not None:
            H_values = np.tile(self.H, Nx)

        in_frame = self.frame

        # We need to (or should) work in barycentric coordinates in rebound
        if in_frame == 'heliocentric':
            self.to_bary()

        # Rehash as a dataframe for easy access
        df = self.pandas_df()
        df['epoch'] = df['epoch'].apply(lambda idx: idx.jd)

        # Integrate all particles to the same obsdate
        pickup_times = df.epoch
        sim = self.set_simulation(np.min(pickup_times))
        sim.t = np.min(df.epoch)

        for time in np.sort(np.unique(pickup_times)):
            ps = df[df.epoch == time]
            for p in ps.itertuples():
                sim.add(x=p.x, y=p.y, z=p.z,
                        vx=p.vx, vy=p.vy, vz=p.vz,
                        hash=p.name)
                sim.integrate(time, exact_finish_time=1)

        for ii, time in enumerate(np.sort(Propagate.obsdates)):
            sim.integrate(time, exact_finish_time=1)
            for jj, name in enumerate(self.name):
                x_values[ii, jj] = sim.particles[name].x
                y_values[ii, jj] = sim.particles[name].y
                z_values[ii, jj] = sim.particles[name].z
                vx_values[ii, jj] = sim.particles[name].vx
                vy_values[ii, jj] = sim.particles[name].vy
                vz_values[ii, jj] = sim.particles[name].vz
                obsdate_values[ii, jj] = sim.t

        self.x = x_values.flatten() * u.au
        self.y = y_values.flatten() * u.au
        self.z = z_values.flatten() * u.au
        self.vx = vx_values.flatten() * (u.au / u.day)
        self.vy = vy_values.flatten() * (u.au / u.day)
        self.vz = vz_values.flatten() * (u.au / u.day)
        self.name = name_values.flatten()
        self.epoch = Time(obsdate_values.flatten(), format='jd', scale='utc')

        self.xyz_to_kep(mu_bary)

        self.t_peri = np.zeros(Nx * Ny)
        lp = self.M < np.pi * u.rad
        self.t_peri[lp] = self.epoch.jd[lp] * u.day - self.M[lp] / np.sqrt(mu_bary / self.a[lp]**3)
        self.t_peri[~lp] = self.epoch.jd[~lp] * u.day + (2*np.pi * u.rad - self.M[~lp]) / np.sqrt(mu_bary / self.a[~lp]**3)
        self.t_peri = Time(self.t_peri, format='jd', scale='utc')

        self.xyz_to_equa()
        self.varpi = (self.arg + self.node).wrap_at(2 * np.pi * u.rad)

        if self.H is not None:
            self.H = H_values
            self.mag = self.estimate_mag()

        # calculate new hpix
        if self.NSIDE is not None:
            for value in self.NSIDE:
                setattr(self, 'HPIX_{}'.format(value), self.radec_to_hpix(value))

        # be polite and return orbital parameters in the input frame.
        if in_frame == 'heliocentric':
            self.to_helio()

        return self


    def set_simulation(self, startdate):

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

        # TODO implement obliquity

        #R_sun = (696000 * u.km).to(u.au).value # JPLHorizons
        #R_mercury = (2440 * u.km).to(u.au).value # JPLHorizons
        #R_venus = (6051.893 * u.km).to(u.au).value # JPLHorizons
        #R_earth = (6378.137 * u.km).to(u.au).value # JPLHorizons
        #R_moon = (1737.4 * u.km).to(u.au).value # JPLHorizons
        #R_mars = (3396.19 * u.km).to(u.au).value # JPLHorizons
        #R_jupiter = (71492 * u.km).to(u.au).value # JPLHorizons
        #R_saturn = (60268 * u.km).to(u.au).value # JPLHorizons
        #R_uranus = (25559 * u.km).to(u.au).value # JPLHorizons
        #R_neptune = (24766 * u.km).to(u.au).value # JPLHorizons
        #R_pluto = (1188 * u.km).to(u.au).value # NASA fact sheet

        #J2_sun = -6.13e-7
        #J2_mercury = 60e-6 # Murray & Dermott
        #J2_venus = 4e-6 # Murray & Dermott
        #J2_earth = 0.00108262545 # JPLHorizons
        #J2_moon = 202.7e-6
        #J2_mars = 1960e-6 # Murray & Dermott
        #J2_jupiter = 14736e-6 # Murray & Dermott
        #J2_saturn = 16298e-6 # Murray & Dermott
        #J2_uranus = 3343e-6 # Murray & Dermott
        #J2_neptune = 3411e-6 # Murray & Dermott
        #J2_pluto = 0

        #J4_sun = 2.8e-12
        #J4_mercury = 0 # Murray & Dermott
        #J4_venus = 2e-6 # Murray & Dermott
        #J4_earth = -2e-6 # Murray & Dermott
        #J4_moon = 0
        #J4_mars = -19e-6 # Murray & Dermott
        #J4_jupiter = -587e-6 # Murray & Dermott
        #J4_saturn = -915e-6 # Murray & Dermott
        #J4_uranus = -29e-6 # Murray & Dermott
        #J4_neptune = -35e-6 # Murray & Dermott
        #J4_pluto = 0

        if Propagate.model == 0:
            active_bodies = [sun]
            names = ['Sun']
            M_sun += M_mercury + M_venus + M_earth + M_mars \
                     + M_jupiter + M_saturn + M_uranus + M_neptune
            masses = [M_sun]
            #J2 = [J2_sun]
            #J4 = [J4_sun]
            #R_eq = [R_sun]

        elif Propagate.model == 1:
            active_bodies = [sun, jupiter, saturn, uranus, neptune]
            names = ['Sun', 'Jupiter', 'Saturn', 'Uranus', 'Neptune']
            M_sun += M_mercury + M_venus + M_earth + M_mars
            masses = [M_sun, M_jupiter, M_saturn, M_uranus, M_neptune]
            #J2 = [J2_sun, J2_jupiter, J2_saturn, J2_uranus, J2_neptune]
            #J4 = [J4_sun, J4_jupiter, J4_saturn, J4_uranus, J4_neptune]
            #R_eq = [R_sun, R_jupiter, R_saturn, R_uranus, R_neptune]

        elif Propagate.model == 2:
            active_bodies = [sun, earth, jupiter, saturn, uranus, neptune]
            names = ['Sun', 'Earth', 'Jupiter', 'Saturn', 'Uranus', 'Neptune']
            M_sun += M_mercury + M_venus + M_mars
            masses = [M_sun, M_earth, M_jupiter, M_saturn, M_uranus, M_neptune]
            #J2 = [J2_sun, J2_earth, J2_jupiter, J2_saturn, J2_uranus, J2_neptune]
            #J4 = [J4_sun, J4_earth, J4_jupiter, J4_saturn, J4_uranus, J4_neptune]
            #R_eq = [R_sun, R_earth, R_jupiter, R_saturn, R_uranus, R_neptune]

        elif Propagate.model == 3:
            active_bodies = [sun, earth, mars, jupiter, saturn, uranus, neptune]
            names = ['Sun', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune']
            M_sun += M_mercury + M_venus
            masses = [M_sun, M_earth, M_mars, M_jupiter, M_saturn, M_uranus, M_neptune]
            #J2 = [J2_sun, J2_earth, J2_mars, J2_jupiter, J2_saturn, J2_uranus, J2_neptune]
            #J4 = [J4_sun, J4_earth, J4_mars, J4_jupiter, J4_saturn, J4_uranus, J4_neptune]
            #R_eq = [R_sun, R_earth, R_mars, R_jupiter, R_saturn, R_uranus, R_neptune]

        elif Propagate.model == 4:
            active_bodies = [sun, venus, earth, mars, jupiter, saturn, uranus, neptune]
            names = ['Sun', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune']
            M_sun += M_mercury
            masses = [M_sun, M_venus, M_earth, M_mars, M_jupiter, M_saturn, M_uranus, M_neptune]
            #J2 = [J2_sun, J2_venus, J2_earth, J2_mars, J2_jupiter, J2_saturn, J2_uranus, J2_neptune]
            #J4 = [J4_sun, J4_venus, J4_earth, J4_mars, J4_jupiter, J4_saturn, J4_uranus, J4_neptune]
            #R_eq = [R_sun, R_venus, R_earth, R_mars, R_jupiter, R_saturn, R_uranus, R_neptune]

        elif Propagate.model == 5:
            active_bodies = [sun, mercury, venus, earth, mars, jupiter, saturn, uranus, neptune]
            names = ['Sun', 'Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune']
            masses = [M_sun, M_mercury, M_venus, M_earth, M_mars, M_jupiter, M_saturn, M_uranus, M_neptune]
            #J2 = [J2_sun, J2_mercury, J2_venus, J2_earth, J2_mars, J2_jupiter, J2_saturn, J2_uranus, J2_neptune]
            #J4 = [J4_sun, J4_mercury, J4_venus, J4_earth, J4_mars, J4_jupiter, J4_saturn, J4_uranus, J4_neptune]
            #R_eq = [R_sun, R_mercury, R_venus, R_earth, R_mars, R_jupiter, R_saturn, R_uranus, R_neptune]

        elif Propagate.model == 6:
            active_bodies = [sun, mercury, venus, earth, moon, mars, jupiter, saturn, uranus, neptune]
            names = ['Sun', 'Mercury', 'Venus', 'Earth', 'Moon', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune']
            masses = [M_sun, M_mercury, M_venus, M_earth, M_moon, M_mars, M_jupiter, M_saturn, M_uranus, M_neptune]
            #J2 = [J2_sun, J2_mercury, J2_venus, J2_earth, J2_moon, J2_mars, J2_jupiter, J2_saturn, J2_uranus, J2_neptune]
            #J4 = [J4_sun, J4_mercury, J4_venus, J4_earth, J4_moon, J4_mars, J4_jupiter, J4_saturn, J4_uranus, J4_neptune]
            #R_eq = [R_sun, R_mercury, R_venus, R_earth, R_moon, R_mars, R_jupiter, R_saturn, R_uranus, R_neptune]

        else:
            raise ValueError('Model not recognized. Check the documentation.')

        if Propagate.add_pluto == True:
            pluto = planets['pluto barycenter']
            active_bodies.append(pluto)
            names.append('Pluto')
            masses.append(M_pluto)
            #J2.append(J2_pluto)
            #J4.append(J4_pluto)
            #R_eq.append(R_pluto)

        startdate = Time(startdate, scale='utc', format='jd')
        t = ts.ut1(jd=startdate.ut1.jd)

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
        ss['a'] = 1 / (2 / norm([ss.x, ss.y, ss.z]) - norm([ss.vx, ss.vy, ss.vz])**2 / mu_bary.value)
        ss['hill_radius'] = ss.a * pow(ss.mass / (3 * M_sun), 1/3)
        ss['name'] = names
        #ss['J2'] = J2
        #ss['J4'] = J4
        #ss['R_eq'] = R_eq

        sim = rebound.Simulation()
        sim.units = ('day', 'AU', 'Msun')

        for p in ss.itertuples():
            sim.add(x=p.x, y=p.y, z=p.z,
                    vx=p.vx, vy=p.vy, vz=p.vz,
                    m=p.mass, hash=p.name, r=p.hill_radius)


        #if Propagate.J2 == True:

        #    rebx = reboundx.Extras(sim)
        #    harmonics = rebx.load_force("gravitational_harmonics")
        #    rebx.add_force(harmonics)
        #    bodies = sim.particles

        #    for body in ss.name:
        #        bodies[body].params['J2'] = ss[ss.name == body].J2
        #        bodies[body].params['R_eq'] = ss[ss.name == body].R_eq
        #        if Propagate.J4 == True:
        #            bodies[body].params['J4'] = ss[ss.name == body].J4

        sim.N_active = len(ss)

        if Propagate.LargeAsteroids == True:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                sim.add('Ceres', date='JD{}'.format(startdate))
                sim.particles[-1].hash = 'Ceres'
                sim.particles['Ceres'].m = 4.504e-10

                sim.add('Vesta', date='JD{}'.format(startdate))
                sim.particles[-1].hash = 'Vesta'
                sim.particles['Vesta'].m = 1.302e-10

                sim.add('Pallas', date='JD{}'.format(startdate))
                sim.particles[-1].hash = 'Pallas'
                sim.particles['Pallas'].m = 1.06e-10

                sim.N_active += 3

        if Propagate.gr_fast == True:
            bodies = sim.particles
            rebx = reboundx.Extras(sim)
            gr = rebx.load_force('gr_potential')
            gr.params['c'] = constants.C
            rebx.add_force(gr)
            bodies['Sun'].params['gr_source'] = 1

        elif Propagate.gr == True:
            bodies = sim.particles
            rebx = reboundx.Extras(sim)
            gr = rebx.load_force('gr')
            gr.params['c'] = constants.C
            rebx.add_force(gr)
            bodies['Sun'].params['gr_source'] = 1

        elif Propagate.gr_full == True:
            rebx = reboundx.Extras(sim)
            gr = rebx.load_force('gr_full')
            gr.params["c"] = constants.C
            rebx.add_force(gr)

        sim.testparticle_type = 0

        if Propagate.CloseEncounter == True:
            sim.integrator = 'ias15'
            sim.ri_ias15.epsilon = 1e-10
            rebx = reboundx.Extras(sim)
            harmonics = rebx.load_force("gravitational_harmonics")
            rebx.add_force(harmonics)

            sim.particles['Earth'].params['J2'] = 0.00108262545 # JPLHorizons
            sim.particles['Earth'].params['R_eq'] = (6378.137 * u.km).to(u.au).value # JPLHorizons


        else:
            sim.integrator = 'mercurius'
            sim.dt = 1 # one day
            sim.ri_ias15.min_dt = sim.dt / 1440 # one minute
            sim.ri_mercurius.hillfac = 3

        sim.move_to_com()

        return sim
