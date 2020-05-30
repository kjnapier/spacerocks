from spacerocks import SpaceRock
import numpy as np
import pandas as pd
import rebound
import reboundx
from reboundx import constants

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
sun = planets['sun']

from .linalg3d import *
from .constants import *

# Return the planet positions at each time.

class Propagate(SpaceRock):

    def __init__(self, rocks, obsdates, model=0, gr=False,
                 gr_fast=False, gr_full=False, add_pluto=False):

        for key in list(rocks.__dict__.keys()):
            setattr(self, key, rocks.__dict__[key])

        self.frame = rocks.frame
        self.NSIDE = rocks.NSIDE
        Propagate.gr = gr
        Propagate.gr_fast = gr_fast
        Propagate.gr_full = gr_full
        Propagate.add_pluto = add_pluto

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
        df['obsdate'] = df['obsdate'].apply(lambda idx: idx.jd)

        # Integrate all particles to the same obsdate
        pickup_times = df.obsdate
        sim = self.set_simulation(np.min(pickup_times))
        sim.t = np.min(df.obsdate)

        for time in np.sort(np.unique(pickup_times)):
            ps = df[df.obsdate == time]
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
        self.obsdate = Time(obsdate_values.flatten(), format='jd', scale='utc')

        self.xyz_to_kep(mu_bary)

        self.t_peri = np.zeros(Nx * Ny)
        lp = self.M < np.pi * u.rad
        self.t_peri[lp] = self.obsdate.jd[lp] * u.day - self.M[lp] / np.sqrt(mu_bary / self.a[lp]**3)
        self.t_peri[~lp] = self.obsdate.jd[~lp] * u.day + (2*np.pi * u.rad - self.M[~lp]) / np.sqrt(mu_bary / self.a[~lp]**3)
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

        mercury = planets['mercury']
        venus = planets['venus']
        earth = planets['earth']
        mars = planets['mars']
        jupiter = planets['jupiter barycenter']
        saturn = planets['saturn barycenter']
        uranus = planets['uranus barycenter']
        neptune = planets['neptune barycenter']

        M_mercury = 1.6601367952719304e-7
        M_venus = 2.4478383396645447e-6
        M_earth = 3.040432648022642e-6 # Earth-Moon Barycenter
        M_mars = 3.2271560375549977e-7 # Mars Barycenter
        M_sun = 1
        M_jupiter = 9.547919384243222e-4
        M_saturn = 2.858859806661029e-4
        M_uranus = 4.3662440433515637e-5
        M_neptune = 5.151389020535497e-5
        M_pluto = 7.361781606089469e-9

        if Propagate.model == 0:
            active_bodies = [sun]
            names = ['Sun']
            M_sun += M_mercury + M_venus + M_earth + M_mars \
                     + M_jupiter + M_saturn + M_uranus + M_neptune
            masses = [M_sun]

        elif Propagate.model == 1:
            active_bodies = [sun, jupiter, saturn, uranus, neptune]
            names = ['Sun', 'Jupiter', 'Saturn', 'Uranus', 'Neptune']
            M_sun += M_mercury + M_venus + M_earth + M_mars
            masses = [M_sun, M_jupiter, M_saturn, M_uranus, M_neptune]

        elif Propagate.model == 2:
            active_bodies = [sun, earth, jupiter, saturn, uranus, neptune]
            names = ['Sun', 'Earth', 'Jupiter', 'Saturn', 'Uranus', 'Neptune']
            M_sun += M_mercury + M_venus + M_mars
            masses = [M_sun, M_earth, M_jupiter, M_saturn, M_uranus, M_neptune]

        elif Propagate.model == 3:
            active_bodies = [sun, earth, mars, jupiter, saturn, uranus, neptune]
            names = ['Sun', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune']
            M_sun += M_mercury + M_venus
            masses = [M_sun, M_earth, M_mars,
                      M_jupiter, M_saturn, M_uranus, M_neptune]

        elif Propagate.model == 4:
            active_bodies = [sun, venus, earth, mars, jupiter, saturn, uranus, neptune]
            names = ['Sun', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune']
            M_sun += M_mercury
            masses = [M_sun, M_venus, M_earth, M_mars,
                      M_jupiter, M_saturn, M_uranus, M_neptune]

        elif Propagate.model == 5:
            active_bodies = [sun, mercury, venus, earth, mars, jupiter, saturn, uranus, neptune]
            names = ['Sun', 'Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune']
            masses = [M_sun, M_mercury, M_venus, M_earth, M_mars, M_jupiter, M_saturn, M_uranus, M_neptune]

        if Propagate.add_pluto == True:
            pluto = planets['pluto barycenter']
            active_bodies.append(pluto)
            names.append('Pluto')
            masses.append(M_pluto)

        t = ts.tt(jd=startdate)

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

        sim = rebound.Simulation()
        sim.units = ('day', 'AU', 'Msun')

        for p in ss.itertuples():
            sim.add(x=p.x, y=p.y, z=p.z,
                    vx=p.vx, vy=p.vy, vz=p.vz,
                    m=p.mass, hash=p.name, r=p.hill_radius)

        sim.move_to_com()

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

        sim.N_active = len(ss)
        sim.testparticle_type = 0

        sim.integrator = 'mercurius'
        sim.dt = 1 # one day

        sim.ri_ias15.min_dt = sim.dt / 1440 # one minute
        sim.ri_mercurius.hillfac = 3

        return sim
