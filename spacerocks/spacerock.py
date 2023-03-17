import uuid

from astropy import units as u
from astropy.coordinates import Angle, Distance
from astropy.time import Time

import datetime

import numpy as np
import pandas as pd

import asdf
import copy

from .constants import mu_bary, c, epsilon
from .keplerorbit import KeplerOrbit
from .convenience import Convenience
from .utils import time_handler
from .units import Units
from .vector import Vector
from .ephemerides import Ephemerides 
from .observer import Observer
from .cbindings import kepM_to_xyz, correct_for_ltt, correct_for_ltt_single_observer, kepf_to_xyz
from .spice import SpiceBody
from .paths import OBSERVATORIES_PATH

observatories = pd.read_csv(OBSERVATORIES_PATH)


class SpaceRock(KeplerOrbit, Convenience):
    '''
    SpaceRock objects provide an interface to work with Keplerian orbits. 
    '''
    def __init__(self, origin='ssb', frame='eclipJ2000', units=Units(), *args, **kwargs):

        coords = self.detect_coords(kwargs)
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

        self.units = units
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
                self.t_peri = time_handler(kwargs.get('t_peri'), units)

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

            if kwargs.get('name') is not None:
                self.name = np.atleast_1d(kwargs.get('name'))
            else:
                # produces random, non-repeting integers between 0 and 1e10 - 1
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

            if kwargs.get('name') is not None:
                self.name = np.atleast_1d(kwargs.get('name'))
            else:
                # produces random, non-repeting identifiers
                self.name = np.array([uuid.uuid4().hex for _ in range(len(self.x))])

        if kwargs.get('epoch') is not None:
            self.epoch = time_handler(kwargs.get('epoch'), units)

        # brightness information

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

        # Physical properties

        if kwargs.get('radius') is not None:
            self.radius = Distance(kwargs.get('radius'), units.size, allow_negative=False)

        if kwargs.get('mass') is not None:
            self.mass = kwargs.get('mass') * units.mass

        if kwargs.get('density') is not None:
            self.density = kwargs.get('density') * units.density

        # Uncertainties
        if kwargs.get('cov') is not None:
            self.cov = kwargs.get('cov')

    def clone(self, N: int):
        epochs = np.repeat(self.epoch.tdb.jd, N)
        names = [self.name[0] + f'-clone{idx}' for idx in range(N)]

        mean = [self.x.au[0], 
                self.y.au[0], 
                self.z.au[0], 
                self.vx.value[0], 
                self.vy.value[0], 
                self.vz.value[0]]
        xs, ys, zs, vxs, vys, vzs = np.random.multivariate_normal(mean, self.cov[0], N).T

        return SpaceRock(x=xs, y=ys, z=zs, vx=vxs, vy=vys, vz=vzs, name=names,
                         origin=self.origin, frame=self.frame, units=self.units, epoch=epochs)

    # @classmethod
    # def from_spice(cls, spiceid: str, epoch=datetime.date.today(), kernel=SpiceKernel(), units=Units()):
    #     body = SpiceBody(spiceid=spiceid, kernel=kernel)
    #     if isinstance(epoch, datetime.date):
    #         epoch = Time(datetime.date.today().isoformat(), format='iso', scale='utc')
    #     else:
    #         epoch = self.detect_timescale(epoch, units.timescale)
    #     r = body.at(epoch)
    #     rock = SpaceRock(x=r.x.au[0], 
    #                      y=r.y.au[0], 
    #                      z=r.z.au[0], 
    #                      vx=r.vx.value[0], 
    #                      vy=r.vy.value[0], 
    #                      vz=r.vz.value[0], 
    #                      epoch=epoch, 
    #                      origin='ssb', 
    #                      mass=r.mass,
    #                      frame='ECLIPJ2000')
    #     return rock

    @classmethod
    def from_mpc(cls, catalog: str, download_data=False, metadata='Orbit_type'):

        metadata = np.atleast_1d(metadata)

        from .paths import MPC_PATH
        import pathlib
        pathlib.Path(MPC_PATH).mkdir(exist_ok=True)
        path = pathlib.Path(MPC_PATH + f'/{catalog}.feather')
        
        if path.exists() and download_data == False:
            df = pd.read_feather(path)
        else: 
            from .downloader import download
            datafile = ['https://minorplanetcenter.net/Extended_Files/' + catalog + '.json.gz']
            download(datafile, MPC_PATH)
            path = MPC_PATH + f'/{catalog}.json.gz'
            df = pd.read_json(path)
            df.to_feather(f'{MPC_PATH}/{catalog}.feather')

         

        units = Units()
        units.timescale = 'tt'
        units.timeformat = 'jd'

        

        H = df.H
        G = df.G
        name = df.Principal_desig
        epoch = df.Epoch
        M = df.M
        arg = df.Peri
        node = df.Node
        inc = df.i
        e = df.e
        a = df.a

        rocks = cls(a=a, e=e, inc=inc, arg=arg, node=node, M=M, 
                    epoch=epoch, name=name, H=H, G=G, units=units, origin='sun')

        for data_item in metadata:
            setattr(rocks, data_item, df[data_item].values)        

        return rocks

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
        covs = []

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
                      'OUT_UNITS':   quote("AU-D"),
                      'CSV_FORMAT':  quote("YES"),
                      'VEC_DELTA_T': quote("NO"),
                      'VEC_TABLE':   quote("2x"),
                      'VEC_LABELS':  quote("NO")
                     }
    
            try:
                url = "https://ssd.jpl.nasa.gov/api/horizons.api?" + urlencode(params)
                with urlopen(url) as f:
                    body = f.read().decode()

                lines = body.split("$$SOE")[-1].split("\n")
                #spl = [float(i) for i in lines[2].split(',')[2:-1]]

                spl_1 = [float(i) for i in lines[2].split(',')[2:-7]]
                x, y, z, vx, vy, vz = spl_1
                
                try:
                    spl_2 = [float(i) for i in lines[2].split(',')[8:-1]]
                    dx, dy, dz, dvx, dvy, dvz = spl_2
                except ValueError as e:
                    dx, dy, dz, dvx, dvy, dvz = 0, 0, 0, 0, 0, 0
                
                epoch = float(lines[2].split(',')[0])
                #x, y, z, vx, vy, vz = spl[:6]
                #dx, dy, dz, dvx, dvy, dvz = spl[6:]

                cov = np.diag([dx**2, dy**2, dz**2, dvx**2, dvy**2, dvz**2])
        
                try:
                    pps = body.split('physical parameters')[-1].split('\n')[2].split()
                    H = float(pps[1])
                    G = float(pps[3])
                except Exception as E:
                    print(f'Magnitude information not available for {n}. Is this a comet? A planet? A moon?')
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
                epochs.append(epoch)
                names.append(n)
                covs.append(cov)

            except Exception as E:
                #print(E)
                #raise ValueError(f'Body {n} not found.')
                raise ValueError(f'{E}')
            
        units = Units()
        units.distance = u.au
        units.speed = u.au/u.day
        units.timescale = 'tdb'

        return cls(x=xs, y=ys, z=zs, vx=vxs, vy=vys, vz=vzs, name=names,
                   origin='sun', frame='eclipJ2000', units=units, H=Hs, G=Gs, epoch=epochs, cov=covs)
        
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

    def propagate(self, epochs, model='GIANTS', units=Units(), gr=False, progress=True, in_memory=True, outfile=None, **kwargs):
        from .simulation import Simulation
        epochs = time_handler(epochs, units)
        origin = copy.copy(self.origin)
        frame = copy.copy(self.frame)

        # We need to integrate in barycentric coordinates
        self.to_bary()

        # Put the simulation into ecliptic coordinates
        self.change_frame('eclipJ2000')

        units = Units()
        units.timescale = 'tdb'
        units.timeformat = 'jd'
        sim = Simulation(model=model, epoch=self.epoch.tdb.jd.min(), units=units)
        sim.add_spacerocks(self)
        sim.N_active = len(sim.model.perturbers)
        sim.testparticle_type = 0
        sim.integrator = 'ias15'

        if gr == True:
            import reboundx
            rebx = reboundx.Extras(sim)
            gr_force = rebx.load_force('gr_full')
            gr_force.params["c"] = c.to(u.au / u.day).value
            rebx.add_force(gr_force)

        units = Units()
        units.timescale = 'tdb'
        units.timeformat = 'jd'
        prop, planets, sim = sim.propagate(epochs=epochs.tdb.jd, units=units, progress=progress, 
                                           in_memory=in_memory, outfile=outfile, **kwargs)

        # be polite and return orbital parameters using the input origin.
        if origin != 'ssb':
            self.change_origin(origin)
           
        if frame != 'eclipJ2000':
            self.change_frame(frame)

        if in_memory == True:

            if frame != 'eclipJ2000':
                prop.change_frame(frame)
                planets.change_frame(frame)

            N = len(epochs)
            if hasattr(self, 'G'):
                prop.G = np.repeat(self.G, N)
    
            if hasattr(self, 'H_func'):
                prop.H_func = np.repeat(self.H_func, N)
            elif hasattr(self, '_H'):
                prop.H = np.repeat(self.H, N)
    
            if hasattr(self, 'mag_func'):
                prop.mag_func = np.repeat(self.mag_func, N)

        return prop, planets, sim


    def calc_H(self, obscode):
        return self.__calc_H_from_mag(obscode=obscode)

    def __calc_H_from_mag(self, obscode):
        obs = self.observe(obscode=obscode)

        earth = SpiceBody(spiceid='Earth')
        sun = SpiceBody(spiceid='Sun')

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

        The James Webb Space Telescope will be supported as soon as
        NASA provides the necessary spk files.
        '''
        if kwargs.get('obscode') is not None:
            x, y, z, vx, vy, vz = self.xyz_to_tel(obscode=kwargs.get('obscode'))
        elif kwargs.get('spiceid') is not None:
            x, y, z, vx, vy, vz = self.xyz_to_tel(spiceid=kwargs.get('spiceid'))
        elif kwargs.get('observer') is not None:
            x, y, z, vx, vy, vz = self.xyz_to_tel(observer=kwargs.get('observer'))
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
            # observer = Observer(obscode=kwargs.get('obscode'), epoch=self.epoch.utc.jd)
            observer = Observer.from_obscode(obscode=kwargs.get('obscode'))
            observer = observer.at(self.epoch.utc.jd)
            observer.change_frame('eclipJ2000')
        # elif kwargs.get('spiceid') is not None:
        #     observer = Observer(spiceid=kwargs.get('spiceid'), epoch=self.epoch.utc.jd)
        elif kwargs.get('observer') is not None:
            observer = kwargs.get('observer')
        else:
            raise ValueError('Must pass either an obscode, spiceid, or Observer object.')

        in_origin = copy.copy(self.origin)
        self.to_bary()

        in_frame = copy.copy(self.frame)
        self.change_frame('eclipJ2000') 

        if len(observer.x) == 1:
            dx, dy, dz, dvx, dvy, dvz = correct_for_ltt_single_observer(self, observer)
        else:
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

    def analytic_propagate(self, epoch, units=Units()):
        '''
        analytically propagate all rocks to a common epoch
        '''
        epoch = time_handler(epoch, units)
        #if not isinstance(epoch, Time):
        #    epoch = self.detect_timescale(np.atleast_1d(epoch), units.timescale)
        dt = (epoch.utc.jd - self.epoch.utc.jd) * u.day
        dM = self.n * dt

        new_units = copy.deepcopy(self.units)
        new_units.distance = u.au
        new_units.angle = u.deg
        if hasattr(self, '_mass'):
            rock =  self.__class__(a=self.a.au, 
                                  e=self.e, 
                                  inc=self.inc.deg, 
                                  node=self.node.deg, 
                                  arg=self.arg.deg, 
                                  M=(self.M + dM).deg, 
                                  name=self.name, 
                                  mass=self.mass.value,
                                  epoch=np.repeat(epoch.utc.jd, len(self)), 
                                  origin=self.origin, 
                                  frame=self.frame, 
                                  units=new_units)
        else:
            rock = self.__class__(a=self.a.au, 
                                  e=self.e, 
                                  inc=self.inc.deg, 
                                  node=self.node.deg, 
                                  arg=self.arg.deg, 
                                  M=(self.M + dM).deg, 
                                  name=self.name, 
                                  epoch=np.repeat(epoch.utc.jd, len(self)), 
                                  origin=self.origin, 
                                  frame=self.frame, 
                                  units=new_units)

        N = 1
        if hasattr(self, 'G'):
            rock.G = np.tile(self.G, N)

        if hasattr(self, 'H_func'):
            rock.H_func = np.tile(self.H_func, N)
        elif hasattr(self, '_H'):
            rock.H = np.tile(self.H, N)

        if hasattr(self, 'mag_func'):
            rock.mag_func = np.tile(self.mag_func, N)

        return rock
        

    def orbits(self, N=1000):

        f = Angle(np.linspace(0, 2*np.pi, N), u.rad)

        # xs = []
        # ys = []
        # zs = []

        vectors = []

        for r in self:
            x, y, z, _, _, _ = kepf_to_xyz(np.repeat(r.a, N),
                                           np.repeat(r.e, N),
                                           np.repeat(r.inc.rad, N),
                                           np.repeat(r.arg.rad, N),
                                           np.repeat(r.node.rad, N),
                                           f.rad)
            vectors.append(Vector(x, y, z))
            # xs.append(x)
            # ys.append(y)
            # zs.append(z)

        return vectors

    def write_to(self, path, compression='zlib'):
        uniquenames = np.unique(self.name)
        if len(uniquenames) == 1:
            name = uniquenames.tolist()
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

    def groupby(self, parameter):
        return ([x, self[getattr(self, parameter) == x]] for x in sorted(set(getattr(self, parameter))))
