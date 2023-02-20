from astropy import units as u
from astropy.coordinates import Distance, Angle

from .constants import epsilon
from .spice import SpiceKernel
from .vector import Vector

import numpy as np
import pandas as pd
import spiceypy as spice
import copy

from .paths import OBSERVATORIES_PATH
observatories = pd.read_csv(OBSERVATORIES_PATH)


EQUAT_RAD = 6378137
FLATTEN = 1 / 298.257223563
O_M_FLATTEN = 1 - FLATTEN
DEG_TO_RAD = np.pi / 180


class Observer:

    def __init__(self, origin='ssb', frame='ECLIPJ2000', kernel=SpiceKernel(), **kwargs):
        
        self.origin = origin
        self.frame = frame
        
        if kwargs.get('epoch') is not None:
            self.epoch = np.atleast_1d(kwargs.get('epoch'))

        if kwargs.get('obscode') is not None:
            self.obscode = np.atleast_1d(kwargs.get('obscode'))
            self.spiceid = np.atleast_1d('Earth')

        elif kwargs.get('spiceid') is not None:
            self.spiceid = np.atleast_1d(kwargs.get('spiceid'))

        else:
            raise ValueError('Must specify either a spiceid or an obscode')

        self.__get_all_state_vectors()
        
    def __len__(self):
        '''
        This method allows you to use the len() function on a SpaceRocks object.
        '''
        return len(self.epoch)

    def __getitem__(self, idx):
        '''
        This method allows you to index a SpaceRocks object.
        '''
        p = copy.copy(self)
        for attr in self.__dict__.keys():
            if (attr != 'mu') and (attr != 'frame') and (attr != 'origin') and (attr != 'units'):
                if isinstance(getattr(self, attr), Vector):
                    setattr(p, attr, getattr(self, attr)[idx])
                else:
                    setattr(p, attr, getattr(self, attr)[idx])

        return p

    def __get_all_state_vectors(self):
        '''
        Very optimized way to get all state vectors from spice.

        The code is sort of convoluted, but it works and is several 
        times faster than any alternative without resorting to 
        using CSpice directly (don't want to deal with that).

        1. set(list(zip(...))) gets unique (body, time) pairs.

        2. These pairs are used as the keys in a dictionary storing
           the state vectors as numpy arrays. The key is also passed 
           to the __state_from_spice function to calculate the state 
           vector.

        3. Use a list comprehension to populate a numpy array with 
           all state vectors by reading from the dictionary. This is faster 
           than numpy indexing. The time complexity is O(1).

        4. Set attributes using the usual astropy units.

        TODO: Allow for mixing between terrestrial and non-terrestrial observers 
        '''

        arr = np.empty((len(self.epoch), 6))
        
        if len(self.spiceid) == 1:
            unique_times = np.unique(self.epoch)
            sid = self.spiceid[0]
            unique_dict = {key:self.__state_from_spice((sid, key)) for key in unique_times}
            for idx, k in enumerate(self.epoch):
                arr[idx] = unique_dict[k]
            # for key, value in unique_dict.items():
            #     arr[self.epoch == key] = value
        else:
            unique_rocks_and_times = set(list(zip(self.spiceid, self.epoch)))
            unique_dict = {key:self.__state_from_spice(key) for key in unique_rocks_and_times}
            for key, value in unique_dict.items():
                arr[(self.spiceid == key[0]) * (self.epoch == key[1])] = value

        x, y, z, vx, vy, vz = arr.T 

        if hasattr(self, 'obscode'):
            
            icrs_arr = np.empty((len(self.epoch), 3))
            
            if len(self.obscode) == 1:
                lat, lon, elevation = self.__get_observatory_params(self.obscode[0])
                unique_obs_dict = {key:self.__compute_topocentric_correction((lon, lat, elevation, key)) for key in unique_times}
                for idx, k in enumerate(self.epoch):
                    icrs_arr[idx] = unique_obs_dict[k]
                
            else:
                unique_codes = set(self.obscode)
                unique_codes_dict = {key:self.__get_observatory_params(key) for key in unique_codes}

                lats, lons, elevations = np.array([unique_codes_dict[obscode] for obscode in self.obscode]).T

                self.lat = Angle(lats, u.deg)
                self.lon = Angle(lons, u.deg)
                self.elevation = elevations
            
                unique_locations_and_times = set(list(zip(self.lon.deg, self.lat.deg, self.elevation, self.epoch)))
                unique_obs_dict = {key:self.__compute_topocentric_correction(key) for key in unique_locations_and_times}
                for key, value in unique_dict.items():
                    icrs_arr[(self.lon == key[0]) * (self.lat == key[1]) * (self.elevation == key[2]) * (self.epoch == key[3])] = value
            

            dx_icrs, dy_icrs, dz_icrs = icrs_arr.T
            # Transform from icrs to ecliptic
            dx = dx_icrs
            dy = dy_icrs * np.cos(epsilon) + dz_icrs * np.sin(epsilon)
            dz = -dy_icrs * np.sin(epsilon) + dz_icrs * np.cos(epsilon)

            x += (dx * u.m).to(u.km).value
            y += (dy * u.m).to(u.km).value
            z += (dz * u.m).to(u.km).value
       
        self.x = Distance(x, u.km, allow_negative=True).to(u.au)
        self.y = Distance(y, u.km, allow_negative=True).to(u.au)
        self.z = Distance(z, u.km, allow_negative=True).to(u.au)

        self.vx = (vx * u.km/u.s).to(u.au / u.day)
        self.vy = (vy * u.km/u.s).to(u.au / u.day)
        self.vz = (vz * u.km/u.s).to(u.au / u.day)

    def __get_observatory_params(self, obscode):
        obs = observatories[observatories.obscode == obscode]
        return np.array([obs.lat.values[0], obs.lon.values[0], obs.elevation.values[0]])

    def __compute_topocentric_correction(self, x):
        observer_lon, observer_lat, observer_elevation, epoch = x
        
        observer_lat *= DEG_TO_RAD
        observer_lon *= DEG_TO_RAD

        sin_lat = np.sin(observer_lat)
        cos_lat = np.cos(observer_lat)

        lon = self.__compute_local_sidereal_time(epoch, observer_lon)
        
        sin_lon = np.sin(lon)
        cos_lon = np.cos(lon)

        denom = O_M_FLATTEN * sin_lat
        denom = cos_lat * cos_lat + denom*denom

        C_geo = 1 / np.sqrt(denom)
        S_geo = O_M_FLATTEN * O_M_FLATTEN * C_geo
        C_geo = C_geo * EQUAT_RAD + observer_elevation
        S_geo = S_geo * EQUAT_RAD + observer_elevation
        dx = C_geo * cos_lat * cos_lon
        dy = C_geo * cos_lat * sin_lon
        dz = S_geo * sin_lat

        return dx, dy, dz

    def __compute_local_sidereal_time(self, epoch, lon):
        T = (epoch - 2451545.0) / 36525
        theta = 280.46061837 + 360.98564736629 * (epoch - 2451545.0) + (0.000387933 * T * T) - (T * T * T / 38710000.0)
        theta *= DEG_TO_RAD
        return theta + lon
    
    def __compute_ephemeris_time(self, epoch):
        '''
        Wrapper for spiceypy's str2et function.
        '''
        return spice.str2et('JD{} UTC'.format(epoch))

    def __state_from_spice(self, x):
        spiceid, epoch = x
        ephemeris_time = self.__compute_ephemeris_time(epoch)
        state, _ = spice.spkezr(spiceid, ephemeris_time, self.frame, 'none', self.origin)
        return state
