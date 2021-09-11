from astropy import units as u
from astropy.coordinates import Distance, Angle

from .constants import epsilon

import numpy as np
from numpy import sin, cos, sqrt
import pandas as pd
import spiceypy as spice

import os
import pkg_resources

DATA_PATH = pkg_resources.resource_filename('spacerocks', 'data/observatories.csv')
observatories = pd.read_csv(DATA_PATH)

SPICE_PATH = pkg_resources.resource_filename('spacerocks', 'data/spice')
spice.furnsh(os.path.join(SPICE_PATH, 'latest_leapseconds.tls'))
spice.furnsh(os.path.join(SPICE_PATH, 'de440s.bsp'))
spice.furnsh(os.path.join(SPICE_PATH, 'hst.bsp'))
spice.furnsh(os.path.join(SPICE_PATH, 'nh.bsp'))


class Observer:

    def __init__(self, origin='ssb', frame='ECLIPJ2000', **kwargs):
        
        self.origin = origin
        self.frame = frame

        if kwargs.get('epoch') is not None:
            self.epoch = np.atleast_1d(kwargs.get('epoch'))

        if kwargs.get('obscode') is not None:
            if len(kwargs.get('obscode')) < len(self.epoch):
                self.obscode = np.repeat(kwargs.get('obscode'), len(self.epoch)) 
            else:
                self.obscode = np.atleast_1d(kwargs.get('obscode'))
            self.spiceid = np.repeat('Earth', len(self.obscode))
            
            unique_codes = set(self.obscode)
            unique_codes_dict = {key:self.__get_observatory_params(key) for key in unique_codes}

            lats, lons, elevations = np.array([unique_codes_dict[obscode] for obscode in self.obscode]).T

            self.lat = Angle(lats, u.deg)
            self.lon = Angle(lons, u.deg)
            self.elevation = elevations

        elif kwargs.get('spiceid') is not None:
            if len(kwargs.get('spiceid')) < len(self.epoch):
                self.spiceid = np.repeat(kwargs.get('spiceid'), len(self.epoch))
            else:
                self.spiceid = np.atleast_1d(kwargs.get('spiceid'))

        else:
            raise ValueError('Must specify either a spiceid or an obscode')

        self.__get_all_state_vectors()


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

        unique_rocks_and_times = set(list(zip(self.spiceid, self.epoch)))
        unique_dict = {key:self.__state_from_spice(key) for key in unique_rocks_and_times}
        x, y, z, vx, vy, vz = np.array([unique_dict[(body, time)] for body, time in zip(self.spiceid, self.epoch)]).T

        if hasattr(self, 'obscode'):
            unique_locations_and_times = set(list(zip(self.lon.deg, self.lat.deg, self.elevation, self.epoch)))
            unique_obs_dict = {key:self.__compute_topocentric_correction(key) for key in unique_locations_and_times}
            dx_icrs, dy_icrs, dz_icrs = np.array([unique_obs_dict[(lon, lat, elev, time)] for lon, lat, elev, time in zip(self.lon.deg, self.lat.deg, self.elevation, self.epoch)]).T

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

        observer_lat = Angle(observer_lat, u.deg)
        observer_lon = Angle(observer_lon, u.deg)

        EQUAT_RAD = 6378137
        FLATTEN = 1 / 298.257223563

        lon = self.__compute_local_sidereal_time(epoch, observer_lon)
        
        denom = (1 - FLATTEN) * sin(observer_lat)
        denom = cos(observer_lat) * cos(observer_lat) + denom*denom

        C_geo = 1 / sqrt(denom)
        S_geo = (1 - FLATTEN) * (1 - FLATTEN) * C_geo
        C_geo = C_geo * EQUAT_RAD + observer_elevation
        S_geo = S_geo * EQUAT_RAD + observer_elevation
        dx = C_geo * cos(observer_lat) * cos(lon)
        dy = C_geo * cos(observer_lat) * sin(lon)
        dz = S_geo * sin(observer_lat)
        return dx, dy, dz

    def __compute_local_sidereal_time(self, epoch, lon):
        T = (epoch - 2451545.0) / 36525
        theta = Angle(280.46061837 + 360.98564736629 * (epoch - 2451545.0) + (0.000387933 * T * T) - (T * T * T / 38710000.0), u.deg)
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