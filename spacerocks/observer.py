from .paths import OBSERVATORIES_PATH
from .units import Units
from .spice import SpiceKernel

import pandas as pd
import spiceypy as spice
import numpy as np
from astropy import units as u

import copy

observatories = pd.read_csv(OBSERVATORIES_PATH)
EQUAT_RAD = 6378137
FLATTEN = 1 / 298.257223563
O_M_FLATTEN = 1 - FLATTEN
DEG_TO_RAD = np.pi / 180

def get_observatory_params(obscode):
    obs = observatories[observatories.obscode == obscode]
    return np.array([obs.lat.values[0], obs.lon.values[0], obs.elevation.values[0]])

def compute_local_sidereal_time(epoch, lon):
    T = (epoch - 2451545.0) / 36525
    theta = 280.46061837 + 360.98564736629 * (epoch - 2451545.0) + (0.000387933 * T * T) - (T * T * T / 38710000.0)
    theta *= DEG_TO_RAD
    return theta + lon

def sidereal_rate(epoch):
    T = (epoch - 2451545.0) / 36525
    Tprime = 1 / 36525
    theta_dot = 360.98564736629 + 2 * 0.000387933 * T * Tprime + 3 * T * T * Tprime / 38710000.0
    return theta_dot * DEG_TO_RAD

def compute_topocentric_correction(lon, lat, elevation, epoch):
        
    observer_lat = lat * DEG_TO_RAD
    observer_lon = lon * DEG_TO_RAD
    sin_lat = np.sin(observer_lat)
    cos_lat = np.cos(observer_lat)

    phi = compute_local_sidereal_time(epoch, observer_lon)
    phi_prime = sidereal_rate(epoch)
    
    sin_lon = np.sin(phi)
    cos_lon = np.cos(phi)

    sin_lon_prime = cos_lon * phi_prime
    cos_lon_prime = -sin_lon * phi_prime
        
    denom = O_M_FLATTEN * sin_lat
    denom = cos_lat * cos_lat + denom*denom

    C_geo = 1 / np.sqrt(denom)
    S_geo = O_M_FLATTEN * O_M_FLATTEN * C_geo
    C_geo = C_geo * EQUAT_RAD + elevation
    S_geo = S_geo * EQUAT_RAD + elevation
    dx = C_geo * cos_lat * cos_lon
    dy = C_geo * cos_lat * sin_lon
    dz = S_geo * sin_lat

    # C_geo_prime = S_geo_prime = 0

    dvx = cos_lat * (cos_lon_prime * C_geo)
    dvy = cos_lat * (sin_lon_prime * C_geo)
    dvz = 0
    return dx, dy, dz, dvx, dvy, dvz

# function to get the earth's position at any given epoch with caching
def get_earth_state(epoch):
    if epoch not in get_earth_state.cache:
        et = spice.str2et('JD{} UTC'.format(epoch))
        state, _ = spice.spkezr('Earth', et, 'J2000', 'none', 'ssb')
        get_earth_state.cache[epoch] = state
    return get_earth_state.cache[epoch]
get_earth_state.cache = {}


class Observer:

    def __init__(self, lat, lon, elevation, kernel=SpiceKernel()):
        self.lat = lat
        self.lon = lon
        self.elevation = elevation

    @classmethod
    def from_obscode(cls, obscode):
        obscode = np.atleast_1d(obscode)
        unique_codes = set(obscode)

        unique_codes_dict = {key:get_observatory_params(key) for key in unique_codes}
        lats, lons, elevations = np.array([unique_codes_dict[obscode] for obscode in obscode]).T

        lat = lats
        lon = lons
        elevation = elevations
        return cls(lat, lon, elevation)

    @classmethod
    def from_coordinates(cls, lat, lon, elevation):
        return cls(lat, lon, elevation)

    def at(self, epochs):
        from spacerocks.spacerock import SpaceRock
        epochs = np.atleast_1d(epochs)

        states = []
        for epoch in epochs:
            state = get_earth_state(epoch)
            dx, dy, dz, dvx, dvy, dvz = compute_topocentric_correction(self.lon, self.lat, self.elevation, epoch)

            st = copy.deepcopy(state)

            st[0] += (dx * u.m).to(u.km).value
            st[1] += (dy * u.m).to(u.km).value
            st[2] += (dz * u.m).to(u.km).value

            st[3] += (dvx * u.m/u.day).to(u.km/u.s).value
            st[4] += (dvy * u.m/u.day).to(u.km/u.s).value
            st[5] += (dvz * u.m/u.day).to(u.km/u.s).value

            states.append(st)
        states = np.array(states).T

        units = Units()
        units.distance = u.km
        units.speed = u.km / u.s
        units.timescale = 'utc'
        units.timeformat = 'jd'

        x, y, z, vx, vy, vz = states
        
        rocks = SpaceRock(x=x, y=y, z=z, vx=vx, vy=vy, vz=vz, epoch=epochs, units=units, frame='J2000', origin='ssb')
        return rocks

        
        