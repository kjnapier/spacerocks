import warnings
from healpy.pixelfunc import ang2pix

import numpy as np
from astropy import units as u
import pandas as pd

from .constants import *

from astropy.time import Time
from astropy.coordinates import SkyCoord
from skyfield.api import Topos, Loader
# Load in planets for ephemeride calculation.
load = Loader('./Skyfield-Data', expire=False, verbose=False)
ts = load.timescale()
planets = load('de423.bsp')
sun = planets['sun']
earth = planets['earth']

class Transformations:

    #def equa_to_ecl(self, x, y, z):

    #    R = np.array([[1, 0, 0],
    #                  [0, np.cos(epsilon), np.sin(epsilon)],
    #                  [0, -np.sin(epsilon), np.cos(epsilon)]])

    #    return R @ np.array([x, y, z])

    #def ecl_to_tel(self, x, y, z):

    #    #l =
    #    #b =

    #    R = np.array([[-np.sin(l), -np.cos(l) * np.sin(b), np.cos(l) * np.cos(b)],
    #                  [np.cos(l), -np.sin(l) * np.sin(b), np.sin(l) * np.cos(b)],
    #                  [0, np.cos(b), np.sin(b)]])

    #    return R @ np.array([x, y, z])

    #def xyz_to_abg(self):

    #    earth = planets['earth']

    #    if self.__class__.obscode != 500:
    #        earth += Topos(latitude_degrees=self.__class__.obslat,
    #                        longitude_degrees=self.__class__.obslon,
    #                        elevation_m=self.__class__.obselev) # topocentric calculation

    #    t = ts.tt(jd=self.epoch.tt.jd)
    #    x_earth, y_earth, z_earth = earth.at(t).position.au * u.au # earth ICRS position
    #    vx_earth, vy_earth, vz_earth = earth.at(t).velocity.au_per_d * u.au / u.day # earth ICRS position

    #    x, y, z = self.equa_to_ecl(self.x.value, self.y.value, self.z.value) * u.au
    #    vx, vy, vz = self.equa_to_ecl(self.vx.value, self.vy.value, self.vz.value) * u.au / u.day

    #    #x0 = self.x - x_earth
    #    #y0 = self.y - y_earth
    #    #z0 = self.z - z_earth
    #    #vx0 = self.vx - vx_earth
    #    #vy0 = self.vy - vy_earth
    #    #vz0 = self.vz - vz_earth

    #    x0 = x - x_earth
    #    y0 = y - y_earth
    #    z0 = z - z_earth
    #    vx0 = vx - vx_earth
    #    vy0 = vy - vy_earth
    #    vz0 = vz - vz_earth

    #    self.alpha = x0 / z0
    #    self.beta = y0 / z0
    #    self.gamma = 1 / z0
    #    self.dalpha = vx0 / z0
    #    self.dbeta = vy0 / z0
    #    self.dgamma = vz0 / z0

    #    return self

    #def abg_to_xyz(self):

    #    earth = planets['earth']

    #    if self.__class__.obscode != 500:
    #        earth += Topos(latitude_degrees=self.__class__.obslat,
    #                        longitude_degrees=self.__class__.obslon,
    #                        elevation_m=self.__class__.obselev) # topocentric calculation

    #    t = ts.tt(jd=self.epoch.tt.jd)
    #    x_earth, y_earth, z_earth = earth.at(t).position.au * u.au # earth ICRS position
    #    vx_earth, vy_earth, vz_earth = earth.at(t).velocity.au_per_d * u.au / u.day # earth ICRS position

    #    self.x = self.alpha / self.gamma + x_earth
    #    self.y = self.beta / self.gamma + y_earth
    #    self.z = 1 / self.gamma + z_earth
    #    self.vx = self.alpha_dot / self.gamma
    #    self.vy = self.beta_dot / self.gamma
    #    self.vz = self.gamma_dot / self.gamma

    #    return self

    def radec_to_hpix(self, NSIDE):
        '''
        Convert (ra, dec) into healpix.
        '''
        return hp.pixelfunc.ang2pix(NSIDE, np.pi/2 - self.dec.rad, self.ra.rad, nest=True)


    def kep_to_xyz_temp(self, a, e, inc, arg, node, M):
        '''
        Just compute the xyz position of an object. Used for iterative equatorial
        calculation.
        '''
        # compute eccentric anomaly E
        E = self.calc_E(e.value, M.value) * u.rad

        # compute true anomaly ν
        ν = 2 * np.arctan2((1 + e)**0.5*np.sin(E/2), (1 - e)**0.5*np.cos(E/2))

        # compute the distance to the central body r
        r = a * (1 - e * np.cos(E))

        # obtain the position vector o
        o = [r * np.cos(ν), r * np.sin(ν), np.zeros(len(ν))]
        ov = [(mu_bary * a)**0.5 / r * (-np.sin(E))/ u.rad,
              (mu_bary * a)**0.5 / r * ((1 - e**2)**0.5 * np.cos(E))/ u.rad,
              np.zeros(len(ν))/ u.rad]

        # Rotate o to the inertial frame
        x, y, z = euler_rotation(arg, inc, node, o) #* u.au
        vx, vy, vz = euler_rotation(arg, inc, node, ov) #* u.au / u.day

        return x, y, z, vx, vy, vz
