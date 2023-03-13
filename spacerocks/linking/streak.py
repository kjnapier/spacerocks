from ..spacerock import SpaceRock
from ..units import Units
from ..observer import Observer
from ..vector import Vector
from ..utils import time_handler
from ..cbindings import kepM_to_xyz
from ..constants import c as speed_of_light
from ..constants import mu_bary

from astropy.coordinates import Angle
from astropy.time import Time
from astropy import units as u

import numpy as np

import copy

class Streak:
 
    def __init__(self, ra, dec, ra_rate, dec_rate, epoch, obscode, units: Units = Units()):
 
        self.ra = Angle(np.atleast_1d(ra), units.ra).rad
        self.dec = Angle(np.atleast_1d(dec), units.dec).rad
        self.ra_rate = (np.atleast_1d(ra_rate) * units.ra_rate).to(u.rad/u.day).value
        self.dec_rate = (np.atleast_1d(dec_rate) * units.dec_rate).to(u.rad/u.day).value
        self.epoch = time_handler(epoch, units)
        self.obscode = obscode

    def __getitem__(self, idx):
        '''
        This method allows you to index a SpaceRocks object.
        '''
        p = copy.copy(self)
        for attr in self.__dict__.keys():
            if (attr != 'obscode'):
                setattr(p, attr, getattr(self, attr)[idx])
        return p
 
    @property
    def pointing_vector(self):
        return Vector(np.cos(self.dec) * np.cos(self.ra), 
                      np.cos(self.dec) * np.sin(self.ra), 
                      np.sin(self.dec))

    @property
    def pointing_vector_rate(self):
        return Vector(-np.cos(self.dec) * np.sin(self.ra) * self.ra_rate - np.sin(self.dec) * np.cos(self.ra) * self.dec_rate, 
                       np.cos(self.dec) * np.cos(self.ra) * self.ra_rate - np.sin(self.dec) * np.sin(self.ra) * self.dec_rate, 
                       np.cos(self.dec) * self.dec_rate)

    @property
    def solar_elongation(self):
        return np.arccos(-self.pointing_vector.unit.dot(self.observer.position.unit))

    @property
    def solar_elongation_rate_times_sin_solar_elongation(self):
        o = self.observer
        xx = (o.r * o.velocity.x - o.position.x * o.rdot) / o.r**2
        yy = (o.r * o.velocity.y - o.position.y * o.rdot) / o.r**2
        zz = (o.r * o.velocity.z - o.position.z * o.rdot) / o.r**2
        r_hat_dot = Vector(xx, yy, zz)
        elongation_rate = (self.pointing_vector.dot(r_hat_dot) + o.position.unit.dot(self.pointing_vector_rate)) #/ np.sin(self.solar_elongation)
        return elongation_rate

    @property
    def observer(self):
        o = Observer(obscode=self.obscode, epoch=self.epoch.utc.jd, frame='J2000')
        o.position = Vector(o.x.au, o.y.au, o.z.au)
        o.velocity = Vector(o.vx.value, o.vy.value, o.vz.value)
        o.rdot = o.position.unit.dot(o.velocity)
        o.r = o.position.norm
        return o

    def calculate_rho(self, r: float) -> float:
        A = self.observer.r * np.cos(self.solar_elongation)
        B = np.sqrt(self.observer.r**2 * (np.cos(self.solar_elongation)**2 - 1) + r**2)
        return A + B

    def calculate_rho_rate(self, r, r_rate):
        o = self.observer
        A = o.rdot * np.cos(self.solar_elongation)
        B = o.r * self.solar_elongation_rate_times_sin_solar_elongation
        C = r * r_rate - o.r * (np.sin(self.solar_elongation)**2 * o.rdot \
            + np.cos(self.solar_elongation) * o.r * self.solar_elongation_rate_times_sin_solar_elongation)
        D = np.sqrt(r**2 - o.r**2 * np.sin(self.solar_elongation)**2)
        return A - B + (C / D)

    def solve_for_range_and_rate(self, r, r_rate, t0):
        dt = self.epoch.utc.jd - t0
        r_guess = r + r_rate * dt + 0.5 * (mu_bary.value / r**3) * dt**2
        r_rate_guess = r_rate + mu_bary.value / r**3 * dt

        converged = False
        N_iter = 0
        while not converged:    
    
            rho = self.calculate_rho(r_guess)
            rho_rate = self.calculate_rho_rate(r_guess, r_rate_guess)

            ltt = rho / speed_of_light.value
            
            x = self.observer.position.x + rho * self.pointing_vector.x
            y = self.observer.position.y + rho * self.pointing_vector.y
            z = self.observer.position.z + rho * self.pointing_vector.z
            vx = self.observer.velocity.x + rho_rate * self.pointing_vector.x + rho * self.pointing_vector_rate.x
            vy = self.observer.velocity.y + rho_rate * self.pointing_vector.y + rho * self.pointing_vector_rate.y
            vz = self.observer.velocity.z + rho_rate * self.pointing_vector.z + rho * self.pointing_vector_rate.z

            units = Units()
            units.timescale = 'utc'
            units.timeformat = 'jd'
            rock = SpaceRock(x=x,
                             y=y, 
                             z=z, 
                             vx=vx, 
                             vy=vy, 
                             vz=vz, 
                             epoch=self.epoch.utc.jd - ltt, 
                             name=list(range(len(x))),
                             frame='J2000', 
                             origin='ssb', 
                             units=units)

            prop = rock.analytic_propagate(epoch=t0, units=units)

            dr = r - prop.r.au
            dr_rate = r_rate - prop.position.unit.dot(prop.velocity).value

            
            if abs(max(dr)) < 1e-5 and abs(max(dr_rate)) < 1e-8:
                converged = True
            else:
                N_iter += 1
                print(N_iter)
                if N_iter > 10:
                    return r_guess, r_rate_guess                    
                else:
                    r_guess += dr
                    r_rate_guess += dr_rate

        return r_guess, r_rate_guess

    def generate_orbits(self, r: float, r_rate: float) -> SpaceRock:

        rho = self.calculate_rho(r)
        rho_rate = self.calculate_rho_rate(r, r_rate)

        x = self.observer.position.x + rho * self.pointing_vector.x
        y = self.observer.position.y + rho * self.pointing_vector.y
        z = self.observer.position.z + rho * self.pointing_vector.z
        vx = self.observer.velocity.x + rho_rate * self.pointing_vector.x + rho * self.pointing_vector_rate.x
        vy = self.observer.velocity.y + rho_rate * self.pointing_vector.y + rho * self.pointing_vector_rate.y
        vz = self.observer.velocity.z + rho_rate * self.pointing_vector.z + rho * self.pointing_vector_rate.z

        ltt = rho / speed_of_light.value
        units = Units()
        units.timescale = 'utc'
        units.timeformat = 'jd'
        rocks = SpaceRock(x=x,
                          y=y, 
                          z=z, 
                          vx=vx, 
                          vy=vy, 
                          vz=vz, 
                          epoch=self.epoch.utc.jd - ltt, 
                          frame='J2000', 
                          origin='ssb', 
                          units=units)
        return rocks