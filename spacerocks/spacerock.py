################################################################################
# SpaceRocks, version 1.0.0
#
# Author: Kevin Napier kjnapier@umich.edu
################################################################################
import random

from astropy import units as u
from astropy.coordinates import Angle, Distance
from astropy.time import Time

from numpy import sin, cos, arctan2, sqrt, array, pi

from .constants import *
from .orbitfuncs import OrbitFuncs
from .convenience import Convenience
from .units import Units
from .vector import Vector

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
