from ..units import Units
from ..constants import epsilon
from ..vector import Vector
from ..observer import Observer

import numpy as np
from astropy import units as u
from astropy.time import Time

from functools import cached_property

class Detection:

    def __init__(self, units=Units(), **kwargs):

        self.ra = kwargs.get('ra', None) * (units.ra.to(u.rad))
        self.dec = kwargs.get('dec', None) * (units.dec.to(u.rad))
        self.ra_rate = kwargs.get('ra_rate', None) * (units.ra_rate.to(u.rad / u.day))
        self.dec_rate = kwargs.get('dec_rate', None) * (units.dec_rate.to(u.rad / u.day))
        self.epoch = Time(kwargs.get('epoch', None), format=units.timeformat, scale=units.timescale).utc.jd

        self.observatory = kwargs.get('observatory', 500)
        self.objid = kwargs.get('objid', None)
        self.mag = kwargs.get('mag', None)
        self.mag_err = kwargs.get('mag_err', None)
        self.filter = kwargs.get('filter', None)
        self.cov = kwargs.get('cov', None)

    @property
    def b(self):
        return np.arcsin(np.cos(epsilon) * np.sin(self.dec) - np.sin(epsilon) * np.cos(self.dec) * np.sin(self.ra))

    @property
    def l(self):
        return np.arctan2((np.cos(epsilon) * np.cos(self.dec) * np.sin(self.ra) + np.sin(epsilon) * np.sin(self.dec)), (np.cos(self.dec) * np.cos(self.ra)))
    
    @property
    def b_rate(self):
        num = self.dec_rate * (np.cos(epsilon) * np.cos(self.dec) + np.sin(epsilon) * np.sin(
            self.dec) * np.sin(self.ra)) - self.ra_rate * np.cos(self.dec) * np.cos(self.ra) * np.sin(epsilon)
        denom = np.sqrt(1 - (np.cos(epsilon) * np.sin(self.dec) -
                     np.cos(self.dec) * np.sin(self.ra) * np.sin(epsilon))**2)
        return num / denom

    @property
    def l_rate(self):
        num = self.ra_rate * (1 / np.cos(self.ra)) * (np.cos(epsilon)
                                                   * (1 / np.cos(self.ra)) + np.sin(epsilon) * np.tan(self.ra))
        denom = 1 + ((1 / np.cos(self.ra)) * np.sin(epsilon) +
                     np.cos(epsilon) * np.tan(self.ra))**2
        return num / denom
    
    @cached_property
    def pointing_vector(self):
        return Vector(np.cos(self.dec) * np.cos(self.ra), 
                      np.cos(self.dec) * np.sin(self.ra), 
                      np.sin(self.dec))
    
    @cached_property
    def observer(self):
        o = Observer(obscode=self.observatory, epoch=self.epoch, frame='J2000')
        o.position = Vector(o.x.au[0], o.y.au[0], o.z.au[0])
        # o.velocity = Vector(o.vx.value, o.vy.value, o.vz.value)
        # o.rdot = o.position.unit.dot(o.velocity)
        # o.r = o.position.norm
        return o




