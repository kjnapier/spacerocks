from ..units import Units
from ..constants import epsilon
import numpy as np

class Detection:

    def __init__(self, units=Units(), **kwargs):

        self.ra = kwargs.get('ra', None)
        self.dec = kwargs.get('dec', None)
        self.ra_rate = kwargs.get('ra_rate', None)
        self.dec_rate = kwargs.get('dec_rate', None)
        self.epoch = kwargs.get('epoch', None)

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




