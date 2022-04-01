from astropy.coordinates import Angle
from astropy.time import Time
from astropy import units as u

import numpy as np

import dateutil

from .convenience import Convenience


class Prediction(Convenience):

    def __init__(self, **kwargs):

        units = kwargs.get('units')

        if units.timeformat is None:
            if isinstance(kwargs.get('epoch')[0], str):
                self.epoch = [Time(dateutil.parser.parse(date, fuzzy_with_tokens=True)[0], format='datetime', scale=units.timescale) for date in kwargs.get('epoch')]
            elif np.all(kwargs.get('epoch') > 100000):
                self.epoch = Time(kwargs.get('epoch'),
                                  format='jd', scale=units.timescale)
            else:
                self.epoch = Time(kwargs.get('epoch'),
                                  format='mjd', scale=units.timescale)
        else:
            self.epoch = Time(kwargs.get('epoch'),
                              format=units.timeformat, scale=units.timescale)

        self.ra = Angle(kwargs.get('ra'), u.rad)
        self.dec = Angle(kwargs.get('dec'), u.rad)

        self.err_a = Angle(kwargs.get('err_a'), u.arcsec)
        self.err_b = Angle(kwargs.get('err_b'), u.arcsec)
        self.err_pa = Angle(kwargs.get('err_pa'), u.deg)

        self.elong = Angle(kwargs.get('elong'), u.deg)
        self.opp = Angle(kwargs.get('opp'), u.deg)
