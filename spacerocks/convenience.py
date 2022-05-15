import numpy as np
from astropy.time import Time
import copy

from .vector import Vector
from .units import Units

import dateutil


def infer_time_format(d, units):
    if isinstance(d, Time):
        epoch = d
    elif isinstance(d[0], str):
        dates = [dateutil.parser.parse(x, fuzzy_with_tokens=True)[0] for x in d]
        epoch = Time(dates, format='datetime', scale=units.timescale)  
    elif np.isscalar(d[0]):
        if np.all(d > 100000):
            epoch = Time(d, format='jd', scale=units.timescale)
        else:
            epoch = Time(d, format='mjd', scale=units.timescale)
            
    return epoch

def time_handler(d, units=Units()):
    d = np.atleast_1d(d)
    if units.timeformat is not None:
        epoch = Time(d, format=units.timeformat, scale=units.timescale)
    else:
        epoch = infer_time_format(d, units)
    
    return epoch


class Convenience:

    def __len__(self):
        '''
        This method allows you to use the len() function on a SpaceRocks object.
        '''
        return len(self.name)

    def __getitem__(self, idx):
        '''
        This method allows you to index a SpaceRocks object.
        '''
        p = copy.copy(self)
        for attr in self.__dict__.keys():
            if (attr != 'mu') and (attr != 'frame') and (attr != 'origin') and (attr != 'units'):# and (attr != '_x') and (attr != '_y') and (attr != '_z') and (attr != '_vx') and (attr != '_vy') and (attr != '_vz'):
                if isinstance(getattr(self, attr), Vector):
                    setattr(p, attr, getattr(self, attr)[idx])
                else:
                    setattr(p, attr, getattr(self, attr)[idx])

        return p

    # def detect_timescale(self, timevalue, timescale):
    #     if isinstance(timevalue[0], str):
    #         dates = [dateutil.parser.parse(d, fuzzy_with_tokens=True)[
    #             0] for d in timevalue]
    #         return Time(dates, format='datetime', scale=timescale)
    #     elif np.all(timevalue > 100000):
    #         return Time(timevalue, format='jd', scale=timescale)
    #     else:
    #         return Time(timevalue, format='mjd', scale=timescale)

    def detect_coords(self, kwargs):

        kep = ['a', 'e', 'q', 'inc', 'node', 'arg', 'M', 'f', 'E',
               'varpi', 't_peri', 'mean_longitude', 'true_longitude', 'v_inf', 'b', 'Q']
        xyz = ['x', 'y', 'z', 'vx', 'vy', 'vz']

        input_coords = list(kwargs.keys())

        if np.in1d(input_coords, kep).sum() == 6:
            coords = 'kep'

        elif np.in1d(input_coords, xyz).sum() == 6:
            coords = 'xyz'

        else:
            raise ValueError(
                'Invalid input coordinates. Please see the documentation for accepted input.')

        return coords