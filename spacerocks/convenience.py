import numpy as np
from astropy.time import Time
import copy

from .vector import Vector
from .units import Units


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
            if (attr != 'mu') and (attr != 'frame') and (attr != 'origin') and (attr != 'units'):
                if isinstance(getattr(self, attr), Vector):
                    setattr(p, attr, getattr(self, attr)[idx])
                else:
                    setattr(p, attr, getattr(self, attr)[idx])

        return p

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