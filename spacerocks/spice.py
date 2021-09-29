#import spacerocks.spacerock as sr
#from .spacerock import SpaceRock

from astropy import units as u
from astropy.coordinates import Distance

import numpy as np

from astropy.constants import G as GravitationalConstant

import spiceypy as spice
import os

import pkg_resources

SPICE_PATH = pkg_resources.resource_filename('spacerocks', 'data/spice')
spice.furnsh(os.path.join(SPICE_PATH, 'latest_leapseconds.tls'))
spice.furnsh(os.path.join(SPICE_PATH, 'de440s.bsp'))
spice.furnsh(os.path.join(SPICE_PATH, 'hst.bsp'))
spice.furnsh(os.path.join(SPICE_PATH, 'nh.bsp'))
#spice.furnsh(os.path.join(SPICE_PATH, 'gm_de431.tpc'))
spice.furnsh(os.path.join(SPICE_PATH, 'gm_Horizons.pck'))

class SpiceBody:

    def __init__(self, spiceid, frame='ECLIPJ2000', origin='ssb'):

        self.frame = frame
        self.origin = origin
        self.spiceid = spiceid

    
    def at(self, epoch):
        '''
        Return a SpaceRock object at the specified epoch(s).
        '''
        from spacerocks.spacerock import SpaceRock
        epoch = epoch.utc.jd
        x, y, z, vx, vy, vz = self.__get_all_state_vectors(epoch)
        return SpaceRock(x=x, y=y, z=z, vx=vx, vy=vy, vz=vz, epoch=epoch, name=self.spiceid)

    @property
    def mass(self):
        '''
        Return mass of the specified body.
        '''
        return (spice.bodvrd(self.spiceid, 'GM', 1)[1][0] * u.km**3 * u.s**(-2) / GravitationalConstant).to(u.Msun)

    @property
    def mu(self):
        '''
        Return GM of the specified body.
        '''
        return (spice.bodvrd(self.spiceid, 'GM', 1)[1][0] * u.km**3 * u.s**(-2)).to(u.au**3 / u.day**2) * u.radian**2

    
    def __get_all_state_vectors(self, epoch):
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

        '''

        epoch = np.atleast_1d(epoch)
        spiceid = np.atleast_1d(self.spiceid)

        unique_rocks_and_times = set(list(zip(spiceid, epoch)))
        unique_dict = {key:self.__state_from_spice(key) for key in unique_rocks_and_times}
        x, y, z, vx, vy, vz = np.array([unique_dict[(body, time)] for body, time in zip(spiceid, epoch)]).T
       
        x = Distance(x, u.km, allow_negative=True).to(u.au)
        y = Distance(y, u.km, allow_negative=True).to(u.au)
        z = Distance(z, u.km, allow_negative=True).to(u.au)

        vx = (vx * u.km/u.s).to(u.au / u.day)
        vy = (vy * u.km/u.s).to(u.au / u.day)
        vz = (vz * u.km/u.s).to(u.au / u.day)

        return x.au, y.au, z.au, vx.value, vy.value, vz.value


    def __compute_ephemeris_time(self, epoch):
        '''
        Wrapper for spiceypy's str2et function.
        '''
        return spice.str2et('JD{} UTC'.format(epoch))


    def __state_from_spice(self, x):
        '''
        Routine for actually getting the state vector from spice.
        '''
        spiceid, epoch = x
        ephemeris_time = self.__compute_ephemeris_time(epoch)
        state, _ = spice.spkezr(spiceid, ephemeris_time, self.frame, 'none', self.origin)
        return state