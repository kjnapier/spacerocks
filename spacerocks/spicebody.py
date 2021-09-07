from .units import Units
from .spacerock import SpaceRock

import spiceypy as spice
from astropy import units as u

class SpiceBody:

    def __init__(self, frame='ECLIPJ2000', origin='ssb'):

        self.frame = frame
        self.origin = origin

        if kwargs.get('epoch') is not None:
            self.epoch = np.atleast_1d(kwargs.get('epoch'))

        if len(kwargs.get('spiceid')) < len(self.epoch):
                self.spiceid = np.repeat(kwargs.get('spiceid'), len(self.epoch))
            else:
                self.spiceid = np.atleast_1d(kwargs.get('spiceid'))

    
    def __get_all_state_vectors(self):
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

        unique_rocks_and_times = set(list(zip(self.spiceid, self.epoch)))
        unique_dict = {key:self.__state_from_spice(key) for key in unique_rocks_and_times}
        x, y, z, vx, vy, vz = np.array([unique_dict[(body, time)] for body, time in zip(self.spiceid, self.epoch)]).T
       
        self.x = Distance(x, u.km, allow_negative=True).to(u.au)
        self.y = Distance(y, u.km, allow_negative=True).to(u.au)
        self.z = Distance(z, u.km, allow_negative=True).to(u.au)

        self.vx = (vx * u.km/u.s).to(u.au / u.day)
        self.vy = (vy * u.km/u.s).to(u.au / u.day)
        self.vz = (vz * u.km/u.s).to(u.au / u.day)


    def __compute_ephemeris_time(self, epoch):
        '''
        Wrapper for spiceypy's str2et function.
        '''
        return spice.str2et('JD{} UTC'.format(epoch))

    def __state_from_spice(self, x):
        spiceid, epoch = x
        ephemeris_time = self.__compute_ephemeris_time(epoch)
        state, _ = spice.spkezr(spiceid, ephemeris_time, self.frame, 'none', self.origin)
        return state