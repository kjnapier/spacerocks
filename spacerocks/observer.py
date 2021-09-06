from . import spacerock
from .units import Units

from astropy import units as u
import numpy as np
import pandas as pd
import spiceypy as spice

import pkg_resources

DATA_PATH = pkg_resources.resource_filename('spacerocks', 'data/observatories.csv')
observatories = pd.read_csv(DATA_PATH)

SPICE_PATH = pkg_resources.resource_filename('spacerocks', 'data/spice-test.txt')
spice.furnsh(SPICE_PATH)

class Observer(SpaceRock, Earth, SpiceBody):

    def __init__(self, origin='ssb', frame='ECLIPJ2000', **kwargs):
        
        self.origin = origin
        self.frame = frame

        if kwargs.get('obscode') is not None:
            self.obscode = kwargs.get('obscode')

        elif kwargs.get('spiceid') is not None:
            self.spiceid = kwargs.get('spiceid')

        else:
            raise ValueError('Must specify either a spiceid or an obscode')

        if kwargs.get('epoch') is not None:
            self.epoch = kwargs.get('epoch')

        ephemeris_time = spice.str2et(['JD{} UTC'.format(ep) for ep in self.epoch.utc.jd])

        state, _ = spice.spkezr(self.obscode, ephemeris_time, self.frame, 'none', self.origin)
        x, y, z, vx, vy, vz = np.vstack(state).T

        


    def at(self, epoch):

        ephemeris_time = spice.str2et(['JD{} UTC'.format(ep) for ep in epoch.utc.jd])

        if hasattr(self, 'obscode'):

            state, _ = spice.spkezr(self.obscode, ephemeris_time, self.frame, 'none', self.origin)
            x, y, z, vx, vy, vz = np.vstack(state).T

            pass

        elif hasattr(self, 'spiceid'):
    
            state, _ = spice.spkezr(self.spiceid, ephemeris_time, self.frame, 'none', self.origin)
            x, y, z, vx, vy, vz = np.vstack(state).T
        
            units = Units()
            units.distance = u.km
            units.speed = u.km / u.s
            units.timescale = 'tdb'
        
            rocks = spacerock.SpaceRock(x=x,
                                        y=y,
                                        z=z,
                                        vx=vx,
                                        vy=vy,
                                        vz=vz,
                                        units=units,
                                        epoch=epoch.tdb.jd,
                                        origin=self.origin,
                                        frame=self.frame)
            
            return rocks

        else:
            raise ValueError('Must provide a known spiceid or obscode.')

    def __decode_obscode(self):
        pass