from .units import Units
from .spacerock import SpaceRock

import spiceypy as spice


class SpiceBody:

    def __init__(self, spiceid, frame='ECLIPJ2000', origin='ssb'):
        self.spiceid = spiceid
        self.frame = frame
        self.origin = origin
        
    def at(self, epoch):

        ephemeris_times = spice.str2et(['JD{} TDB'.format(epoch)])
        state, _ = spice.spkezr(self.spiceid, ephemeris_times, self.frame, 'none', self.origin)
        state = np.vstack(state).T
        x, y, z, vx, vy, vz = state 

        units = Units()
        units.distance = u.km
        units.speed = u.km / u.s
        
        rocks = SpaceRock(x=x, 
                          y=y, 
                          z=z, 
                          vx=vx, 
                          vy=vy, 
                          vz=vz, 
                          units=units, 
                          epoch=epoch, 
                          origin=self.origin,
                          frame=self.frame)
    
class Earth(SpiceBody):

    def topos(self, latitude, longitude, elevation):
        pass

    def rotation_model(self):
        pass

# class Planet(SpiceBody):

#     # self.mass = mass
#     # self.radius = radius
#     # self.albedo = albedo
#     pass
