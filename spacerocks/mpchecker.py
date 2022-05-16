from spacerocks import SpaceRock, Units

import numpy as np

import warnings

class MPChecker:
    def __init__(self, update = False):
        self.rocks = self.load_rocks(update = update)
    
    def load_rocks(self, update):
        from .paths import MPC_PATH
        import pathlib
        rocksfile = pathlib.Path(MPC_PATH+'/mpcorb_extended.json.gz')
        #rocksfile = pathlib.Path(MPC_PATH+'mpcorb_extended.rocks')
        if rocksfile.is_file() and update == False:
            rocks = SpaceRock.from_mpc('mpcorb_extended', download_data=False, metadata='Orbit_type')
            #rocks = SpaceRock.from_file(rocksfile)
        else:
            rocks = SpaceRock.from_mpc('mpcorb_extended', download_data=True, metadata='Orbit_type')
            #rocks.write_to(MPC_PATH+'mpcorb_extended.rocks', compression=None)
        return rocks
    
    def __Great_Circle_Distances(self, obsra, obsdec, ra, dec):
        return np.arccos(np.sin(obsdec.radian)*np.sin(dec.radian) + 
                         np.cos(obsdec.radian)*np.cos(dec.radian)*np.cos(obsra.radian-ra.radian))
    
    def checker(self, ra, dec, radius, epoch, lim_mag, code):
        prop = self.rocks.analytic_propagate(epoch=epoch)
        obs = prop.observe(obscode=code)
        arc_dis = self.__Great_Circle_Distances(obs.ra, obs.dec, ra, dec)
        if radius.deg > 3:
            warnings.warn('Exceeed maximum search radius, use radius = 3 degrees instead.')
        in_field_radius = max([3*np.pi/180, 5*radius.radian, 1*np.pi/180])
        in_field = arc_dis < in_field_radius
        
        if in_field.sum() >0:
            rocks = self.rocks[in_field]
            units = Units()
            units.timescale = 'utc'
            prop, _, _ = rocks.propagate(epochs=epoch, model='PLANETS', units=units)
            obs = prop.observe(obscode=code)
            arc_dis = self.__Great_Circle_Distances(obs.ra, obs.dec, ra, dec)
            in_field = arc_dis < radius.radian
            bright = obs.mag < lim_mag
            gotcha = in_field * bright
            if gotcha.sum() > 0: 
                return rocks[gotcha], obs[gotcha]
            else:
                warnings.warn('No known object! Return None.')
                return None, None        
        else:
            warnings.warn('No known object! Return None.')
            return None, None