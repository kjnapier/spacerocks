import rebound
import os
import spiceypy as spice
import pkg_resources
from astropy.time import Time
import numpy as np

from .spice import SpiceBody

SPICE_PATH = pkg_resources.resource_filename('spacerocks', 'data/spice')

class Simulation(rebound.Simulation):   

    def __init__(self):
        self.units = ('day', 'AU', 'Msun')

    def add_perturbers(self, names, epoch):
        names = np.atleast_1d(names)
        epoch = Time(epoch, scale='tdb', format='jd')
        for name in names:
            body = SpiceBody(spiceid=name)
            b = body.at(epoch)
            self.add(x=b.x.au, 
                     y=b.y.au, 
                     z=b.z.au,
                     vx=b.vx.value, 
                     vy=b.vy.value, 
                     vz=b.vz.value,
                     m=body.mass.value, 
                     hash=name)

    def add_spacerocks(self, rocks):
        rocks.position
        rocks.velocity
        for rock in rocks:
            self.add(x=rock.x.au, 
                     y=rock.y.au, 
                     z=rock.z.au,
                     vx=rock.vx.value, 
                     vy=rock.vy.value, 
                     vz=rock.vz.value,
                     m=0, 
                     hash=rock.name)
