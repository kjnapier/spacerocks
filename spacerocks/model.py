from .spice import SpiceKernel, SpiceBody
from .spacerock import SpaceRock

import pkg_resources
SPICE_PATH = pkg_resources.resource_filename('spacerocks', 'data/spice')

HORIZONS_PERTURBERS = ['Sun', 
                       'Mercury Barycenter',
                       'Venus Barycenter',
                       'Earth',
                       'Moon',
                       'Mars Barycenter',
                       'Jupiter Barycenter', 
                       'Saturn Barycenter', 
                       'Uranus Barycenter', 
                       'Neptune Barycenter', 
                       'Pluto Barycenter', 'Ceres', 'Vesta', 'Pallas', 
                       '2000003', '2000007', '2000010', '2000015', '2000016', 
                       '2000031', '2000052', '2000052', '2000065', '2000087', 
                       '2000088', '2000107', '2000511', '2000704']

builtin_models = {'SUN': [['Sun'], SpiceKernel()], 
                  'GIANTS': [['Sun', 
                             'Jupiter Barycenter', 
                             'Saturn Barycenter', 
                             'Uranus Barycenter', 
                             'Neptune Barycenter'], SpiceKernel()],
                  'PLANETS': [['Sun', 
                              'Mercury Barycenter',
                              'Venus Barycenter',
                              'Earth',
                              'Moon',
                              'Mars Barycenter',
                              'Jupiter Barycenter', 
                              'Saturn Barycenter', 
                              'Uranus Barycenter', 
                              'Neptune Barycenter', 'Pluto Barycenter'], SpiceKernel()],
                  'HORIZONS': [HORIZONS_PERTURBERS, SpiceKernel(spk=['de423.bsp', 'sb441-n16s.bsp'])]}

class PerturberModel:
    
    def __init__(self, kernel=SpiceKernel(), **kwargs):
        
        self.perturbers = {}
        self.kernel = kernel
        
        if kwargs.get('spiceids') is not None:
            kernel.furnsh()
            for spiceid in kwargs.get('spiceids'):
                self.perturbers[spiceid] = SpiceBody(spiceid=spiceid)
                
        if kwargs.get('rocks') is not None:
            for rock in kwargs.get('rocks'):
                self.perturbers[rock.name[0]] = rock

    @classmethod
    def from_builtin(cls, model: str):
        spiceids, kernel = builtin_models[model]
        return cls(spiceids=spiceids, kernel=kernel)

    def add(self, body) -> None:
        if isinstance(body, SpaceRock):
            self.perturbers[body.name[0]] = body
        elif isinstance(body, str):
            self.perturbers[body] = SpiceBody(spiceid=body)
            
