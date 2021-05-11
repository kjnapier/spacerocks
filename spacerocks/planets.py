'''
planets
'''

from spacerocks.planets import earth, jupiter

earth.mass.msun
jupiter.mass.msun
saturn.mass.msun

class Planet:

    def __init__(self):
        pass

    def at(self, jd):
        '''compute position at time'''
        pass


earth = Planet(name='earth', mass=..., radius=..., oblateness=..., axial_tilt=..., epoch=...)
