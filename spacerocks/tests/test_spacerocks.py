import unittest
from spacerocks import SpaceRock, Units

units = Units()
units.timescale = 'tdb'
BP = SpaceRock(a=4.487673062316562E+02,
               e=9.214543710796702E-01,
               inc=5.411068217470999E+01,
               arg=3.480587931444684E+02,
               node=1.352131434907198E+02,
               t_peri=2.473015776611103E+06,
               epoch=2458982.5,
               H0=4.4,
               name='2015 BP519',
               frame='barycentric',
               units=units)

class TestSpacerocks(unittest.TestCase):

    def test_ephemerides(self):
        self.assertAlmostEqual(BP.x.au[0], 1.580639409220872E+01)

if __name__ == '__main__':
    unittest.main()
