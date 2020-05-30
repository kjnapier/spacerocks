import unittest
from ..spacerocks import SpaceRock

BP = SpaceRock(a=4.487673062316562E+02,
               e=9.214543710796702E-01,
               inc=5.411068217470999E+01,
               arg=3.480587931444684E+02,
               node=1.352131434907198E+02,
               t_peri=2.473015776611103E+06,
               obsdate=2458982.5,
               H=4.5,
               name='BP',
               input_coordinates='keplerian',
               input_frame='barycentric',
               input_angles='degrees',
               obscode=304)

class TestSpacerocks(unittest.TestCase):

    def test_BP_ephemerides(self):
        self.assertAlmostEqual(BP.ra.rad[0], 1.1921539698)
        self.assertAlmostEqual(BP.dec.rad[0], -0.5464390416)
        self.assertAlmostEqual(BP.mag.value[0], 21.6360133068)
        self.assertAlmostEqual(BP.M.value[0], -0.0254097208)

if __name__ == '__main__':
    unittest.main()
