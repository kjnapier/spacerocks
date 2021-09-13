import unittest
from spacerocks import SpaceRock, Units

class TestSpacerocks(unittest.TestCase):

    def test_kep_to_xyz(self):

        units = Units()
        units.timescale = 'tdb'
        BP = SpaceRock(a=448.7673062316562,
                       e=0.9214543710796702,
                       inc=54.11068217470999,
                       arg=348.0587931444684,
                       node=135.2131434907198,
                       t_peri=2473015.776611103,
                       epoch=2458982.5,
                       name='2015 BP519',
                       origin='ssb',
                       units=units)

        self.assertAlmostEqual(BP.x.au[0], 15.80639409220872)
        self.assertAlmostEqual(BP.y.au[0], 26.34085915326679)
        self.assertAlmostEqual(BP.z.au[0], -41.22486401689469)
        self.assertAlmostEqual(BP.vx.value[0], -0.002654346366438451)
        self.assertAlmostEqual(BP.vy.value[0], 0.0008305911427892275)
        self.assertAlmostEqual(BP.vz.value[0], 0.001769516671619466)

    def test_circular_nonplanar(self):
        # rock = SpaceRock(a=40,
        #                  e=0,
        #                  inc=5,
        #                  arg=0,
        #                  node=0,
        #                  M=90,
        #                  origin='ssb')

        pass
        


    def test_circular_planar(self):
        pass

    def test_noncircular_planar(self):
        pass

    def test_hyperbolic(self):
        pass


if __name__ == '__main__':
    unittest.main()
