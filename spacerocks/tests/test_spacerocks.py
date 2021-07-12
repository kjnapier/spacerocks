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
                       H0=4.4,
                       name='2015 BP519',
                       frame='barycentric',
                       units=units)

        self.assertAlmostEqual(BP.x.au[0], 15.80639409220872)
        self.assertAlmostEqual(BP.y.au[0], 26.34085915326679)
        self.assertAlmostEqual(BP.z.au[0], -41.22486401689469)
        self.assertAlmostEqual(BP.vx.value[0], -0.002654346366438451)
        self.assertAlmostEqual(BP.vy.value[0], 0.0008305911427892275)
        self.assertAlmostEqual(BP.vz.value[0], 0.001769516671619466)


    def test_kep_to_xyz(self):

        units = Units()
        units.timescale = 'tdb'
        BP = SpaceRock(x=15.80639409220872,
                       y=26.34085915326679,
                       z=-41.22486401689469,
                       vx=-0.002654346366438451,
                       vy=0.0008305911427892275,
                       vz=0.001769516671619466,
                       epoch=2458982.5,
                       H0=4.4,
                       name='2015 BP519',
                       frame='barycentric',
                       units=units)

        self.assertAlmostEqual(BP.a.au[0], 448.7673062316562)
        self.assertAlmostEqual(BP.e[0], 0.9214543710796702)
        self.assertAlmostEqual(BP.inc.deg[0], 54.11068217470999)
        self.assertAlmostEqual(BP.arg.deg[0], 348.0587931444684)
        self.assertAlmostEqual(BP.node.deg[0], 135.2131434907198)
        self.assertAlmostEqual(BP.M.deg[0], 358.5441302379153)
        self.assertAlmostEqual(BP.true_anomaly.deg[0], 290.1471338857519)
        self.assertAlmostEqual(BP.t_peri.jd[0], 2473015.776611103)
        self.assertAlmostEqual(BP.q.au[0], 35.24871030684755)


if __name__ == '__main__':
    unittest.main()
