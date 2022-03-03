import unittest
from spacerocks import SpaceRock, Units

import numpy as np

from astroquery.jplhorizons import Horizons

from astropy.time import Time
from astropy import units as u
from astropy.coordinates import SkyCoord

class TestSpacerock(unittest.TestCase):

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

        a = 448.7673062316562
        e = 0.9214543710796702
        self.assertAlmostEqual(BP.q.au[0], a * (1 - e))
        self.assertAlmostEqual(BP.b.au[0], a * (1 - e*e)**0.5)


        x, y, z = BP.orbits(N=100)
        self.assertTrue(len(x[0]) == len(y[0]) == len(z[0]) == 100)
        BP.change_frame('J2000')
        BP.to_helio()

        BP.a 
        BP.e 
        BP.inc 
        BP.arg 
        BP.node 
        BP.varpi 
        BP.M 
        BP.E 
        BP.f 
        BP.true_longitude 
        BP.mean_longitude 
        BP.q 
        BP.t_peri 
        BP.b 
        BP.p 
        BP.n 
        BP.Q 
        BP.r 
        BP.ovec 
        BP.vovec 
        BP.position 
        BP.velocity 
        BP.rrdot 
        BP.x 
        BP.y 
        BP.z 
        BP.vx 
        BP.vy 
        BP.vz
        BP.v_inf

        del BP.a 
        del BP.e 
        del BP.inc 
        del BP.arg 
        del BP.node 
        del BP.varpi 
        del BP.M 
        del BP.E 
        del BP.f 
        del BP.true_longitude 
        del BP.mean_longitude 
        del BP.q 
        del BP.t_peri 
        del BP.b 
        del BP.p 
        del BP.n 
        del BP.Q 
        del BP.r 
        del BP.ovec 
        del BP.vovec 
        del BP.position 
        del BP.velocity 
        del BP.v_inf
        
        BP.position
        BP.velocity

        del BP.x 
        del BP.y 
        del BP.z 
        del BP.vx 
        del BP.vy 
        del BP.vz

    def test_from_horizons(self):

        startdate = Time('2000-01-01', scale='utc', format='iso')
        enddate   = Time('2050-01-01', scale='utc', format='iso')
        testdates = Time(np.arange(startdate.jd, enddate.jd, 30), scale='utc', format='jd')
        
        units = Units()
        units.timescale = 'utc'

        rock = SpaceRock.from_horizons(name='Ceres')
        prop, _, _ = rock.propagate(epochs=testdates.jd, model=2, units=units)
        obs = prop.observe(obscode='W84')
    
        ephem_Horizons = Horizons(id='Ceres', location='W84',
                                epochs={'start':testdates[0].iso, 
                                        'stop':testdates[-1].iso, 
                                        'step':'30d'}).ephemerides()
    
        pos_Horizons = SkyCoord(ephem_Horizons['RA'], ephem_Horizons['DEC'], frame='icrs', unit=(u.deg, u.deg))
        pos_pred = SkyCoord(obs.ra.deg, obs.dec.deg, frame='icrs', unit=(u.deg, u.deg))
        sep = pos_pred.separation(pos_Horizons)

        self.assertTrue(sep.arcsec.max() < 1)
        self.assertTrue(obs.mag.max() < 10.5)
        self.assertTrue(obs.mag.min() > 6)


    def test_multiple_constructors(self):

        def gen_rock():
            rock = SpaceRock(x=15.80639409220872, 
                             y=26.34085915326679, 
                             z=-41.22486401689469, 
                             vx=-0.002654346366438451, 
                             vy=0.0008305911427892275, 
                             vz=0.001769516671619466, epoch='2022 March 2')
            return rock

        rock = gen_rock()
        prop, planets, sim = rock.propagate(epochs='3 March 2022', model=0)
        self.assertEqual(len(planets), 1)

        rock = gen_rock()
        prop, planets, sim = rock.propagate(epochs='3 March 2022', model=1)
        self.assertEqual(len(planets), 5)

        rock = SpaceRock(e=0, q=40, inc=5, node=14, M=10, arg=0, epoch='2 March 2022', frame='J2000')
        prop, planets, sim = rock.propagate(epochs='3 March 2022', model=3)
        self.assertEqual(len(planets), 27)

        rock = gen_rock()
        rock.position
        rock.velocity
        rock.inc.deg

        rock = gen_rock()
        rock.e

        rock = gen_rock()
        rock.a

        rock = gen_rock()
        rock.node

        rock = gen_rock()
        rock.varpi

        rock = gen_rock()
        rock.M

        rock = gen_rock()
        rock.r

        rock = gen_rock()
        rock.arg

        rock = gen_rock()
        rock.true_longitude

        rock = gen_rock()
        rock.f

        rock = SpaceRock(e=0, q=40, inc=5, node=14, M=10, arg=0)
        self.assertAlmostEqual(rock.a.au, 40)
        self.assertAlmostEqual(rock.Q.au, 40)
        self.assertAlmostEqual(rock.varpi.deg, 14)

        rock = SpaceRock(Q=80, q=40, inc=5, node=14, M=10, varpi=9)
        self.assertAlmostEqual(rock.a.au[0], 60)
        self.assertAlmostEqual(rock.e, 1/3)
        self.assertAlmostEqual(rock.arg.deg[0], 355)

        rock = SpaceRock(Q=80, q=40, inc=5, arg=14, M=10, varpi=9)
        self.assertAlmostEqual(rock.node.deg[0], 355)

        rock = SpaceRock(a=80, q=40, inc=5, node=14, M=10, varpi=9)
        self.assertAlmostEqual(rock.e, 0.5)
        
        rock = SpaceRock(a=-100, e=1.5, inc=35, node=265, arg=67, M=356)
        rock.v_inf

        rock = SpaceRock(b=30, e=0.8, inc=35, node=265, arg=67, M=356)
        rock.a

        rock = SpaceRock(b=30, e=0.8, inc=35, node=265, arg=67, mean_longitude=356)
        rock.M

        rock = SpaceRock(b=30, e=0.8, inc=35, node=265, arg=67, true_longitude=356)
        rock.f

        rock = SpaceRock(b=30, e=0.8, inc=35, node=265, arg=67, E=356)
        rock = SpaceRock(b=30, e=0.8, inc=35, node=265, arg=67, f=356)

        rock = SpaceRock(a=80, e=0.1, inc=5, node=14, M=10, varpi=9, mass=1, radius=0.001)
        rock.density

        rock = SpaceRock(a=80, e=0.1, inc=5, node=14, M=10, varpi=9, mass=1, density=0.001)
        rock.radius

        rock = SpaceRock(a=80, e=0.1, inc=5, node=14, M=10, varpi=9, density=1, radius=0.001)
        rock.mass

        rock = SpaceRock(a=80, e=0.1, inc=5, node=14, M=10, varpi=9, mag=23, epoch='2 March 2022')
        rock.calc_H(obscode='W84')
        obs = rock.observe(spiceid='Earth')

        rock = SpaceRock(a=80, e=0.1, inc=5, node=14, M=10, varpi=9, H=4, epoch='2 March 2022')
        obs = rock.observe(spiceid='Earth')
        
        units = Units()
        units.timeformat = 'mjd'
        units.speed = u.km/u.s
        rock = SpaceRock(b=10, v_inf=-888, inc=135, node=87, arg=35, f=30, epoch=59581)
        rock.a.au
        rock.e
        
        with self.assertRaises(ValueError): 
            SpaceRock(a=40, e=0, inc=5, node=14, M=10, arg=100)
        with self.assertRaises(ValueError): 
            SpaceRock(Q=-1, a=40, inc=5, node=14, M=10, arg=9)
        with self.assertRaises(ValueError): 
            SpaceRock(Q=30, q=60, inc=5, node=14, M=10, arg=0)
        with self.assertRaises(ValueError): 
            SpaceRock(Q=30, q=-2, inc=5, node=14, M=10, arg=0)
        with self.assertRaises(ValueError): 
            SpaceRock(a=-40, e=0.5, inc=5, node=14, M=10, arg=100)
        with self.assertRaises(ValueError): 
            SpaceRock(a=40, e=-0.5, inc=5, node=14, M=10, arg=100)
        with self.assertRaises(ValueError): 
            SpaceRock(a=40, e=1.5, inc=5, node=14, M=10, arg=100)

        


if __name__ == '__main__':
    unittest.main()
