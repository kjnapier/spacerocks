import unittest
from spacerocks import SpaceRock, Units

import numpy as np

from astroquery.jplhorizons import Horizons

from astropy.time import Time
from astropy import units as u
from astropy.coordinates import SkyCoord

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

        assert sep.arcsec.max() < 1
        assert obs.mag.max() < 23
        assert obs.mag.min() > 20


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
