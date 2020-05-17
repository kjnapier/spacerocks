import unittest
from spacerocks import SpaceRock
from astropy import units as u

# Mars orbital parameters on obs_date. From JPL Horizons.
obs_date = 2458982.5

mars_X = 3.629001676703048E-01 * u.au
mars_Y = -1.370988251351513E+00 * u.au
mars_Z = -3.784524466405902E-02 * u.au
mars_VX = 1.404127833852117E-02 * u.au / u.day.to(u.year)
mars_VY = 4.807854901674583E-03 * u.au / u.day.to(u.year)
mars_VZ = -2.436304874926014E-04 * u.au / u.day.to(u.year)

mars_a = 1.501181855304217E+00 * u.au
mars_e = 8.937143360621642E-02 * u.dimensionless_unscaled
mars_inc = 1.855182144321076E+00 * u.degree.to(u.rad)
mars_node = 4.935245879721101E+01 * u.degree.to(u.rad)
mars_omega = 2.917160972199084E+02 * u.degree.to(u.rad)
mars_epoch = 2.459072070701163E+06 * u.day

class TestSpacerocks(unittest.TestCase):

    def test_xyz_to_kep(self):
        mars = SpaceRock(X=mars_X, Y=mars_Y, Z=mars_Z,
                         VX=mars_VX, VY=mars_VY, VZ=mars_VZ,
                         coordinates='cartesian')
        self.assertAlmostEqual(mars.a[0], mars_a)
        self.assertAlmostEqual(mars.e[0], mars_e)
        self.assertAlmostEqual(mars.inc[0], mars_inc)
        self.assertAlmostEqual(mars.node[0], mars_node)
        self.assertAlmostEqual(mars.omega[0], mars_omega)
        self.assertAlmostEqual(mars.epoch[0], mars_epoch)


    #def test_kep_to_xyz(self):
#        self.assertAlmostEqual()

    #def test_xyz_to_equa(self):
    #    self.assertAlmostEqual()

    #def test_radec_to_hpix(self):
    #    self.assertAlmostEqual()

    #def test_helio_to_bary(self):
    #    self.assertAlmostEqual()

    #def test_bary_to_helio(self):
    #    self.assertAlmostEqual()

if __name__ == '__main__':
    unittest.main()
