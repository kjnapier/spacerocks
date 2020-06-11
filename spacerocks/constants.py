from astropy.constants import c
from astropy import units as u
from astropy.coordinates import Angle

#M_mercury = 1.6601367952719304e-7
#M_venus = 2.4478383396645447e-6
#M_earth = 3.040432648022642e-6 # Earth-Moon Barycenter
## Moon_mass = 3.69396868e-8
#M_mars = 3.2271560375549977e-7 # Mars Barycenter
##M_inner = 1 + Mercury_mass + Venus_mass + Earth_mass + Mars_mass
#M_sun = 1
#M_jupiter = 9.547919384243222e-4
#M_saturn = 2.858859806661029e-4
#M_uranus = 4.3662440433515637e-5
#M_neptune = 5.151389020535497e-5
#M_pluto = 7.361781606089469e-9

#J4_sun = 0
#J4_mercury = 0 # Murray & Dermott
#J4_venus = 2e-6 # Murray & Dermott
#J4_earth = -2e-6 # Murray & Dermott
#J4_moon = 0
#J4_mars = -19e-6 # Murray & Dermott
#J4_jupiter = -587e-6 # Murray & Dermott
#J4_saturn = -915e-6 # Murray & Dermott
#J4_uranus = -29e-6 # Murray & Dermott
#J4_neptune = -35e-6 # Murray & Dermott

# standard gravitational parameter, GM. This value comes directly from Horizons.
mu_bary = 0.00029630927492415936 * u.radian**2 * u.au**3 / u.day**2
mu_helio = 0.00029591220828559093 * u.radian**2 * u.au**3 / u.day**2

# obliquity of the Earth. From Horizons.
epsilon = Angle(84381.448, u.arcsec).radian

# speed of light
c = c.to(u.au / u.year)

# G parameter for magnitude estimation
G = 0.15
