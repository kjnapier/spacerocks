from astropy.constants import c
from astropy import units as u
from astropy.coordinates import Angle

Mercury_mass = 1.6601367952719304E-07
Venus_mass = 2.4478383396645447E-06
Earth_mass = 3.0404326462685257E-06
Mars_mass = 3.2271514450538743E-07
Sun_mass = 1 + Mercury_mass + Venus_mass + Earth_mass + Mars_mass
Jupiter_mass = 9.547919384243222E-04
Saturn_mass = 2.858859806661029E-04
Uranus_mass = 4.3662440433515637E-05
Neptune_mass = 5.151389020535497E-05
Pluto_mass = 7.361781606089469e-09

# standard gravitational parameter, GM. This value comes directly from Horizons.
mu_bary = 0.00029630927492415936 * u.radian**2 * u.au**3 / u.day**2
mu_helio = 0.00029591220828559093 * u.radian**2 * u.au**3 / u.day**2

# obliquity of the Earth. From Horizons.
epsilon = Angle(84381.448, u.arcsec).radian

# speed of light
c = c.to(u.au / u.year)

# G parameter for magnitude estimation
G = 0.15
