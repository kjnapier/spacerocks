from astropy.constants import c
from astropy import units as u
from astropy.coordinates import Angle


# standard gravitational parameter, GM. This value comes directly from Horizons.
#mu_bary = 0.00029630926160203097 * u.radian**2 * u.au**3 / u.day**2 # Planets only
#mu_bary = 0.00029630927492415936 * u.radian**2 * u.au**3 / u.day**2
mu_bary = 0.00029630927493457475 * u.radian**2 * u.au**3 / u.day**2
mu_helio = 0.00029591220828411951 * u.radian**2 * u.au**3 / u.day**2

# obliquity of the Earth. From Horizons.
epsilon = Angle(84381.448, u.arcsec).rad

# speed of light
c = c.to(u.au / u.day)

# G parameter for magnitude estimation
G = 0.15
