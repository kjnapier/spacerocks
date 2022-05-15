from astropy.units import solMass
import numpy as np

SunGM       =  4. * np.pi * np.pi / 1.000037773533  # solar gravitation
MercuryGM   =  6.55371264e-06
VenusGM     =  9.66331433e-05
EarthMoonGM =  1.20026937e-04
MarsGM      =  1.27397978e-05
JupiterGM   =  3.76844407e-02 + 7.80e-6
SaturnGM    =  1.12830982e-02 + 2.79e-6
UranusGM    =  1.72348553e-03 + 0.18e-6
NeptuneGM   =  2.03318556e-03 + 0.43e-6

SolarSystemGM = SunGM + MercuryGM + VenusGM + EarthMoonGM + MarsGM + JupiterGM + SaturnGM + UranusGM + NeptuneGM


orbitspp_masses = {'Sun': (SunGM + MercuryGM + VenusGM + EarthMoonGM + MarsGM)/(4 * np.pi * np.pi) * solMass,
                   'Jupiter Barycenter' : JupiterGM/(4 * np.pi * np.pi) * solMass,
                   'Saturn Barycenter' : SaturnGM/(4 * np.pi * np.pi) * solMass,
                   'Uranus Barycenter' : UranusGM/(4 * np.pi * np.pi) * solMass,
                   'Neptune Barycenter' : NeptuneGM/(4 * np.pi * np.pi) * solMass}
