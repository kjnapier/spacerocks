__version__ = '1.1.0'
__author__ = 'Kevin Napier'

# Find suffix
import sysconfig
suffix = sysconfig.get_config_var('EXT_SUFFIX')
if suffix is None:
    suffix = ".so"

# Import shared library
import os
import warnings
pymodulepath = os.path.dirname(__file__)
from ctypes import cdll
__libpath__ = pymodulepath+"/../libspacerocks"+suffix
clibspacerocks = cdll.LoadLibrary(__libpath__)

__libpath__ = pymodulepath+"/../_pyOrbfit"+suffix
cliborbfit = cdll.LoadLibrary(__libpath__)

os.environ['ORBIT_EPHEMERIS'] = pymodulepath + '/data/binEphem.423'
os.environ['ORBIT_OBSERVATORIES'] = pymodulepath + '/data/observatories.dat'


from .spacerock import SpaceRock
from .units import Units
from .orbfit import Orbfit
from .pyOrbfit import *
