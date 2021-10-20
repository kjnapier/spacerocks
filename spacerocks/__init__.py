from ctypes import cdll
import os
__version__ = '1.1.6'
__author__ = 'Kevin Napier'

# Find suffix
import sysconfig
suffix = sysconfig.get_config_var('EXT_SUFFIX')
if suffix is None:
    suffix = ".so"

# Import shared libraries
pymodulepath = os.path.dirname(__file__)

__libpath__ = pymodulepath + "/../libspacerocks" + suffix
clibspacerocks = cdll.LoadLibrary(__libpath__)

__libpath__ = pymodulepath + "/../_pyOrbfit" + suffix
cliborbfit = cdll.LoadLibrary(__libpath__)


'''
Set the environment variables for pyOrbfit. This is clunky,
but rewriting the code would be a bear.
'''
os.environ['ORBIT_EPHEMERIS'] = pymodulepath + '/data/binEphem.423'
os.environ['ORBIT_OBSERVATORIES'] = pymodulepath + '/data/observatories.dat'

from .spacerock import SpaceRock
from .units import Units
from .orbfit import Orbfit
from .gauss import gauss
from .observer import Observer
from .orbitfit import OrbitFitter