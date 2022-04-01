from ctypes import cdll
import os
__version__ = '2.1.0'
__author__ = 'Kevin Napier'

# Find suffix
import sysconfig
suffix = sysconfig.get_config_var('EXT_SUFFIX')
if suffix is None:
    suffix = ".so"

# Import shared libraries
pymodulepath = os.path.dirname(__file__)

__libspacerockspath__ = pymodulepath + "/../libspacerocks" + suffix
clibspacerocks = cdll.LoadLibrary(__libspacerockspath__)

import os
pymodulepath = os.path.dirname(__file__)
os.environ['ORBIT_EPHEMERIS'] = pymodulepath + '/data/binEphem.423'
os.environ['ORBIT_OBSERVATORIES'] = pymodulepath + '/data/observatories.dat'

from .spacerock import SpaceRock
from .units import Units