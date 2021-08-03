__version__ = '1.0.7'
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

from .spacerock import SpaceRock
from .units import Units
#from .propagate import Propagate
#from .observe import Observe
#from .neptune_resonances import *
#from .linalg3d import *
#from .constants import *
#from .jacobians import *
