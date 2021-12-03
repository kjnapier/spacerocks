from ctypes import cdll
import os
__version__ = '2.0.0'
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

from .spacerock import SpaceRock
from .units import Units