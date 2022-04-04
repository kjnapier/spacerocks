from ctypes import cdll
import pathlib
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

import pkg_resources
SPICE_PATH = pkg_resources.resource_filename('spacerocks', 'data/spice')

p = pathlib.Path(SPICE_PATH)
p.mkdir(parents=True, exist_ok=True)

required_spice = {'de440s.bsp': 'https://ssd.jpl.nasa.gov/ftp/eph/planets/bsp/de440s.bsp', 
                  'sb441-n16s.bsp': 'https://ssd.jpl.nasa.gov/ftp/xfr/sb441-n16s.bsp',
                  'gm_Horizons.pck': 'https://ssd.jpl.nasa.gov/ftp/xfr/gm_Horizons.pck',
                  'latest_leapseconds.tls': 'https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/latest_leapseconds.tls',
                  'pck00010.tpc': 'https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00010.tpc'}

to_download = []
for filename, url in required_spice.items():
    if not p.joinpath(filename).exists():
        to_download.append(url)

if len(to_download) > 0:
    from .downloader import download
    download(to_download, SPICE_PATH)

from .spacerock import SpaceRock
from .units import Units