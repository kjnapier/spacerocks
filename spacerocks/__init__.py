from ctypes import cdll
import pathlib
import os
import shutil

__version__ = '2.3.4'
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

if os.path.exists(pymodulepath + '/pyOrbfit.py'):
    shutil.move(pymodulepath + '/pyOrbfit.py', pymodulepath + '/pyorbfit/pyOrbfit.py')

os.environ['ORBIT_EPHEMERIS'] = pymodulepath + '/data/pyOrbfit/binEphem.423'
os.environ['ORBIT_OBSERVATORIES'] = pymodulepath + '/data/pyOrbfit/observatories.dat'

from .paths import SPICE_PATH
p = pathlib.Path(SPICE_PATH)
p.mkdir(parents=True, exist_ok=True)

required_spice = {'de423.bsp': 'https://ssd.jpl.nasa.gov/ftp/eph/planets/bsp/de423.bsp',
                  'sb441-n16s.bsp': 'https://ssd.jpl.nasa.gov/ftp/xfr/sb441-n16s.bsp',
                  'gm_Horizons.pck': 'https://ssd.jpl.nasa.gov/ftp/xfr/gm_Horizons.pck',
                  'latest_leapseconds.tls': 'https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/latest_leapseconds.tls',
                  'pck00010.tpc': 'https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00010.tpc',
                  'nh_pred_alleph_od151.bsp': 'https://naif.jpl.nasa.gov/pub/naif/NEWHORIZONS/misc/pds_rel0005/nh_pred_alleph_od151.bsp'}

#'de440s.bsp': 'https://ssd.jpl.nasa.gov/ftp/eph/planets/bsp/de440s.bsp',
#'de441.bsp': 'https://ssd.jpl.nasa.gov/ftp/eph/planets/bsp/de441.bsp',

to_download = []
for filename, url in required_spice.items():
    if not p.joinpath(filename).exists():
        to_download.append(url)

if len(to_download) > 0:
    from .downloader import download
    download(to_download, SPICE_PATH)

from .spacerock import SpaceRock
from .units import Units