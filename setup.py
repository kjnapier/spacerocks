import numpy
from setuptools import setup, Extension

import os
import sys

import sysconfig
suffix = sysconfig.get_config_var('EXT_SUFFIX')
if suffix is None:
    suffix = ".so"

extra_link_args = []
if sys.platform == 'darwin':
    from distutils import sysconfig
    vars = sysconfig.get_config_vars()
    vars['LDSHARED'] = vars['LDSHARED'].replace('-bundle', '-shared')
    extra_link_args = ['-Wl,-install_name,@rpath/libspacerocks' + suffix]

libspacerocksmodule = Extension('libspacerocks',
                                sources=['src/speedy.c'],
                                include_dirs=['src'],
                                language='c',
                                extra_compile_args=[
                                    '-O3', '-fPIC', '-Wno-unknown-pragmas', '-std=c99'],
                                extra_link_args=extra_link_args
                                )


# Third-party modules - we depend on numpy for everything

# Obtain the numpy include directory.  This logic works across numpy versions.
try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()

extra_link_args = []
if sys.platform == 'darwin':
    from distutils import sysconfig
    vars = sysconfig.get_config_vars()
    vars['LDSHARED'] = vars['LDSHARED'].replace('-bundle', '-shared')
    extra_link_args = ['-Wl,-install_name,@rpath/_pyOrbfit' + suffix]

_pyOrbfitModule = Extension('_pyOrbfit',
                            ['src/orbfit/fit_radec.c',
                             'src/orbfit/pyOrbfit.i',
                             'src/orbfit/orbfit1.c',
                             'src/orbfit/nrutil.c',
                             'src/orbfit/ephem_earth.c',
                             'src/orbfit/aeiderivs.c',
                             'src/orbfit/gasdev.c',
                             'src/orbfit/abg_to_xyz.c',
                             'src/orbfit/gaussj.c',
                             'src/orbfit/orbfit2.c',
                             'src/orbfit/mrqmin_orbit.c',
                             'src/orbfit/abg_to_aei.c',
                             'src/orbfit/ludcmp.c',
                             'src/orbfit/dms.c',
                             'src/orbfit/covsrt.c',
                             'src/orbfit/ran1.c',
                             'src/orbfit/lubksb.c',
                             'src/orbfit/transforms.c',
                             'src/orbfit/mrqcof_orbit.c'],
                            include_dirs=[numpy_include],
                            language='c',
                            extra_compile_args=[
                                '-O3', '-Wno-implicit-function-declaration', '-Wno-unknown-pragmas'],
                            extra_link_args=extra_link_args
                            )


setup(
    name='spacerocks',
    version='1.1.1',
    description='A Python Package for Solar System Ephemerides and Dynamics.',
    author='Kevin J. Napier',
    author_email='kjnapier@umich.edu',
    url="https://github.com/kjnapier/spacerocks",
    packages=['spacerocks', 'spacerocks.data'],
    package_data={'spacerocks.data': ['*.csv', '*.dat', '*.423'], 'spacerocks.data.spice': ['*']},
    data_files=['spacerocks/data/observatories.csv',
                'spacerocks/data/spice/2000001.bsp',
                'spacerocks/data/spice/2000004.bsp',
                'spacerocks/data/spice/de440s.bsp',
                'spacerocks/data/spice/gm_de431.tpc',
                'spacerocks/data/spice/hst.bsp',
                'spacerocks/data/spice/latest_leapseconds.tls',
                'spacerocks/data/spice/nh.bsp',
                'spacerocks/data/observatories.dat',
                'spacerocks/data/binEphem.423'],
    include_package_data=True,
    install_requires=['healpy',
                      'numpy',
                      'spiceypy',
                      'astropy',
                      'pandas',
                      'rebound',
                      'reboundx'],
    ext_modules=[libspacerocksmodule, _pyOrbfitModule],
    zip_safe=False
)