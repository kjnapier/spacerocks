from setuptools import setup, Extension

import os
import sys

import sysconfig
suffix = sysconfig.get_config_var('EXT_SUFFIX')
if suffix is None:
    suffix = ".so"

extra_link_args=[]
if sys.platform == 'darwin':
    from distutils import sysconfig
    vars = sysconfig.get_config_vars()
    vars['LDSHARED'] = vars['LDSHARED'].replace('-bundle', '-shared')
    extra_link_args=['-Wl,-install_name,@rpath/libspacerocks'+suffix]

libspacerocksmodule = Extension('libspacerocks',
                                sources=['src/speedy.c'],
                                include_dirs=['src'],
                                language='c',
                                extra_compile_args=['-O3', '-fPIC', '-Wno-unknown-pragmas', '-std=c99'],
                                extra_link_args=extra_link_args
                                )


# Third-party modules - we depend on numpy for everything
import numpy

# Obtain the numpy include directory.  This logic works across numpy versions.
try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()

extra_link_args=[]
if sys.platform == 'darwin':
    from distutils import sysconfig
    vars = sysconfig.get_config_vars()
    vars['LDSHARED'] = vars['LDSHARED'].replace('-bundle', '-shared')
    extra_link_args=['-Wl,-install_name,@rpath/_pyOrbfit'+suffix]

_pyOrbfitModule = Extension('_pyOrbfit',
                            ['src/fit_radec.c',
                             'src/pyOrbfit.i',
                             'src/orbfit1.c',
                             'src/nrutil.c',
                             'src/ephem_earth.c',
                             'src/aeiderivs.c',
                             'src/gasdev.c',
                             'src/abg_to_xyz.c',
                             'src/gaussj.c',
                             'src/orbfit2.c',
                             'src/mrqmin_orbit.c',
                             'src/abg_to_aei.c',
                             'src/ludcmp.c',
                             'src/dms.c',
                             'src/covsrt.c',
                             'src/ran1.c',
                             'src/lubksb.c',
                             'src/transforms.c',
                             'src/mrqcof_orbit.c'],
                             include_dirs=[numpy_include],
                             language='c',
                             extra_compile_args=['-O3', '-Wno-implicit-function-declaration'],
                             extra_link_args=extra_link_args
                            )


setup(
   name='spacerocks',
   version='1.0.23',
   description='A Python Package for Solar System Ephemerides and Dynamics.',
   author='Kevin J. Napier',
   author_email='kjnapier@umich.edu',
   url="https://github.com/kjnapier/spacerocks",
   packages=['spacerocks', 'spacerocks.data'],
   package_data={'spacerocks.data': ['*.csv'], 'spacerocks.data': ['*.dat', '*.423']},
   data_files=['spacerocks/data/observatories.csv', 'spacerocks/data/observatories.dat', 'spacerocks/data/binEphem.423'],
   include_package_data=True,
   install_requires=['healpy',
                     'numpy',
                     'skyfield',
                     'astropy',
                     'pandas',
                     'rebound',
                     'reboundx'],
   ext_modules=[libspacerocksmodule, _pyOrbfitModule],
   zip_safe=False
)
