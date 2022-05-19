import os
from setuptools import setup, Extension
import sys
import sysconfig
import glob
import numpy

suffix = sysconfig.get_config_var('EXT_SUFFIX')
if suffix is None:
    suffix = ".so"

extra_link_args = ['-lgomp']
#extra_compile_args = ['-O3', '-fPIC', '-std=gnu++2a', '-march=native', '-fopenmp']
extra_compile_args = ['-O3', '-fPIC', '-march=native', '-fopenmp']

if sys.platform == 'darwin':
    from distutils import sysconfig
    vars = sysconfig.get_config_vars()
    vars['LDSHARED'] = vars['LDSHARED'].replace('-bundle', '-shared')
    extra_link_args = ['-Wl,-lomp,-install_name,@rpath/libspacerocks' + suffix]
    extra_compile_args = ['-O3', '-fPIC', '-std=c++2a', '-march=native', '-Xclang', '-fopenmp']
    

libspacerocksmodule = Extension('libspacerocks',
                                sources=['src/speedy.cpp'],
                                include_dirs=['src'],
                                language='c++',
                                extra_compile_args=extra_compile_args,
                                extra_link_args=extra_link_args
                                )

extra_compile_args = ['-O3', '-fPIC', '-std=c99', '-march=native', '-fno-stack-protector']
extra_link_args = []
if sys.platform == 'darwin':
    from distutils import sysconfig
    vars = sysconfig.get_config_vars()
    vars['LDSHARED'] = vars['LDSHARED'].replace('-bundle', '-shared')
    extra_link_args = ['-Wl,-install_name,@rpath/_pyOrbfit' + suffix]
    extra_compile_args = ['-O3', '-fPIC', '-std=c99', '-march=native', '-w', '-fno-stack-protector']
    

_pyOrbfit = Extension('_pyOrbfit',
                      sources=['src/pyOrbfit/pyOrbfit.i', 
                               'src/pyOrbfit/fit_radec.c', 
                               'src/pyOrbfit/orbfit1.c', 
                               'src/pyOrbfit/nrutil.c', 
                               'src/pyOrbfit/ephem_earth.c', 
                               'src/pyOrbfit/aeiderivs.c', 
                               'src/pyOrbfit/gasdev.c', 
                               'src/pyOrbfit/abg_to_xyz.c', 
                               'src/pyOrbfit/gaussj.c',
                               'src/pyOrbfit/orbfit2.c', 
                               'src/pyOrbfit/mrqmin_orbit.c', 
                               'src/pyOrbfit/abg_to_aei.c', 
                               'src/pyOrbfit/ludcmp.c', 
                               'src/pyOrbfit/dms.c', 
                               'src/pyOrbfit/covsrt.c', 
                               'src/pyOrbfit/ran1.c',
                               'src/pyOrbfit/lubksb.c', 
                               'src/pyOrbfit/transforms.c', 
                               'src/pyOrbfit/mrqcof_orbit.c'],
                      include_dirs=['src/pyOrbfit', numpy.get_include()],
                      language='c',
                      extra_compile_args=extra_compile_args,
                      extra_link_args=extra_link_args
                      )

data_files = []
dirs = ['spacerocks/data/pyOrbfit/*', 
        'spacerocks/data/observatories.csv', 
        'src/pyOrbfit/*']

for dir in dirs:
   for filename in glob.glob(dir):
      if os.path.isfile(filename):
         data_files.append(filename)

setup(
    name='spacerocks',
    version='2.1.11',
    description='A Python Package for Solar System Ephemerides and Dynamics.',
    author='Kevin J. Napier',
    author_email='kjnapier@umich.edu',
    url="https://github.com/kjnapier/spacerocks",
    packages=['spacerocks'],
    # package_data={'spacerocks.data': ['*.csv'], 
    #               'spacerocks.data.spice': ['*'], 
    #               'spacerocks.data.pyOrbfit': ['*']},
    package_data={'spacerocks.data': ['observatories.csv'], 
                  'spacerocks.data.pyOrbfit': ['*']},
    data_files=data_files,
    include_package_data=True,
    install_requires=['healpy',
                      'rich',
                      'asdf',
                      'numpy',
                      'spiceypy',
                      'astropy',
                      'pandas',
                      'rebound', 
                      'astroquery'],
    ext_modules=[libspacerocksmodule, _pyOrbfit],
    zip_safe=False
)
