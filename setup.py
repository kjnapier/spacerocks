import os
import numpy
from setuptools import setup, Extension
import sys
import sysconfig
import glob

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

data_files = []
dirs = ['spacerocks/data/spice/*', 'spacerocks/data/spice/asteroids/*', 'spacerocks/data/*']
for dir in dirs:
   for filename in glob.glob(dir):
      if os.path.isfile(filename):
         data_files.append(filename)

setup(
    name='spacerocks',
    version='2.0.0',
    description='A Python Package for Solar System Ephemerides and Dynamics.',
    author='Kevin J. Napier',
    author_email='kjnapier@umich.edu',
    url="https://github.com/kjnapier/spacerocks",
    packages=['spacerocks'], #, 'spacerocks.data'],
    package_data={'spacerocks.data': ['*.csv'], 
                  'spacerocks.data.spice': ['*'], 
                  'spacerocks.data.spice.asteroids': ['*.bsp']},
    data_files=data_files,
    include_package_data=True,
    install_requires=['healpy',
                      'asdf',
                      'numpy',
                      'spiceypy',
                      'astropy',
                      'pandas',
                      'rebound'],
    ext_modules=[libspacerocksmodule],
    zip_safe=False
)