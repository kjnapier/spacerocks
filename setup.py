from setuptools import setup, Extension

from codecs import open
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
                                sources = ['src/speedy.cpp'],
                                include_dirs = ['src'],
                                #define_macros=[ ('LIBSPACEROCKS', None) ],
                                extra_compile_args=['-O3', '-fPIC'],
                                extra_link_args=extra_link_args
                                )


setup(
   name='spacerocks',
   version='1.0.7',
   description='A Python Package for Solar System Ephemerides and Dynamics.',
   author='Kevin J. Napier',
   author_email='kjnapier@umich.edu',
   url="https://github.com/kjnapier/spacerocks",
   packages=['spacerocks'],
   # I'm not sure package_data is working.
   package_data={'spacerocks': ['data/observatories.csv', 'sr_cpp.so']},
   install_requires=['healpy',
                     'numpy',
                     'skyfield',
                     'astropy',
                     'pandas',
                     'rebound',
                     'reboundx'],
   include_package_data=True
   ext_modules = [libspacerocksmodule]
)
