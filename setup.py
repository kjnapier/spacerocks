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


setup(
   name='spacerocks',
   version='1.0.23',
   description='A Python Package for Solar System Ephemerides and Dynamics.',
   author='Kevin J. Napier',
   author_email='kjnapier@umich.edu',
   url="https://github.com/kjnapier/spacerocks",
   packages=['spacerocks', 'spacerocks.data'],
   package_data={'spacerocks.data': ['*.csv']},
   data_files=['spacerocks/data/observatories.csv'],
   include_package_data=True,
   # I'm not sure package_data is working.
   install_requires=['healpy',
                     'numpy',
                     'skyfield',
                     'astropy',
                     'pandas',
                     'rebound',
                     'reboundx'],
   ext_modules=[libspacerocksmodule],
   zip_safe=False
)
