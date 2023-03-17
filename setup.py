import os
from setuptools import setup, Extension
import sys
import sysconfig
import glob
import numpy
import subprocess

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
    #extra_link_args = ['-Wl,-install_name,@rpath/libspacerocks' + suffix]
    
    omp_path = subprocess.run(['brew', '--prefix', 'libomp'], stdout=subprocess.PIPE).stdout.decode("utf-8").split('\n')[0]
    llvm_path = subprocess.run(['brew', '--prefix', 'llvm'], stdout=subprocess.PIPE).stdout.decode("utf-8").split('\n')[0]
    extra_compile_args = ['-O3', '-fPIC', '-std=c++2a']

    extra_compile_args += [f'-I{llvm_path}/include', f'-I{omp_path}/include']

    from sysconfig import get_config_vars
    compiler = get_config_vars('CXX')
    if (compiler[0] == 'clang++') or (compiler[0] == 'clang'):
        extra_compile_args += ['-Xpreprocessor', '-fopenmp']
    else:
        extra_compile_args += ['-fopenmp']

libspacerocksmodule = Extension('libspacerocks',
                                sources=['src/spacerocks/pybindings.cpp', 
                                         'src/spacerocks/calc_E_from_M.cpp',
                                         'src/spacerocks/calc_E_from_f.cpp',
                                         'src/spacerocks/calc_f_from_E.cpp',
                                         'src/spacerocks/calc_M_from_E.cpp',
                                         'src/spacerocks/calc_kep_from_xyz.cpp',
                                         'src/spacerocks/calc_vovec_from_kep.cpp',
                                         'src/spacerocks/correct_for_ltt.cpp',
                                         'src/spacerocks/kepE_to_xyz.cpp',
                                         'src/spacerocks/kepM_to_xyz.cpp',
                                         'src/spacerocks/kepf_to_xyz.cpp', 
                                         'src/spacerocks/compute_topocentric_correction.cpp', 
                                         'src/spacerocks/compute_lst.cpp'],
                                include_dirs=['src/spacerocks'],
                                language='c++',
                                extra_compile_args=extra_compile_args,
                                extra_link_args=extra_link_args
                                )

extra_compile_args = ['-O', '-fPIC', '-std=c99', '-march=native']#, '-fno-stack-protector']
extra_link_args = []
if sys.platform == 'darwin':
    from distutils import sysconfig
    vars = sysconfig.get_config_vars()
    vars['LDSHARED'] = vars['LDSHARED'].replace('-bundle', '-shared')
    extra_link_args = ['-Wl,-lm,-install_name,@rpath/_pyOrbfit' + suffix]
    extra_compile_args = ['-O', '-fPIC', '-std=c99', '-w']#, '-fno-stack-protector']
    

_pyOrbfit = Extension('_pyOrbfit',
                      sources=['src/pyOrbfit/pyOrbfit.i', 
                               'src/pyOrbfit/fit_radec.c', 
                               'src/pyOrbfit/orbfit1.c', 
                               'src/pyOrbfit/nrutil.c', 
                               'src/pyOrbfit/ephem_earth.c', 
                               'src/pyOrbfit/aeiderivs.c', 
                               'src/pyOrbfit/gasdev.c', 
                               'src/pyOrbfit/gaussj.c',
                               'src/pyOrbfit/orbfit2.c', 
                               'src/pyOrbfit/mrqmin_orbit.c', 
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
    version='2.3.4',
    description='A Python Package for Solar System Ephemerides and Dynamics.',
    author='Kevin J. Napier',
    author_email='kjnapier@umich.edu',
    url="https://github.com/kjnapier/spacerocks",
    packages=['spacerocks', 'spacerocks.durin', 'spacerocks.survey', 
              'spacerocks.pyorbfit', 'spacerocks.orbfit', 'spacerocks.linking'],
    package_data={'spacerocks.data': ['observatories.csv'], 
                  'spacerocks.data.pyOrbfit': ['*']},
    data_files=data_files,
    include_package_data=True,
    install_requires=['rich',
                      'asdf',
                      'numpy',
                      'spiceypy',
                      'astropy',
                      'pandas',
                      'rebound >= 3.19.2', 
                      'astroquery'],
    ext_modules=[libspacerocksmodule, _pyOrbfit],
    zip_safe=False
)
