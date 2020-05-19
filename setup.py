from setuptools import setup

setup(
   name='spacerocks',
   version='0.5.2',
   description='Calculate solar system ephemerides from orbital elements.',
   author='Kevin Napier',
   author_email='kjnapier@umich.edu',
   url="https://github.com/kjnapes/spacerocks",
   packages=['spacerocks'],
   install_requires=['healpy', 'numpy', 'matplotlib', 'scipy', 'skyfield',
                     'astropy', 'numba', 'pandas']
)
