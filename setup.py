from setuptools import setup

setup(
   name='spacerocks',
   version='1.0.7',
   description='A Python Package for Solar System Ephemerides and Dynamics.',
   author='Kevin J. Napier',
   author_email='kjnapier@umich.edu',
   url="https://github.com/kjnapier/spacerocks",
   packages=['spacerocks'],
   # I'm not sure package_data is working.
   package_data={'spacerocks': ['data/observatories.csv']},
   install_requires=['healpy',
                     'numpy',
                     'skyfield',
                     'astropy',
                     'pandas',
                     'rebound',
                     'reboundx'],
   include_package_data=True
)
