from setuptools import setup

setup(
   name='spacerocks',
   version='0.7.0',
   description='Calculate solar system ephemerides from orbital elements.',
   author='Kevin Napier',
   author_email='kjnapier@umich.edu',
   url="https://github.com/kjnapes/spacerocks",
   packages=['spacerocks'],
   # I'm not sure package_data is working.
   package_data={'spacerocks': ['data/observatories.csv']},
   install_requires=['healpy', 'numpy', 'matplotlib', 'skyfield',
                     'astropy', 'numba', 'pandas', 'rebound',
                     'reboundx', 'ephem'],
   include_package_data=True
)
