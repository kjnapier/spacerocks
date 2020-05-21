from setuptools import setup

setup(
   name='spacerocks',
   version='0.6.3',
   description='Calculate solar system ephemerides from orbital elements.',
   author='Kevin Napier',
   author_email='kjnapier@umich.edu',
   url="https://github.com/kjnapes/spacerocks",
   packages=['spacerocks'],
   # I'm not sure package_data is working.
   package_data={'spacerocks.data': ['observatories.csv']},
   install_requires=['healpy', 'numpy', 'matplotlib', 'scipy', 'skyfield',
                     'astropy', 'numba', 'pandas'],
   include_package_data=True
)
