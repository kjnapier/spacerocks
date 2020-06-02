# SpaceRocks
![spacerocks](https://github.com/kjnapes/spacerocks/workflows/spacerocks/badge.svg?branch=master)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

SpaceRocks is an open-source software written in pure Python. It is a suite of tools useful for performing a number of tasks relevant to solar system dynamics, both observational and theoretical. 

#### Currently implemented:

- Calculate the ephemerides of solar system objects on future and past dates.
- Transform between barycentric and heliocentric coordinates.

#### In the works:

- Positional uncertainty calculations.
- Orbit fitting.



### Obtaining and using the software

The easiest way to install spacerocks is with pip:

```zsh
pip install spacerocks
```

If you don't use pip, or if you want to be sure that you have the most up-to-date version, you can clone the repository and install from the source.

```zsh
git clone https://github.com/kjnapes/spacerocks.git
python setup.py install
```

The software is currently verified to be stable on macOS Mojave and Catalina (10.14 and 10.15), as well as the latest version of Ubuntu. It is compatible with Python versions 3.6, 3.7, and 3.8. It is intentionally ***not*** compatible with Python 2, which became deprecated on January 1, 2020.

SpaceRocks has an optinal plotting dependency called `cartopy`. Cartopy is the successor to the `Basemap` package, which became deprecated along with Python 2. On macOS Catalina (10.15), you can successfully install `cartopy` with the following recipe. On macOS Mojave (10.14), you should remove the prefix `CFLAGS=‘-stdlib=libc++’` on the final line.

```zsh
pip3 uninstall shapely
brew install proj geos
pip3 install —upgrade cython numpy pyshp six
pip3 install shapely —no-binary shapely
CFLAGS=‘-stdlib=libc++’ pip install cartopy
```
