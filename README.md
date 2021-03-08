# SpaceRocks
![spacerocks](https://github.com/kjnapes/spacerocks/workflows/spacerocks/badge.svg?branch=master)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

SpaceRocks is an open-source software written in pure Python. It is a suite of tools for performing observational and theoretical tasks in solar system dynamics.

#### Currently implemented:

- Calculate the ephemerides of solar system objects on future and past dates.
- Transform between barycentric and heliocentric coordinates.

#### In the works:

- Positional uncertainty calculations.
- Orbit fitting.
- Minor Planet Center Queries.



### Obtaining and using the software

The easiest way to install spacerocks is with pip:

```zsh
pip install spacerocks
```

If you don't use pip, or if you want to be sure that you have the most up-to-date version, you can clone the repository and install from the source.

```zsh
git clone https://github.com/kjnapier/spacerocks.git
python setup.py install
```

The software is currently verified to be stable on macOS Mojave and Catalina (10.14 and 10.15), as well as the latest version of Ubuntu. It is compatible with Python versions 3.6, 3.7, and 3.8. It is intentionally ***not*** compatible with Python 2, which became deprecated on January 1, 2020.
