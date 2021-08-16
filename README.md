# SpaceRocks
![spacerocks](https://github.com/kjnapes/spacerocks/workflows/spacerocks/badge.svg?branch=master)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Documentation Status](https://readthedocs.org/projects/spacerocks/badge/?version=latest)](https://spacerocks.readthedocs.io/en/latest/?badge=latest)

SpaceRocks is an open-source software written in Python. The most expensive operations are written in C, and then wrapped with Python. It is a suite of tools for performing observational and theoretical tasks in solar system dynamics.

#### Currently implemented:
- Coordinate transformations for arbitrary solar system objects.
- Numerical propagation to future and past epochs.
- Precise ephemeride calculation.
- Transformations between barycentric and heliocentric coordinates.
- Outer solar system orbit fitting.
- Sky position predictions with uncertainties.

#### In the works:

- General orbit fitting.
- Minor Planet Center Queries.



### Obtaining and using the software

The easiest way to install spacerocks is with pip:

```zsh
pip install spacerocks
```

If you don't use pip, or if you want to be sure that you have the most up-to-date version, you can clone the repository and install from the source.

```zsh
git clone https://github.com/kjnapier/spacerocks.git
cd spacerocks
python setup.py build_ext -i
pip install .
```

The software is currently verified to be stable on the latest versions of macOS and Ubuntu. It is compatible with Python versions 3.6, and newer. It is intentionally ***not*** compatible with Python 2, which became deprecated on January 1, 2020.
