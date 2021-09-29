![Alt text](assets/logo.png)

![spacerocks](https://github.com/kjnapes/spacerocks/workflows/spacerocks/badge.svg?branch=master)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Documentation Status](https://readthedocs.org/projects/spacerocks/badge/?version=latest)](https://spacerocks.readthedocs.io/en/latest/?badge=latest)
[![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/kjnapier/spacerocks.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/kjnapier/spacerocks/context:python)


SpaceRocks is a Python package for performing observational and theoretical tasks in solar system dynamics.

#### Currently implemented:
- Coordinate transformations for arbitrary solar system objects.
- Numerical propagation to future and past epochs.
- Precise ephemeride calculation.
- Transformations between barycentric and heliocentric coordinates.
- Outer solar system orbit fitting.
- Sky position predictions with uncertainties.

#### Under Development:

- General orbit fitting.
- Minor Planet Center Queries.


### Installation

The easiest way to install spacerocks is with pip:

```zsh
pip install spacerocks
```

If you want to be sure that you have the most up-to-date version, you can clone the repository and install from the source.

```zsh
git clone https://github.com/kjnapier/spacerocks.git
cd spacerocks
python setup.py build_ext -i
pip install .
```

The software is currently verified to be stable on the latest versions of macOS and Ubuntu. It is compatible with Python versions 3.6, and newer.
