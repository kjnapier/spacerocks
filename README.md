![Alt text](assets/logo.png)

![spacerocks](https://github.com/kjnapes/spacerocks/workflows/spacerocks/badge.svg?branch=master)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Documentation Status](https://readthedocs.org/projects/spacerocks/badge/?version=latest)](https://spacerocks.readthedocs.io/en/latest/?badge=latest)
[![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/kjnapier/spacerocks.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/kjnapier/spacerocks/context:python)


`spacerocks` is a Python package that provides high-level abstractions 
for orbital dynamics and solar system observations. Its modern, 
expressive API makes it extremely easy to use.

The primary data structure in `spacerocks` is a class called `SpaceRock`. 
You can instantiate a `SpaceRock` object using any valid set of 6 Keplerian 
elements, or a state vector.

```Python
from spacerocks import SpaceRock
from spacerocks.units import Units

units = Units()
units.timescale = 'utc'

rock = SpaceRock(a=44, 
                 e=0.1, 
                 inc=10, 
                 node=140, 
                 arg=109, 
                 M=98, 
                 H=7, 
                 epoch='1 December 2021', 
                 origin='ssb', 
                 frame='eclipJ2000', 
                 units=units)
```
The `Units` object we instantiated is extremely useful for avoiding bugs, 
as it allows for an explicit set of units to be specified only once. 
You can print the current units with `units.current()`, and 
you can set the individual attributes using either strings or astropy units.
We advise using explicit astropy units where applicable 
(i.e. `units.angle = u.deg` rather than `units.angle = 'deg'`). 
Finally, note that you can also pass a JD or an MJD as an epoch. 
If you provide a string (in any format) the program will use `dateutil` 
to try to parse the date.

You can then access the attributes of your `SpaceRock` object. 
The attributes are lazily computed in the interest of efficiency. 
The attributed all carry `astropy` units, with angles stored as `Angle` 
objects, distances stored as `Distance` objects, and times are 
stored as `Time` objects. This encures a complete lack of ambiguity, 
and very easy unit conversions.

`SpaceRock` objects are vectorized, allowing for the processing of multiple objects at once. 

```Python
rocks = SpaceRock(a=[44, 45], 
                  e=[0.1, 1.2], 
                  inc=[10, 170], 
                  node=[140, 303], 
                  arg=[109, 23], 
                  f=[98, 124], 
                  H=[7, 8.2],
                  epoch=['1 December 2021', '3 December 2021'], 
                  origin='ssb', 
                  units=units)
```
Notice that a `SpaceRock` object can handle both elliptical and hyperbolic orbits 
simultaneously (though parabolic orbits are not yet supported), and it can handle 
nonuniform epochs. 

`SpaceRock` objects have a number of other utilities. 
First, you can slice `SpaceRock` objects in a Pythonic way
```Python
r = rocks[rocks.e < 1]
```

You can easily change the origin of the coordinate system. This is particularly 
useful for converting the Minor Planet Center's heliocentric elements to 
barycentric elements. 
```Python
# convert to heliocentric coordinates
rocks.to_helio()

# convert back to barycentric coordinates
rocks.to_bary()

# convert to New Horizons-centric coordinates
rocks.change_origin(spiceid=-98)
```

You can use the `propagate` method to propagate the rocks to any epochs. 
This method uses `rebound's` `ias15` integrator under the hood, and automatically 
synchronizes the epochs of the rocks so you don't have to. Notice that we are 
specifying the epochs to be in Barycentric Dynamical Time.

```Python
units = Units()
units.timescale = 'tdb'

prop, planets, sim = rock.propagate(epochs=['2 December 2021', '4 December 2021'], model=2, units=units)
```

Here `prop` is a `SpaceRock` object containing the test particles at 
all of the epochs, `planets` is a `SpaceRock` object containing the 
perturbers at all of the epochs, and `sim` is the rebound simulation 
object at the final epoch. The `model` argument sets the perturbers as follows.

| model | Perturbers                                                               |
|:-----:|:-------------------------------------------------------------------------|
|   0   | Sun                                                                      |
|   1   | Sun, Jupiter, Saturn, Uranus, Neptune                                    |
|   2   | Sun, Mercury, Venus, Earth, Moon, Mars, Jupiter, Saturn, Uranus, Neptune |
|   3   | Full set of JPL Horizons perturbers                                      |

You can use the `observe` method compute the objects' ephemerides from an 
arbitrary location in the solar system. Here we compute the ephemerides of
our rocks from DECam.

```Python
obs = rocks.observe(obscode='W84')
```

This method returns an `Ephemerides` object which contains the rocks' state 
vectors with respect to the observer, corrected for light travel time. 
These values allow us to compute the objects' observable properties, which 
are accessible as attriutes to the `Ephemerides` object.

Finally, you can write and read `SpaceRock` objects to and from `asdf` files.
```Python
rocks.to_file('rocks.rocks')

rocks_from_disk = SpaceRock.from_file('rocks.rocks')
```

For more information, see [this link](../cods/spacerocks.md)

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
