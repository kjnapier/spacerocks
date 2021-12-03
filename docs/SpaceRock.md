# SpaceRock

The primary data structure in `spacerocks` is a class called `SpaceRock`. 
You can instantiate a `SpaceRock` object using any valid set of 6 Keplerian 
elements, or a state vector.

## Instantiation

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
You can read more about it [here](./Units.md).

Note that you can also pass a JD or an MJD as an epoch. 
If you provide a string (in any format) the program will use `dateutil` 
to try to parse the date.

You can then access the attributes of your `SpaceRock` object. 
The attributes are lazily computed in the interest of efficiency. 
The attributed all carry `astropy` units, with angles stored as `Angle` 
objects, distances stored as `Distance` objects, and times are 
stored as `Time` objects. This encures a complete lack of ambiguity, 
and very easy unit conversions.

| Orbital Parameter                      | Attribute    |
|:---------------------------------------|:-------------|
| semi-major axis                        | a            |
| eccentricity                           | e            |
| inclination                            | inc          |
| longitude of ascending node            | node         |
| longitude of pericenter                | varpi        |
| argument of pericenter                 | arg          |
| mean anomaly                           | M            |
| true anomaly                           | f            |
| eccentric anomaly                      | E            |
| semi-minor axis                        | b            |
| pericenter distance                    | q            |
| apocenter distance                     | Q            |


| Metadata                               | Attribute    |
|:---------------------------------------|:-------------|
| name                                   | name         |


| Physical Property                      | Attribute  | 
|:---------------------------------------|:-----------|
| absolute magnitude                     | H          |
| phase-slope constant                   | G          |
| albedo                                 | albedo     |
| mass                                   | mass       |
| radius                                 | radius     |
| diameter                               | diameter   |


## Vectorization

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

## Object Slicing

`SpaceRock` objects have a number of other utilities. 
First, you can slice `SpaceRock` objects in a Pythonic way
```Python
r = rocks[rocks.e < 1]
```

## Coordinate Transformations

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

## The `propagate` Method

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

## The `observe` Method

```Python
obs = rocks.observe(obscode='W84')
```

This method returns an [Ephemerides](./Ephemerides.md) object which contains the rocks' state 
vectors with respect to the observer, corrected for light travel time. 
These values allow us to compute the objects' observable properties, which 
are accessible as attriutes to the `Ephemerides` object.

Finally, you can write and read `SpaceRock` objects to and from `asdf` files.
```Python
rocks.to_file('rocks.rocks')

rocks_from_disk = SpaceRock.from_file('rocks.rocks')
```








