# `Origin`

The coordinates of any body must be specified with respect to some origin. In orbital mechanics, specifically the two-body problem, the origin is usually the sun or the solar system barycenter. For planetary satellites, the origin is the center of the planet or the center of mass of the planet-satellite system. `spacerocks` keeps track of the origin of a body using the `Origin` enum. `Origin` contains a string that represents the name of the origin, and a quantity called `mu`, which is the gravitational parameter of the origin. 

`spacerocks` has two built-in origins (sun and solar system barycenter), and supports the creation of custom `Origin` instances.

#### Examples
```python
from spacerocks import Origin

# Sun is built-in
sun = Origin.sun()

# Solar System Barycenter is built-in
ssb = Origin.ssb()

# If you have a different use case, you can create an Origin instance as
mu = 0.000314159
custom = Origin("My Origin", mu)
```

#### Technical Details
`Origin` is implemented under the hood as a `rust` enum so that the built-in cases don't need to carry around a string to specify the name. Cloning and moving strings is expensive, so by implementing `Origin` as a `rust` enum, we can avoid the overhead of cloning strings in the most common use cases. All of this is abstracted away from the user, so you can use `Origin` like any other Python class.