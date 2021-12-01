![Alt text](../assets/logo.png)

`spacerocks` is a Python package that provides high-level abstractions 
for orbital dynamics and solar system observations. Its modern, 
expressive API makes it extremely easy to use.


```Python
from spacerocks import SpaceRock
from spacerocks.units import Units

units = Units()
units.angle = 'deg'

rock = SpaceRock(a=44, e=0.1, inc=10, node=140, arg=109, M=98, epoch='1 December 2021', origin='ssb', units=units)
```

You can then use the `propagate` method to propagate the rock to any epoch. This method used 
`rebound's` `ias15` integrator 

```Python
prop, planets, sim = rock.propagate(epochs=['2 December 2022'], model=2)
```

Here `prop` is a `SpaceRock` object containing the test particles, 
`planets` is a `SpaceRock` object containing the perturbers, and `sim` 
is the rebound simulation object.