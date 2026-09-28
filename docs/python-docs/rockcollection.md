<h1 style="border-bottom: 5px solid white;">RockCollection Module</h1>


### Table of Contents
1. [Overview](#overview)
2. [Primary Classes](#primary-classes)
3. [Methods](#methods)
4. [Properties](#properties)
5. [Examples](#examples)
6. [Notes](#notes)

<h2 style="border-bottom: 3px solid white;">Overview</h2>


The RockCollection module provides tools for managing and manipulating groups of SpaceRock objects. It supports batch operations like observing multiple objects, filtering collections, and fetching orbital data from external sources such as the Minor Planet Center (MPC).

<h2 style="border-bottom: 3px solid white;">Primary Classes</h2>


### RockCollection
```python
class RockCollection:
    """
    Manages collections of SpaceRock objects with batch operations support.
    
    Example:
        collection = RockCollection()
        collection.add(rock)
        observations = collection.observe(observer)
    """
```

<h2 style="border-bottom: 3px solid white;">Methods</h2>


### Constructor Methods
---
**`new()`**

**Returns:**
- New RockCollection instance

*Example:*
```python
collection = RockCollection()
```

**`from_mpc()`**
```python
@classmethod
def from_mpc(cls, mpc_path: str, catalog: str, 
             download_data: bool = True) -> RockCollection
```

**Arguments:**
- `mpc_path`: Directory for MPC data storage
- `catalog`: Name of MPC catalog (e.g., 'mpcorb_extended')
- `download_data`: Whether to download missing data from MPC
- `orbit_type`: Specify only certain types of objects, i.e. "Atira"

**Returns:**
- RockCollection populated with MPC data

*Example:*
```python
trojan_collection = RockCollection.from_mpc(
    mpc_path="/path/to/store/data",
    catalog="mpcorb_extended",
    download_data=True, 
    orbit_type="Jupiter Trojan"
)
```

### Operation Methods
---

**`propagate()`**
```python
def propagate(self, epoch: Time, kernel: SpiceKernel, method: str = "nbody", chunk_size: int = 64) -> None
```
Moves every rock to `epoch`. With `method="nbody"` (default) the rocks are integrated with IAS15
in the field of the Sun, planets, Moon, Pluto and the 16 most massive asteroids from `kernel`;
`method="twobody"` uses Keplerian motion instead.

Rocks are integrated together, up to `chunk_size` per simulation: they are grouped by epoch and
sorted by perihelion distance, so each simulation holds dynamically similar objects, and
simulations run in parallel with the GIL released. Sharing a simulation means the perturber
ephemeris is evaluated once per step for the whole group, which makes this 3–6x faster than
integrating rocks one at a time. Rocks in a group share its step size, which is never longer than
any member's own, so results differ from `SpaceRock.propagate` only at the integrator's tolerance
(typically metres after months) and are at least as accurate. `chunk_size=1` integrates every
rock on its own, exactly like `SpaceRock.propagate`.

Each rock keeps its reference plane and origin (SUN or SSB; a custom origin is returned as SSB).
Raises `ValueError` if a rock cannot be integrated (for example, when the kernel does not cover the
requested epoch).

*Example:*
```python
kernel = SpiceKernel.defaults()
rocks.propagate(Time(2461000.5, "tdb", "jd"), kernel)
```

**`ephemeris()`**
```python
def ephemeris(self, epochs, observer, kernel: SpiceKernel, method: str = "nbody",
              timescale: str = "utc", return_states: bool = False, chunk_size: int = 64) -> dict
```
**Arguments:**
- `epochs`: a list of `Time` objects, or an array of Julian dates in `timescale`
- `observer`: an `Observatory` (its position is computed at each epoch), or a list of
  `Observer` objects, one per epoch (then `epochs` may be `None`)
- `kernel`: the SPICE kernel (perturbers and observer positions)
- `method`: `"nbody"` or `"twobody"`

**Returns:** a dict of NumPy arrays of shape `(len(rocks), n_epochs)`:

| Key | Units | Description |
|---|---|---|
| `ra`, `dec` | radians | Astrometric position (J2000 equator), corrected for light travel time |
| `ra_rate`, `dec_rate` | radians/day | Rates (`ra_rate` is dRA/dt, not multiplied by cos(dec)) |
| `range`, `range_rate` | AU, AU/day | Observer–object distance and its rate |
| `r_helio` | AU | Sun–object distance |
| `phase`, `elong` | radians | Sun–object–observer and Sun–observer–object angles |
| `mag` | mag | H-G apparent magnitude (NaN for rocks without an absolute magnitude) |
| `epoch` | TDB JD | The epochs, shape `(n_epochs,)` |
| `states` | AU, AU/day | Barycentric J2000 states, shape `(len(rocks), n_epochs, 6)` (only with `return_states=True`) |

Each rock is integrated once through all epochs, which can be in any order and before or after the
rocks' epochs, and states are read off the integrator's dense output between steps. This is
typically 40x or more faster than propagating and observing at each epoch. The collection is not
modified.

*Example:*
```python
import numpy as np
kernel = SpiceKernel.defaults()
w84 = Observatory.from_obscode("W84")
nights = np.arange(2461000.6, 2461030.6, 1.0)      # UTC Julian dates
eph = rocks.ephemeris(nights, w84, kernel)
visible = (eph["mag"] < 24) & (np.degrees(eph["elong"]) > 90)
```

**`analytic_propagate()`**
```python
def analytic_propagate(self, epoch: Time) -> RockCollection
```
**Arguments:**
- `epoch`: The Time object representing the epoch that we want to propagate the SpaceRocks in our RockCollection to.

Calculates and converts SpaceRock elements at new epoch, in place, for all rocks in the RockCollection.

*Example:*
```python
rockcollection = RockCollection()
# ...
ten_years = Time.now() + (365.25 * 10)
rockcollection.analytic_propagate(ten_years)
```

**`add()`**
```python
def add(self, rock: SpaceRock) -> None
```

**Arguments:**
- `rock`: SpaceRock object to add

*Example:*
```python
collection = RockCollection()
rock = SpaceRock.from_horizons("Ceres", epoch, "ECLIPJ2000", "SSB")
collection.add(rock)
```

**`observe()`**
```python
def observe(self, observer: Observer) -> List[Observation]
```

**Arguments:**
- `observer`: Observer object representing viewing location

**Returns:**
- List of Observation objects

*Example:*
```python
observatory = Observatory.from_obscode("F51")
observer = observatory.at(epoch)

observations = collection.observe(observer)
```

**`observe_arrays()`**
```python
def observe_arrays(self, observer: Observer) -> dict
```
The same computation as `observe`, but returns a dict of 1-D NumPy arrays (one entry per rock;
same keys as `ephemeris` without `epoch` and `states`) instead of a list of `Observation` objects.
For large collections this is roughly 7x faster, and much faster again than pulling values out of
`Observation` objects one by one. The rocks and the observer must be at the same epoch and in the
same reference plane.

*Example:*
```python
obs = collection.observe_arrays(observatory.at(epoch, kernel, "J2000", "SSB"))
ra, dec, mag = obs["ra"], obs["dec"], obs["mag"]
```

**`filter()`**


```python
def filter(self, indices: List[bool]) -> RockCollection
```

**Arguments**
- indices: List of boolean values. Must have the same length as the number of rocks.

**Returns**
- A new RockCollection containing only selected rocks

*Example*
```python
collection_under_5_au = collection.filter(collection.a() < 5)
```

**`change_reference_plane()`**
```python
def change_reference_plane(self, reference_plane: str) -> None
```

**Arguments:**
- `reference_plane`: New reference plane identifier

*Example:*
```python
collection.change_reference_plane("ECLIPJ2000")
```

<h2 style="border-bottom: 3px solid white;">Properties</h2>


### Position and Velocity Properties
| Property | Type | Description |
|----------|------|-------------|
| `x`, `y`, `z` | `numpy.ndarray` | Position components in AU |
| `vx`, `vy`, `vz` | `numpy.ndarray` | Velocity components in AU/day |
| `epoch` | `List[Time]` | Epochs of all objects |
| `rocks` | `List[SpaceRock]` | All SpaceRock objects |

<h2 style="border-bottom: 3px solid white;">Examples</h2>


### Basic Usage
```python
from spacerocks import RockCollection, SpaceRock
from spacerocks.time import Time

# Create collection and add objects
collection = RockCollection()
epoch = Time.now()

ceres = SpaceRock.from_horizons("Ceres", epoch, "ECLIPJ2000", "SSB")
vesta = SpaceRock.from_horizons("Vesta", epoch, "ECLIPJ2000", "SSB")

collection.add(ceres)
collection.add(vesta)

print(f"Collection size: {len(collection.rocks)}")
```

### Observing Multiple Objects
```python
from spacerocks.observing import Observatory

# Set up observer
observatory = Observatory.from_obscode('F51')
observer = observatory.at(epoch)

# Get observations for all objects
observations = collection.observe(observer)

# Process observations
for rock, obs in zip(collection.rocks, observations):
    print(f"{rock.name}:")
    print(f"  RA: {obs.ra}")
    print(f"  Dec: {obs.dec}")
    print(f"  Magnitude: {obs.mag}")
```
### Filtering Rocks
```python
# Filter rocks with eccentricity greater than 0.1
mask = [rock.e() > 0.05 for rock in collection]
filtered_collection = collection.filter(mask)

print(f"Filtered collection size: {len(filtered_collection.rocks)}")
```

### Loading from MPC
```python
# Create collection from MPC data
collection = RockCollection.from_mpc(
    mpc_path="~/mpc_data",
    catalog="mpcorb_extended",
    download_data=True
)

# Change reference plane for all objects
collection.change_reference_plane("J2000")
```

<h2 style="border-bottom: 3px solid white;">Notes</h2>


- Collections maintain the order of added objects
- All position values are in AU
- All velocity values are in AU/day
- MPC data is stored locally when downloaded
- Compatible with all SpaceRock instantiation methods
- Empty collections can be created with RockCollection()

