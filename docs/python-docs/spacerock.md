<h1 style="border-bottom: 5px solid white;">SpaceRock Module</h1>

### Table of Contents
1. [Overview](#overview)
2. [Primary Classes](#primary-classes)
3. [Constructor Methods](#constructor-methods)
4. [Operation Methods](#operation-methods)
5. [Derived Methods](#derived-methods)
6. [Setter Methods](#setter-methods)
7. [Getter Methods](#getter-methods)
8. [Examples](#examples)
9. [Notes](#notes)

<h2 style="border-bottom: 3px solid white;">Overview</h2>
The SpaceRock module provides functionality for representing and manipulating celestial bodies in spacerocks. It handles object state vectors, orbital elements, physical properties, and observational calculations.

<h2 style="border-bottom: 3px solid white;">Primary Classes</h2>


### SpaceRock
```python
class SpaceRock:
    """
    Core class for representing celestial bodies in spacerocks.
    
    Example:
        rock = SpaceRock.from_horizons("Arrokoth", epoch)
        rock.set_mass(1e-10)  
    """
```

<h2 style="border-bottom: 3px solid white;">Constructor Methods</h2>


**`from_horizons()`**
```python
@classmethod
def from_horizons(cls, name: str, epoch: Time, reference_plane: str = "ECLIPJ2000", 
                 origin: str = "SSB") -> SpaceRock
```
**Arguments:**
- `name`: Object identifier for JPL Horizons
- `epoch`: Time object representing the state epoch
- `reference_plane`: Reference frame (e.g., "ECLIPJ2000")
- `origin`: Origin of coordinate system (e.g., "SSB")

**Returns:**
- New SpaceRock instance from Horizons data

*Example:*
```python
rock = SpaceRock.from_horizons("Ceres", epoch=Time.now())
```

**`from_spice()`**
```python
@classmethod
def from_spice(cls, name: str, epoch: Time, reference_plane: str = "ECLIPJ2000", 
               origin: str = "SSB") -> SpaceRock
```
**Arguments:**
- `name`: Object name in loaded SPICE kernel
- `epoch`: Time object representing the state epoch
- `reference_plane`: Reference frame (e.g., "ECLIPJ2000")
- `origin`: Origin of coordinate system (e.g., "SSB")

**Returns:**
- New SpaceRock instance from SPICE data

*Example:*
```python
jupiter = SpaceRock.from_spice("jupiter barycenter", epoch=Time.now())
```

**`from_xyz()`**
```python
@classmethod
def from_xyz(cls, name: str, x: float, y: float, z: float, vx: float, vy: float, 
            vz: float, epoch: Time, reference_plane: str = "ECLIPJ2000", 
            origin: str = "SSB") -> SpaceRock
```
**Arguments:**
- `name`: Object identifier
- `x, y, z`: Position components in AU
- `vx, vy, vz`: Velocity components in AU/day
- `epoch`: Time object representing the state epoch
- `reference_plane`: Reference frame (e.g., "ECLIPJ2000")
- `origin`: Origin of coordinate system (e.g., "SSB")

**Returns:**
- New SpaceRock instance from Cartesian coordinates

*Example:*
```python
rock = SpaceRock.from_xyz("asteroid", x=1.0, y=0.0, z=0.0, 
                         vx=0.0, vy=1.0, vz=0.0, epoch=epoch)
```

**`from_spherical()`**
```python
@classmethod
def from_spherical(cls, name: str, phi: float, theta: float, r: float, vr: float, 
                  vo: float, psi: float, epoch: Time, reference_plane: str = "ECLIPJ2000",
                  origin: str = "SSB") -> SpaceRock
```
***Note:***
This is intended to be used with a custom spherical coordinate basis (Napier and Holman (2024)). The foundation of this basis can be found [here](https://arxiv.org/abs/2410.03598).

**Arguments:**
- `name`: Object identifier
- `phi`: Longitude in radians
- `theta`: Latitude in radians
- `r`: Distance in AU
- `vr`: Radial velocity in AU/day
- `vo`: Tangential velocity in AU/day
- `psi`: Angle between radial and tangential velocities in radians
- `epoch`: Time object representing the state epoch
- `reference_plane`: Reference frame (e.g., "ECLIPJ2000")
- `origin`: Origin of coordinate system (e.g., "SSB")

**Returns:**
- New SpaceRock instance from spherical coordinates

*Example:*
```python
rock = SpaceRock.from_spherical("object", phi=0.0, theta=0.0, r=1.0,
                               vr=0.0, vo=1.0, psi=0.0, epoch=epoch)
```

**`from_kepler()`**
```python
@classmethod
def from_kepler(cls, name: str, q: float, e: float, inc: float, arg: float, 
                node: float, true_anomaly: float, epoch: Time, 
                reference_plane: str = "ECLIPJ2000", origin: str = "SSB") -> SpaceRock
```
**Arguments:**
- `name`: Object identifier
- `q`: Perihelion distance in AU
- `e`: Eccentricity
- `inc`: Inclination in radians
- `arg`: Argument of perihelion in radians
- `node`: Longitude of ascending node in radians
- `true_anomaly`: True anomaly in radians
- `epoch`: Time object representing the state epoch
- `reference_plane`: Reference frame (e.g., "ECLIPJ2000")
- `origin`: Origin of coordinate system (e.g., "SSB")

**Returns:**
- New SpaceRock instance from orbital elements

*Example:*
```python
rock = SpaceRock.from_kepler("asteroid", q=2.5, e=0.5, inc=0.4, arg=2.1, 
                            node=1.2, true_anomaly=0.0, epoch=epoch)
```

**`random()`**
```python
@classmethod
def random(cls, epoch: Time, reference_plane: str = "ECLIPJ2000", 
           origin: str = "SSB") -> SpaceRock
```
**Arguments:**
- `epoch`: Time object representing the state epoch
- `reference_plane`: Reference frame (e.g., "ECLIPJ2000")
- `origin`: Origin of coordinate system (e.g., "SSB")

**Returns:**
- New SpaceRock instance with random orbital elements

*Example:*
```python
rock = SpaceRock.random(epoch=Time.now())
```

<h2 style="border-bottom: 3px solid white;">Operation Methods</h2>

**`analytic_propagate()`**
```python
def analytic_propagate(self, epoch: Time) -> SpaceRock
```
**Arguments:**
- `epoch`: The Time object representing the epoch that we want to propagate the SpaceRock to.

Calculates and converts SpaceRock elements at new epoch, in place.

*Example:*
```python
rock = SpaceRock.from_horizons("Arrokoth", Time.now())

ten_years = Time.now() + (365.25 * 10)
rock.analytic_propagate(ten_years)
```

**`analytic_at()`**
```python
def analytic_at(self, epoch: Time) -> SpaceRock
```
**Arguments:**
- `epoch`: The Time object representing the epoch that we want to propagate the SpaceRock to in time 

The same as `analytic_propogate()`, except returning a new SpaceRock object

*Example:*
```python
rock = SpaceRock.from_horizons("Arrokoth", Time.now())

ten_years = Time.now() + (365.25 * 10)
rock_in_ten_years = rock.analytic_at(ten_years)
```

**`observe()`**
```python
def observe(self, observer: Observer) -> Observation
```
**Arguments:**
- `observer`: Observer object with position and velocity state

**Returns:**
- Observation object containing calculated ephemeris

*Example:*
```python
observation = rock.observe(observer)
```

**`change_reference_plane()`**
```python
def change_reference_plane(self, reference_plane: str) -> None
```
**Arguments:**
- `reference_plane`: String specifying the new reference frame (e.g., "ECLIPJ2000")

*Example:*
```python
rock.change_reference_plane("ECLIPJ2000")
```

<h2 style="border-bottom: 3px solid white;">Derived Properties</h2>


These methods calculate various properties from the SpaceRock's state. Usage: `rock.property()`

| Property | Returns | Description |
| --- | --- | --- |
| `a` | `float` | Calculate the semi-major axis in AU |
| `e` | `float` | Calculate the eccentricity |
| `inc` | `float` | Calculate the inclination in radians |
| `node` | `float` | Calculate the longitude of ascending node in radians |
| `arg` | `float` | Calculate the argument of periapsis in radians |
| `true_anomaly` | `float` | Calculate the true anomaly in radians |
| `mean_anomaly` | `float` | Calculate the mean anomaly in radians |
| `conic_anomaly` | `float` | Calculate the conic anomaly in radians |
| `q` | `float` | Calculate the perihelion distance in AU |
| `p` | `float` | Calculate the semi-latus rectum in AU |


<h2 style="border-bottom: 3px solid white;">Setter Methods</h2>


These methods modify SpaceRock properties. Usage: `rock.set_property(value)`

| Method | Returns | Description |
| --- | --- | --- |
| `set_absolute_magnitude` | `None` | Set absolute magnitude (H) |
| `set_gslope` | `None` | Set G-slope parameter |
| `set_mass` | `None` | Set mass in Solar Mass (M☉) |
| `set_x` | `None` | Set x-coordinate in AU |
| `set_y` | `None` | Set y-coordinate in AU |
| `set_z` | `None` | Set z-coordinate in AU |
| `set_vx` | `None` | Set x velocity in AU/day |
| `set_vy` | `None` | Set y velocity in AU/day |
| `set_vz` | `None` | Set z velocity in AU/day |

<h2 style="border-bottom: 3px solid white;">Getter Methods</h2>


These methods access SpaceRock properties. Usage: `rock.property`

| Property | Returns | Description |
| --- | --- | --- |
| `absolute_magnitude` | `float` or `None` | Get absolute magnitude (H) |
| `gslope` | `float` or `None` | Get G-slope parameter |
| `mass` | `float` or `None` | Get mass in Solar Mass (M☉)|
| `x` | `float` | Get x-coordinate in AU |
| `y` | `float` | Get y-coordinate in AU |
| `z` | `float` | Get z-coordinate in AU |
| `vx` | `float` | Get x velocity in AU/day |
| `vy` | `float` | Get y velocity in AU/day |
| `vz` | `float` | Get z velocity in AU/day |
| `evec` | `tuple` | Calculate the eccentricity vector |
| `r` | `float` | Calculate the distance from origin in AU |
| `mu` | `float` | Get the gravitational parameter of origin |

<h2 style="border-bottom: 3px solid white;">Examples</h2>


### Basic Usage
```python
from spacerocks import SpaceRock
from spacerocks.time import Time

# Create SpaceRock
epoch = Time.now()
arrokoth = SpaceRock.from_horizons("Arrokoth", epoch)

# Calculate orbital elements
print(f"Semi-major axis: {arrokoth.a()} AU")
print(f"Eccentricity: {arrokoth.e()}")

# Set physical properties
arrokoth.set_mass(1e-10)
arrokoth.set_absolute_magnitude(15.0)
```

### Observation Calculation
```python
from spacerocks.observing import Observatory

# Set up objects
epoch = Time.now()

telescope = Observatory.from_obscode("w84")
observer = telescope.at(epoch)  

asteroid = SpaceRock.from_horizons("Ceres", epoch, "J2000")

# Calculate ephemeris
observation = asteroid.observe(observer)
print(f"RA: {observation.ra} rad")
print(f"Dec: {observation.dec} rad")
```

<h2 style="border-bottom: 3px solid white;">Notes</h2>

- All positions are in Astronomical Units (AU)
- All velocities are in AU/day
- All angles are in radians
- Masses are in Solar Mass M☉
- SPICE kernels must be loaded before using SPICE-dependent features
- Physical properties are optional and default to 0.0
- The G slope parameter defaults to 0.15 if not specified