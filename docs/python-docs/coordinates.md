<h1 style="border-bottom: 5px solid white;">Coordinates Module</h1>

### Table of Contents
1. [Overview](#overview)
2. [Primary Classes](#primary-classes)
3. [Methods](#methods)
4. [Examples](#examples)
5. [Notes](#notes)

<h2 style="border-bottom: 3px solid white;">Overview</h2>

The Coordinates module provides tools for handling different coordinate systems in astrodynamics. It consists of two main components:

- `ReferencePlane`: Defines the orientation of coordinate systems
- `Origin`: Specifies the center point of coordinate systems

To fully specify the position of a celestial body, both a reference plane (orientation) and an origin (center point) are required.

<h2 style="border-bottom: 3px solid white;">Primary Classes</h2>

### ReferencePlane
```python
class ReferencePlane:
    """
    Represents commonly used celestial reference planes.
    Each plane has an associated rotation matrix for coordinate transformations.
    """
```

#### Available Reference Planes
| Name | Description |
|------|-------------|
| `J2000` | Equatorial coordinates based on the J2000 epoch |
| `ECLIPJ2000` | Ecliptic coordinates based on the J2000 epoch (default) |
| `INVARIABLE` | Solar system's invariable plane |
| `GALACTIC` | Galactic coordinates |
| `FK4` | Equatorial coordinates based on B1950 epoch |

### Origin
```python
class Origin:
    """
    Represents the center point of a coordinate system.
    Includes the gravitational parameter (μ) of the central body.
    """
```

#### Built-in Origins
| Name | Description |
|------|-------------|
| `SUN` | Solar center |
| `SSB` | Solar System Barycenter |

An origin is a body id and μ. `Origin(name, mu)` makes any other origin. A name SPICE knows
gets that body's NAIF id (`Origin("earth", mu).id == 399`), and any other name gets an id of
its own. The properties are `name`, `mu` and `id` (10 for the Sun, 0 for the SSB).

<h2 style="border-bottom: 3px solid white;">Methods</h2>

### ReferencePlane Methods
---

**`from_str()`**
```python
@classmethod
def from_str(cls, name: str) -> ReferencePlane
```

**Arguments:**
- `name`: String name of reference plane

**Returns:**
- ReferencePlane instance

*Example:*
```python
plane = ReferencePlane.from_str("ECLIPJ2000")
```

**`get_rotation_matrix()`**
```python
def get_rotation_matrix(self) -> numpy.ndarray
```

**Returns:**
- 3x3 rotation matrix for coordinate transformations

*Example:*
```python
rotation = plane.get_rotation_matrix()
```

### Origin Methods
---

**`sun()`**
```python
@classmethod
def sun(cls) -> Origin
```

**Returns:**
- Origin instance centered on the Sun

*Example:*
```python
sun_origin = Origin.sun()
```

**`ssb()`**
```python
@classmethod
def ssb(cls) -> Origin
```

**Returns:**
- Origin instance centered on the Solar System Barycenter

*Example:*
```python
ssb_origin = Origin.ssb()
```

**`custom()`**
```python
@classmethod
def custom(cls, name: str, mu: float) -> Origin
```

**Arguments:**
- `name`: Name of the origin
- `mu`: Gravitational parameter in AU³/day²

**Returns:**
- Custom Origin instance

*Example:*
```python
custom_origin = Origin.custom("Earth", 8.887e-10)
```

<h2 style="border-bottom: 3px solid white;">Examples</h2>

### Basic Reference Plane Usage
```python
from spacerocks.coordinates import ReferencePlane

# Use predefined reference planes
j2000 = ReferencePlane.J2000
ecliptic = ReferencePlane.ECLIPJ2000  # Default plane
galactic = ReferencePlane.GALACTIC

# Get rotation matrix
rotation = galactic.get_rotation_matrix()
```

### Working with Origins
```python
from spacerocks.coordinates import Origin

# Use built-in origins
sun = Origin.sun()
ssb = Origin.ssb()

# Create custom origin
earth_moon = Origin.custom(
    name="Earth-Moon Barycenter",
    mu=8.997011346712499e-10  # AU³/day²
)
```

### Combined Usage
```python
from spacerocks import SpaceRock
from spacerocks.time import Time

# Create a space rock using specific coordinates
epoch = Time.now()
asteroid = SpaceRock.from_horizons(
    "Ceres",
    epoch=epoch,
    reference_plane="ECLIPJ2000",
    origin="SSB"
)
```

<h2 style="border-bottom: 3px solid white;">Notes</h2>

- The default reference plane is ECLIPJ2000
- Reference plane rotations are implemented as 3x3 matrices
- Gravitational parameters (μ) are in AU³/day²
- Custom origins can be created but require gravitational parameter
- When using planetary origins, use the format "[PLANET] BARYCENTER" (e.g., "MARS BARYCENTER")