<h1 style="border-bottom: 5px solid white;">Simulation Module</h1>

### Table of Contents
1. [Overview](#overview)
2. [Primary Classes](#primary-classes)
3. [Methods](#methods)
4. [Examples](#examples)
5. [Notes](#notes)

<h2 style="border-bottom: 3px solid white;">Overview</h2>


The Simulation module provides functionality for running n-body simulations in spacerocks. It handles orbital dynamics, force interactions, and numerical integration of celestial bodies.

<h2 style="border-bottom: 3px solid white;">Primary Classes</h2>


### Simulation
```python
class Simulation:
    """
    Core class for n-body simulations in spacerocks.
    
    Example:
        sim = Simulation()
        sim.set_epoch(Time.now())
        sim.set_reference_plane("ECLIPJ2000")
        sim.set_origin("ssb")
        sim.set_integrator(Integrator.ias15(timestep=20.0))
    """
```

The Simulation class manages collections of celestial bodies, their interactions, and their evolution over time through numerical integration.

### Force
```python
class Force:
    """
    Base class for gravitational and non-gravitational forces in the simulation.
    
    Available forces:
        - NewtonianGravity: Standard gravitational force between bodies
        - SolarGR: Post-Newtonian corrections for general relativity around the Sun
        - SolarJ2: Solar oblateness perturbation
    
    Example:
        sim.add_force(Force.newtonian_gravity())
        sim.add_force(Force.solar_gr())
        sim.add_force(Force.solar_j2())
    """
```

Forces calculate acceleration on bodies in the simulation. 

### Integrator
```python
class Integrator:
    """
    Numerical integrators for advancing the simulation state.
    
    Available integrators:
        - IAS15: High-accuracy adaptive integrator (default)
        - Leapfrog: Fast fixed-timestep integrator
    
    Example:
        sim.set_integrator(Integrator.ias15(timestep=1.0))
        sim.set_integrator(Integrator.leapfrog(timestep=20.0))
    """
```

Integrators advance the simulation state through numerical integration of equations of motion.

<h2 style="border-bottom: 3px solid white;">Methods</h2>

### Constructor Methods
---
**`new()`**

**Returns:**
- New Simulation instance

*Example:*
```python
sim = Simulation()
```

**`giants()`**
```python
@classmethod
def giants(cls, epoch: Time, reference_plane: str, origin: str) -> Simulation
```
**Arguments:**
- `epoch`: Time object representing the simulation's start time
- `reference_plane`: String specifying the reference frame (e.g., "ECLIPJ2000")
- `origin`: String specifying the origin point (e.g., "SSB")

**Returns:**
- New Simulation object with Sun and giant planets

*Example:*
```python
sim = Simulation.giants(epoch=Time.now(), reference_plane="ECLIPJ2000", origin="SSB")
```

**`planets()`**
```python
@classmethod
def planets(cls, epoch: Time, reference_plane: str, origin: str) -> Simulation
```
**Arguments:**
- `epoch`: Time object representing the simulation's start time
- `reference_plane`: String specifying the reference frame
- `origin`: String specifying the origin point

**Returns:**
- New Simulation object with Sun and all planets

*Example:*
```python
sim = Simulation.planets(epoch=Time.now(), reference_plane="ECLIPJ2000", origin="SSB")
```

**`horizons()`**
```python
@classmethod
def horizons(cls, epoch: Time, reference_plane: str, origin: str) -> Simulation
```
**Arguments:**
- `epoch`: Time object representing the simulation's start time
- `reference_plane`: String specifying the reference frame
- `origin`: String specifying the origin point

**Returns:**
- New Simulation object with full JPL Horizons perturber set

**Example:**
```python
sim = Simulation.horizons(epoch=Time.now(), reference_plane="ECLIPJ2000", origin="SSB")
```

### Operation Methods
---

### Particle Management

**`add()`**
```python
def add(self, rock: SpaceRock) -> None
```
**Arguments:**
- `rock`: SpaceRock object to add to simulation

*Example:*
```python
sun = SpaceRock.from_spice("sun", epoch=epoch, reference_plane="ECLIPJ2000", origin="SSB")
sim.add(sun)
```

**`remove()`**
```python
def remove(self, name: str) -> None
```
**Arguments:**
- `name`: Name of SpaceRock to remove from simulation

*Example:*
```python
sim.remove("custom_body")
```

### Integration Methods
---

**`integrate()`**
```python
def integrate(self, epoch: Time) -> None
```
**Arguments:**
- `epoch`: Target time to integrate to

*Example:*
```python
sim.integrate(epoch + 365.25)  # Integrate one year forward
```

**`step()`**
```python
def step(self) -> None
```
Advances simulation by one timestep.

*Example:*
```python
sim.step()  # Take single integration step
```

### Frame Transformations
---

**`move_to_center_of_mass()`**
```python
def move_to_center_of_mass(self) -> None
```
Centers the simulation on its center of mass.

*Example:*
```python
sim.move_to_center_of_mass()
```

**`change_origin()`**
```python
def change_origin(self, origin: str) -> None
```
**Arguments:**
- `origin`: Name of body to use as new origin

*Example:*
```python
sim.change_origin("sun")
```

### Configuration Methods
---

**`set_epoch()`**
```python
def set_epoch(self, epoch: Time) -> None
```
**Arguments:**
- `epoch`: New simulation epoch

*Example:*
```python
sim.set_epoch(Time.now())
```

**`set_reference_plane()`**
```python
def set_reference_plane(self, reference_plane: str) -> None
```
**Arguments:**
- `reference_plane`: New reference plane identifier

*Example:*
```python
sim.set_reference_plane("ECLIPJ2000")
```

**`set_origin()`**
```python
def set_origin(self, origin: str) -> None
```
**Arguments:**
- `origin`: New origin identifier

*Example:*
```python
sim.set_origin("SSB")
```

**`set_integrator()`**
```python
def set_integrator(self, integrator: Integrator) -> None
```
**Arguments:**
- `integrator`: Integrator object specifying integration method

*Example:*
```python
sim.set_integrator(Integrator.ias15(timestep=1.0))
```

**`add_force`**
```python
def add_force(self, force: Force) -> None
```
**Arguments:**
- `force`: Force object to add to simulation

*Example:*
```python
sim.add_force(Force.solar_gr())
```

### Analysis Methods
---

**`energy()`**
```python
def energy(self) -> float
```
**Returns:**
- Total energy of the system

* Example:*
```python
total_energy = sim.energy()
```

**`get_particle()`**
```python
def get_particle(self, name: str) -> SpaceRock
```
**Arguments:**
- `name`: Name of particle to retrieve

**Returns:**
- SpaceRock object

*Example:*
```python
jupiter = sim.get_particle("jupiter barycenter")
```

<h2 style="border-bottom: 3px solid white;">Examples</h2>

### Basic Integration
```python
from spacerocks.nbody import Simulation, Integrator, Force
from spacerocks.time import Time
from spacerocks.spice import SpiceKernel

# Load required SPICE kernels
kernel = SpiceKernel()
kernel.load("de440s.bsp")
kernel.load("latest_leapseconds.tls")

# Create simulation
epoch = Time.now()
sim = Simulation.giants(epoch=epoch, reference_plane="ECLIPJ2000", origin="SSB")

# Add forces
sim.add_force(Force.solar_gr())
sim.add_force(Force.solar_j2())

# Configure and run
sim.set_integrator(Integrator.ias15(timestep=1.0))
sim.move_to_center_of_mass()
sim.integrate(epoch + 365.25)

# Analyze results
jupiter = sim.get_particle("jupiter barycenter")
print(f"Jupiter position: {jupiter.position}")
```

### Custom Body Integration
```python
from spacerocks.nbody import Simulation
from spacerocks import SpaceRock
from spacerocks.time import Time

# Initialize
sim = Simulation()
epoch = Time.now()

# Add bodies
sun = SpaceRock.from_spice("sun", epoch=epoch, reference_plane="ECLIPJ2000", origin="SSB")
sim.add(sun)

body = SpaceRock.from_kepler(
    name="custom_body",
    q=2.5,           # AU
    e=0.85,          # dimensionless
    inc=0.4,         # radians
    node=1.2,        # radians
    arg=2.1,         # radians
    true_anomaly=0.0,# radians
    epoch=epoch,
    reference_plane="ECLIPJ2000",
    origin="SSB"
)
body.set_mass(1e-10)  # Solar Mass (M☉)
sim.add(body)

# Run simulation
sim.move_to_center_of_mass()
sim.integrate(epoch + 100.0)
```

## Notes
- All positions are in Astronomical Units (AU)
- Velocities are in AU/day
- Angular quantities are in radians
- For planets, use barycenter names (e.g., "jupiter barycenter") except for Earth ("earth") and Moon ("moon")
- Newtonian gravity is included by default; other forces must be added explicitly
- Available integrators:
  - IAS15: High-accuracy, adaptive timestep (default)
  - Leapfrog: Fixed timestep, fast
- Available forces:
  - Newtonian gravity (default)
  - Solar general relativity
  - Solar J2 (oblateness)
- SPICE kernels must be loaded before using SPICE-dependent features
- Automatic unit conversions:
  - SPICE positions: km → AU
  - SPICE velocities: km/s → AU/day