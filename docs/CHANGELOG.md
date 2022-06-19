### 2.1.20
- Add option to Simulation to dump data to file, and not keep in memory (not exclusive)

### 2.1.1
- Consolidated `spice` logic into a `SpiceKernel` class
- Created a `PerturberModel` class for managing the perturbers in simulations
- Changed the `SpaceRock.propagate` method to use the new `Simulation` class
- Added functionality to provide a pre-computed `Observer` to `SpaceRock.observe`

### 2.1.0
- Reintroduced `PyOrbfit`, with the class now titled `Bernstein`
- Fixed a memory leak caused by Numpy not taking ownership of C pointers, and thus being missed by garbage collection.

### 2.0.1
- Added rotation curve hinting to the Units class. This speeds up the case of constant H.
- Parallelized the C++ backend with OpenMP
- Added a `rich` progress bar to the `Simulation` `propagate` method.
- Added more involved tests.
- Added `from_horizons` instantiator for the SpaceRock class.
- Changed a bug in the Observer class that was sometimes comparing the length of a string rather than an array.
- Fixed a bug in [varpi, node, arg] relations for retrograde orbits.
- Large speedup to Observer class internal list comprehension by turning a `Quantity` into a `float`
- Renamed `SpaceRock.writeto` to `SpaceRock.write_to`.
- Added a compression kwarg to `SpaceRock.write_to`.

### 2.0.0
- Dropped support for pyOrbfit, as installing across platforms was difficult.
- Began maintaining documentation.
- Added custom subclass of `Rebound` simulation objects.
- Added support for changing coordinate frames between `J2000` and `eclipJ2000` using `change_frame()` method.
- Began including asteroid ephemerides in the distribution.