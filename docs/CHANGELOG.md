### 2.0.1
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