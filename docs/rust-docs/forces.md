# Force models (Rust)

`spacerocks::assist::forces` ports ASSIST's force models, each with its variational terms:

| Type | ASSIST flag | Constructor |
|---|---|---|
| `NewtonianGravity` | `SUN`, `PLANETS`, `ASTEROIDS` | unit struct |
| `EarthHarmonics` | `EARTH_HARMONICS` | `EarthHarmonics::new(&constants)` |
| `SolarJ2` | `SUN_HARMONICS` | `SolarJ2::new(&constants)` |
| `GrEih` | `GR_EIH` | `GrEih::new(&constants).with_sources(n)` |
| `GrSimple` | `GR_SIMPLE` | `GrSimple::new(&constants)` |
| `GrPotential` | `GR_POTENTIAL` | `GrPotential::new(&constants)` |
| `NonGravitational` | `NON_GRAVITATIONAL` | `NonGravitational::default()` / `::comet()` |

`assist_default_forces(&EphemerisConstants::from_kernel(&kernel))` is ASSIST's default set, used by
`SpiceSimulation::horizons`. `EphemerisConstants` reads AU, CLIGHT, J2E, J3E, J4E, RE, J2SUN and
ASUN from the planetary SPK's comment area (`SpiceKernel::integration_constant`), falling back
to DE440.

## Jacobians and variational particles

`Force::apply_acceleration_and_stm` adds each particle's acceleration partials to
`SimulationParticle::stm` (rows 3..6: d a/d r in columns 0..3, d a/d v in columns 3..6) and
`SimulationParticle::nongrav_partials` (d a/d(A1, A2, A3)). IAS15 computes them only when the
simulation has variational particles. Variational accelerations are then
`J_r dr + J_v dv + J_A dA` (`update_variational_accelerations`), so velocity-dependent forces
(relativity, non-gravitational) and non-gravitational parameter partials are handled.

```rust
let mut sim = SpiceSimulation::horizons(&epoch, &kernel)?;
rock.set_nongrav(0.0, -2.9e-14, 0.0);
sim.add(rock)?;
sim.add_full_variation("apophis")?;                  // x, y, z, vx, vy, vz
for a in ["A1", "A2", "A3"] { sim.add_variation(a, "apophis")?; }
sim.integrate(&t1, &kernel)?;
let v = &sim.state.variational_particles;             // d state(t1) / d (state0, A1, A2, A3)
```

See `docs/python-docs/spice_simulation.md` for validation against ASSIST.
