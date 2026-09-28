# SpiceSimulation: forces and variational equations

`spacerocks.assist.SpiceSimulation` integrates test particles (asteroids, comets, TNOs) with
IAS15 in the field of the Sun, planets, Moon, Pluto and 16 massive asteroids read from a
`SpiceKernel`. Its force model is a port of [ASSIST](https://github.com/matthewholman/assist)
(Holman et al. 2023), including the variational equations of every force.

## Forces

`SpiceSimulation.horizons(epoch, kernel)` uses ASSIST's default model:

| Force | Constructor | Notes |
|---|---|---|
| Non-gravitational | `Force.nongravitational(alpha=1, r0=1, m=2, n=5.093, k=0)` | Marsden et al. (1973) A1, A2, A3 (radial, transverse, normal) times g(r); defaults are JPL's asteroid 1/r² law. Only acts on rocks with `set_nongrav`. `Force.nongravitational_comet()` uses the water-ice g(r). |
| Earth J2, J3, J4 | `Force.earth_harmonics(kernel)` | Pole fixed at the J2000 pole, as in ASSIST and JPL Horizons |
| Solar J2 | `Force.solar_j2(kernel)` | Pole at RA 286.13°, Dec 63.87° |
| Relativity (EIH) | `Force.gr_eih(kernel, sources=1)` | Einstein–Infeld–Hoffmann (PPN β = γ = 1). With `sources=1` only the Sun is a relativistic source; up to 11 (Sun, planets, Moon, Pluto) |
| Newtonian gravity | `Force.newtonian_gravity()` | Point masses of all 27 perturbers |

Also available: `Force.gr_simple(kernel)` (the Sun's one-body post-Newtonian term) and
`Force.gr_potential(kernel)` (a velocity-independent approximation). Constants (J2E, J3E, J4E,
RE, J2SUN, ASUN, AU, CLIGHT) are read from the planetary SPK's comment area, as ASSIST does,
with DE440 values as fallback.

```python
from spacerocks.assist import SpiceSimulation, Force

sim = SpiceSimulation.horizons(epoch, kernel)          # ASSIST defaults
sim.set_forces([Force.gr_eih(kernel, sources=11),     # custom model, smallest first
                Force.earth_harmonics(kernel),
                Force.newtonian_gravity()])
sim.newtonian_only()                                   # point-mass gravity only
Force.assist_defaults(kernel)                          # the default list
```

The extra forces cost about 1.3x (one particle) to 1.6x (64 particles) per step over Newtonian
gravity alone. `SpaceRock.propagate`, `RockCollection.propagate` and `RockCollection.ephemeris`
use the default model.

## Non-gravitational parameters

```python
apophis.set_nongrav(0.0, -2.9e-14, 0.0)   # A1, A2, A3 in AU/day^2 (JPL SBDB values)
apophis.nongrav                           # (0.0, -2.9e-14, 0.0)
```

The parameters are carried by the rock through `propagate`, `RockCollection.propagate` and
`ephemeris`.

## Variational equations

```python
sim = SpiceSimulation.horizons(epoch, kernel)
sim.add(rock)                                  # rock.name == "apophis"
sim.add_full_variation("apophis", nongrav=True)
sim.integrate(t1, kernel)
v = sim.variational_particles                  # shape (9, 6)
stm = v[:6].T                                  # d(state at t1) / d(state at epoch)
d_dA = v[6:].T                                 # d(state at t1) / d(A1, A2, A3)
```

`sim.add_variation(dimension, parent)` adds a single one (`"x"`, `"y"`, `"z"`, `"vx"`, `"vy"`,
`"vz"`, `"A1"`, `"A2"`, `"A3"`).

## Validation

* Every force's analytic partial derivatives agree with finite differences to ~1e-8
  (`tests/forces_tests.rs`), and integrated variational particles agree with finite differences
  of whole propagations.
* `validation/assist/compare.py` integrates five test orbits for 400 days in ASSIST and
  spacerocks with each force switched on in turn. Each force's effect agrees with ASSIST's to
  metres (the integrators' own noise), and the full default model agrees to under 5 m. The
  exception is a 0.002 AU Earth flyby, which agrees to 18 m. The variational particles, including
  the A1–A3 partials, agree to 1e-9–1e-14 relative.

## Integrator settings

```python
sim.epsilon = 1e-11          # IAS15 tolerance (default 1e-9)
sim.min_timestep = 0.001     # step floor in days (default 1e-8)
sim.adaptive_mode = "global" # step-size rule: "prs23" (default), "global" (ASSIST's), "per_particle", "individual"
sim.summation = "kahan"      # round-off control: "kahan" (default), "full" (REBOUND's), "none"
sim.last_timestep, sim.step_epoch
```

IAS15 follows REBOUND's implementation step for step. With the same step rule it reproduces
ASSIST to millimetres on the (5303)–Ceres encounter, and to within a few percent against JPL
over 10^5 days (`validation/assist/examples.py`).

The two rules fail in different places, in ASSIST as in spacerocks. PRS23, the default,
under-resolves encounters with small perturbers: (5303) Parijskij's pass by Ceres comes out
50–17,000 km off at any epsilon. At epsilon 1e-9 it is also about 100 times less accurate than the
global rule against JPL over 10^4 days (Holman: 140 m vs 1.6 m). The last-term "global" rule that
ASSIST uses resolves those, but does not converge through an Earth flyby: without a step floor it
takes thousands of tiny steps and its answer moves by metres to kilometres with epsilon. PRS23
converges there.

`summation` controls round-off. `"kahan"` carries a Kahan compensation term for every position
and velocity from step to step. Where round-off limits the error (tight epsilon, long
integrations), it cuts the error 4–24 times, for well under 1% more time. `"full"` also
compensates the corrector's sums, as REBOUND and ASSIST do. It costs about 10% more and was not
more accurate in the validation runs.

Simulation time is kept as a reference epoch (`state.jd_ref`) plus days since it, as in ASSIST.
Ephemeris lookups are made at the reference epoch plus the offset, so absolute Julian dates (40 µs
resolution) never limit the integration.
