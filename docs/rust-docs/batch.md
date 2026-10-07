# Batch propagation and ephemerides (Rust)

`spacerocks::batch` is the fast path for many rocks and many epochs. Each function takes a
`Population` (one reference plane and origin for all its bodies), and has a `_rocks` / `_batch`
twin for a slice of `SpaceRock`s that may each have their own. Both forms run the same code and
give identical results.

```rust
use spacerocks::batch::{self, BatchOptions, Method};
use spacerocks::{Observatory, Population, SpaceRock, SpiceKernel, Time};

let kernel = SpiceKernel::defaults()?;
let mut pop = Population::from_rocks(rocks)?;

// Move a population to one epoch (in place; its plane and SUN/SSB origin are kept).
batch::propagate(&mut pop, &Time::new(2461000.5, "tdb", "jd")?, &kernel, &BatchOptions::default())?;

// Ephemerides: one observer per epoch, epochs in any order.
let w84 = Observatory::from_obscode("W84")?;
let observers: Vec<_> = epochs.iter().map(|t| w84.at(t, "J2000", "SSB", &kernel)).collect::<Result<_, _>>()?;
let eph = batch::ephemeris(&pop, &observers, &kernel, &BatchOptions::default())?;
let a = eph.get(rock_index, epoch_index);   // Apparent { ra, dec, ra_rate, dec_rate, range, ... }

// Barycentric J2000 states only, row-major [body][target].
let states = batch::states_at(&pop, &[2461000.5, 2461010.5], &kernel, &BatchOptions::default())?;
```

| Population | Slice of rocks |
|---|---|
| `propagate(&mut pop, ...)` | `propagate_batch(&mut rocks, ...)` |
| `ephemeris(&pop, ...)` | `ephemeris_rocks(&rocks, ...)` |
| `states_at(&pop, ...)` | `states_at_rocks(&rocks, ...)` |

`propagate` on a population with a custom origin is an error for N-body (the rocks form returns
such rocks about the SSB instead).

`BatchOptions`:

| Field | Default | Meaning |
|---|---|---|
| `method` | `Method::NBody` | `NBody` (IAS15 with the `SpiceSimulation::horizons` force model: ASSIST's defaults) or `TwoBody` (Keplerian about each rock's origin) |
| `chunk_size` | 64 | Maximum rocks per simulation (smaller groups are used if needed to keep every thread busy) |
| `parallel` | `true` | Run groups with rayon |
| `with_states` | `false` | `ephemeris` also returns barycentric J2000 states |
| `perturber_cache` | `false` | Share a `PerturberCache` between groups |

## How it works

* Rocks are converted to barycentric J2000 with one Sun lookup per distinct epoch, grouped by
  epoch, sorted by perihelion distance, and split into groups. Each group is one
  `SpiceSimulation`, so the 27 perturber states are evaluated once per step for the whole group.
  Grouping by perihelion distance keeps near-Earth objects (short steps) away from TNOs (long
  steps).
* `ephemeris` sorts and de-duplicates the epochs, integrates each group forward through the later
  epochs and backward through the earlier ones, and reads each state off the IAS15 dense output
  (`SpiceSimulation::integrate_jd`), so no integration is restarted.
* Observables come from `observing::apparent`, the allocation-free core of `SpaceRock::observe`
  (same light-time correction and formulas, bit-identical results).

## Measured (2 cores, de440s + sb441-n16, ASSIST's default force model)

| Workload | Before | Batch | Speedup |
|---|---|---|---|
| 2000 TNOs, +100 d | 74 ms (per-rock `propagate`, parallel) | 11.4 ms | 6.4x |
| 2000 main-belt objects, +100 d | 169 ms | 36 ms | 4.7x |
| 2000 mixed NEO/MBA/TNO, +100 d | 185 ms | 59–65 ms | ~3x |
| 1000 rocks x 20 epochs, W84 ephemeris | 800 ms (propagate + observe per epoch) | 29 ms | 28x |
| same, two-body | | 6.2 ms | 129x |

Batch ephemerides agree with propagate-then-observe to ~0.002 mas. Timings on this 2-core
machine vary by ±15% between runs.

## Accuracy

Rocks in a group share its step size. The shared step is never longer than a member's own, so
each rock is integrated at least as accurately as alone; results differ from
`SpaceRock::propagate` at the integrator's tolerance (a few metres after 200 days for main-belt
objects and TNOs). The largest differences come from planetary encounters that the single-rock
step control under-resolves. In one test, a near-Earth object's Mars encounter gave a 34 km
difference; there the batch result was within 7 m of a tight-tolerance (`epsilon = 1e-12`)
reference and the single-rock result was not. `chunk_size: 1` reproduces `SpaceRock::propagate`
exactly.

## PerturberCache

`assist::PerturberCache::build(kernel, ids, start_jd, end_jd)` fits each perturber's barycentric
state with degree-15 Chebyshev polynomials on intervals chosen per body (2 days for the
asteroids, 4–8 days for Mercury, Venus, Earth and Moon, 32 days for the outer planets). The fits
match the kernel to within ~5e-14 AU, which is the kernel's own jitter from Julian-date epochs.
`SpiceSimulation::set_perturber_cache(Arc<PerturberCache>)` makes IAS15 read perturber states from
the fits for epochs inside the span (it falls back to the kernel outside it).

A lookup of all 27 bodies takes 1.0 µs instead of 2.1 µs. That speeds up simulations with one or
a few particles by 1.2–1.6x. Once many rocks share a simulation, the force calculation dominates,
so the batch functions leave the cache off by default. Building one costs ~0.07 ms per day of span.

## Reproducing

```
SPACEROCKS_KERNELS=/path/to/kernels cargo run --release --example bench_batch
SPACEROCKS_KERNELS=/path/to/kernels cargo test --release --test batch_tests
```
