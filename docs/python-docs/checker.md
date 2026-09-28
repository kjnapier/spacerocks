<h1 style="border-bottom: 5px solid white;">Checker Module</h1>

### Table of Contents
1. [Overview](#overview)
2. [Catalog](#catalog)
3. [check()](#check)
4. [Examples](#examples)
5. [Validation](#validation)
6. [Notes](#notes)

<h2 style="border-bottom: 3px solid white;">Overview</h2>

`spacerocks.checker` identifies detections with known objects, like the MPC's MPChecker. For
each detection it predicts every orbit in a catalog at the detection's time, as seen from its
observer. It then asks whether the detection and the prediction agree, taking both
uncertainties into account.

A `Catalog` can hold three kinds of orbit:

- **The MPC's catalog of every minor planet (MPCORB).** It is downloaded to
  `~/.spacerocks/mpc` and holds about 1.4 million orbits. It has no covariances. Each orbit
  instead carries the MPC's uncertainty parameter U, which the check turns into an along-track
  uncertainty.
- **Orbits you fit with `spacerocks.orbfit`.** These carry their full covariance, including any
  fitted non-gravitational parameters.
- **The MPC's `mpc_orb` records** (from the `get-orb` API). These carry the covariance of the
  Cartesian state.

The check runs in two stages:

1. **Coarse, over the whole catalog.** Each object moves by two-body motion from a reference
   state no more than `max_age` days (30 by default) from the detections. The reference is the
   orbit's own epoch, the catalog's snapshot, or a state integrated for the purpose (N-body, the
   whole catalog at once, then kept only for that call). The object becomes a candidate for a
   detection if its position falls inside a gate. The gate is `radius`, or the region
   `nsigma` sigma from the prediction, whichever is larger, plus a margin for the error of
   two-body motion. That margin is 30″ plus 1″/day² × dt² for perihelia inside 1.3 AU, and
   0.02″/day² × dt² beyond. On real MPCORB orbits, two-body motion from an N-body state drifts
   at most 5″, 190″ and 580″ in 5, 20 and 30 days for perihelia inside 1.3 AU, and 0.2″, 0.7″
   and 1.7″ beyond (`examples/checker_twobody_error.rs`).
2. **Refined, for each candidate.** A full N-body prediction is made: ASSIST's force model with
   light time, as in `orbfit.predict`. The object's uncertainty then comes from one of two
   places:
   - its covariance, mapped through the variational equations;
   - for MPCORB orbits, U's in-orbit longitude uncertainty, applied along the object's motion,
     with 10% of it across.

   A `floor` (0.3″ by default) is added in quadrature. With C_pred the prediction's covariance
   and C_det the detection's, the pair's Mahalanobis distance is
   d = sqrt(rᵀ (C_pred + C_det)⁻¹ r), where r is their offset on the tangent plane. The pair is
   **consistent** when d ≤ `nsigma` (3 by default).

Angles going in are radians, as in `orbfit`. Output offsets and uncertainties are in arcseconds.

<h2 style="border-bottom: 3px solid white;">Catalog</h2>

```python
Catalog.mpcorb(kernel, path=None, download=False, update=False) -> Catalog
Catalog.from_fits(fits, names=None, h=None) -> Catalog   # {name: OrbitFit} or a list of OrbitFits
Catalog.from_mpc_orb(orbs, kernel) -> Catalog            # mpc_orb dict(s), JSON text, or a path
Catalog.load(path) -> Catalog
cat.save(path)
cat.snapshot(epoch, kernel, timescale="tdb", chunk_size=64, method="nbody")
cat.select(indices_or_mask) -> Catalog
cat.extend(other)
cat.index(name) -> int
cat.rock(i) -> SpaceRock
len(cat); cat.names; cat.epoch; cat.states; cat.h; cat.u; cat.has_covariance; cat.snapshot_epoch
```

**`Catalog.mpcorb`** reads `mpcorb_extended.json`, `MPCORB.DAT`, or either gzipped. By default
it reads `~/.spacerocks/mpc/mpcorb_extended.json.gz` (or `$SPACEROCKS_MPC_DIR`). With
`download=True` it fetches the file from the MPC first if it is missing, and `update=True`
downloads a fresh copy.

The first read takes about 40 s: parsing, then integrating the ~1% of orbits whose epoch is
older than the rest to the common epoch. The result is cached next to the file
(`<file>.srcat`, about 250 MB) and reloads in under a second.

Numbered objects are named by their number (`"3666"`), and the others by their principal
designation.

**`snapshot`** integrates every orbit to `epoch` (N-body) and keeps the states. The check uses
two-body motion for at most `max_age` days from a reference state. Detections further than that
from both the orbits' epochs and the snapshot therefore make it integrate the whole catalog on
every call. That takes about 4 minutes for MPCORB 16 months out on 2 cores. A snapshot near your
data, saved once, turns this into seconds:

```python
cat = checker.Catalog.mpcorb(kernel, download=True)
cat.snapshot(2461300.5, kernel)            # e.g. the middle of an observing run
cat.save("mpcorb_2461300.srcat")
# later:
cat = checker.Catalog.load("mpcorb_2461300.srcat")
```

**`from_fits`** takes `OrbitFit`s from `orbfit.fit`/`fit_many`. These keep their covariance
(state plus fitted A1–A3) and their non-gravitational parameters. Fits without an orbit are
left out.

**`from_mpc_orb`** takes the MPC's orbit records. It uses the heliocentric ecliptic Cartesian
state (`CAR`) and its covariance, rotated to J2000.

Catalogs of different kinds can be merged with `extend`.

<h2 style="border-bottom: 3px solid white;">check()</h2>

```python
def check(catalog, ra, dec, epoch, observer, kernel,
          sigma_ra=None, sigma_dec=None, correlation=None, timescale="utc",
          nsigma=3.0, radius=60.0, max_age=30.0, max_uncertainty=600.0,
          floor=0.3, cross_track=0.1, epsilon=1e-9) -> dict
```

**Inputs**

- `ra`, `dec`: radians (ICRF).
- `epoch`: Julian dates in `timescale`, or Time objects.
- `observer`: as for `orbfit.fit`. That is an Observatory, an MPC code, a list of codes, a list of
  Observers, or an (n, 3|6) array of barycentric J2000 states. With velocities, rates are
  predicted too.
- `sigma_ra` (along RA·cos Dec), `sigma_dec`: 1σ in radians, 1″ by default. `correlation` is
  their correlation.

**Options**

- `nsigma`: the consistency threshold.
- `radius` (″): objects predicted this close are reported even when inconsistent, as
  MPChecker's search radius does. With 0, only consistent pairs are reported.
- `max_age` (days): the longest two-body step in the coarse search.
- `max_uncertainty` (″): objects whose 1σ uncertainty is larger than this can't usefully be
  matched. They are only reported within `radius`.
- `floor` (″): added in quadrature to every prediction's uncertainty. It covers what neither
  covariance includes, such as element rounding and force-model differences.
- `cross_track`: for orbits without a covariance, the uncertainty across the motion, as a
  fraction of the along-track one.

**Returns** a dict with one value per reported (detection, object) pair, so
`pandas.DataFrame(result)` gives a table. Pairs are sorted by detection, then with consistent
pairs first, then by decreasing `log_likelihood`.

| key | meaning |
|---|---|
| `detection`, `object`, `name` | index of the detection; index and name of the object in the catalog |
| `ra`, `dec` | predicted astrometric position (radians) |
| `dra`, `ddec` | observed − predicted (″, along RA·cos Dec and Dec, on the tangent plane) |
| `separation` | angle between them (″) |
| `distance` | Mahalanobis distance, with both uncertainties |
| `consistent` | `distance <= nsigma` |
| `log_likelihood` | log Gaussian density of the offset (per ″²). It ranks a precise orbit that fits above a vague one that merely does not exclude the detection |
| `sigma_major`, `sigma_minor`, `pa` | the prediction's 1σ ellipse (″, ″, degrees N through E), `floor` included |
| `cov_ra`, `cov_ra_dec`, `cov_dec` | the prediction's covariance (″²) |
| `ra_rate`, `dec_rate` | predicted dRA/dt and dDec/dt (radians/day; NaN without observer velocities) |
| `delta`, `r_helio` | distance from the observer and from the Sun (AU) |
| `mag` | predicted V (H, G) |
| `from_covariance` | whether the uncertainty came from a covariance (True) or from U (False) |

Objects whose prediction fails (for example, outside the ephemeris' span) are skipped with a
warning.

<h2 style="border-bottom: 3px solid white;">Examples</h2>

```python
import pandas as pd
from spacerocks import SpiceKernel, checker

kernel = SpiceKernel.defaults()
cat = checker.Catalog.mpcorb(kernel, download=True)

# ra, dec in radians; epoch UTC JD; station codes; per-detection sigma in radians
r = pd.DataFrame(checker.check(cat, ra, dec, epoch, stations, kernel, sigma_ra=sig, sigma_dec=sig))
best = r[r.consistent].groupby("detection").head(1)   # the most likely identification
```

Orbits you fitted yourself, together with MPCORB:

```python
fits = orbfit.fit_many(...)            # or {name: orbfit.fit(...)}
mine = checker.Catalog.from_fits({"cand1": fit1, "cand2": fit2})
mine.extend(cat)
r = checker.check(mine, ra, dec, epoch, "W84", kernel, sigma_ra=0.1 / 206265)
```

<h2 style="border-bottom: 3px solid white;">Validation</h2>

The scripts are in `validation/checker/`. They use the real MPC astrometry cached by the
orbit-fitting notebook and the MPCORB of February 2025, whose orbit epoch is 2025 May 5.

**Real detections against MPCORB** (`validate_mpcorb.py`). Up to 400 ground-based detections
since 2000 per object were drawn and checked against the five orbits among 30,000 random
MPCORB orbits. Detection σ is the reported rms, or Vereš et al. (2017) where none is reported.

| object | found | consistent (3σ) | right object ranked first | median d | 95% d |
|---|---|---|---|---|---|
| (3666) Holman | 400/400 | 400 | 400 | 0.34 | 1.02 |
| (433) Eros | 400/400 | 400 | 400 | 0.30 | 0.88 |
| (90377) Sedna | 400/400 | 400 | 400 | 0.00 | 0.01 |
| (99942) Apophis | 400/400 | 399 | 400 | 0.51 | 0.79 |
| (101955) Bennu | 400/400 | 396 | 400 | 0.71 | 1.19 |

Only one pair with another object was consistent: an orbit with a 515″ ellipse. The misses are Bennu detections
with reported rms of 0.2″ that lie 1–2.5″ off, and one Apophis detection 4.5″ off. For a
Gaussian at 3σ, 1.1% misses are expected; 0.25% were seen.

Sedna's MPC U is 5. That is because U is built from the uncertainty of the orbital period, and
Sedna's period is about 11,000 years. So its predicted uncertainty is far larger than its
actual error: the prediction lies 0.27″ from the detections.

**Orbits with a covariance** (`validate_covariance.py`, Holman among 30,000 MPCORB orbits):

| orbit | detections | consistent | median d | 95% d | median σ_major |
|---|---|---|---|---|---|
| fitted to 2015–23 → 2024–26 | 814 | 812 | 0.45 | 1.19 | 0.30″ |
| fitted to 2015–23 → 1990–2014 | 300 | 297 | 0.60 | 1.78 | 0.30″ |
| fitted to 30 days of 2023 → later | 300 | 105 (all 300 with `max_uncertainty` lifted) | 0.59 | 0.90 | 31″ |
| MPC `mpc_orb` → 2024–26 | 814 | 812 | 0.44 | 1.18 | 0.30″ |
| MPC `mpc_orb` → 1990–2014 | 300 | 297 | 0.60 | 1.77 | 0.30″ |
| fit, detections moved 10″ | 814 | 0 | 21 | 31 | 0.30″ |

**The full catalog** (`validate_full.py`, 1,439,233 orbits, 2 cores). The test set is 67 detections of Holman
and Eros from September 2026, 16 months after the catalog's epoch.

| step | time |
|---|---|
| without a snapshot (the check integrates the catalog) | 189 s |
| `snapshot` at the detections' median epoch | 197 s |
| `save` / `load` (254 MB) | 0.2 s / 0.6 s |
| check with the snapshot | 3.5 s |

In both runs the right object was ranked first and consistent for all 67 detections. Fifteen
pairs with other objects were also consistent. Those objects have U 3–5 and 1σ uncertainties
of ~4′, so their ellipses cover the detections without pointing to them. `log_likelihood`
ranks them below the right object.

With `max_uncertainty` lifted to 10°, the 30-day fit above is consistent with all 300 later
detections, and it is still ranked first for all of them. But 1,031 pairs with loosely
determined MPCORB orbits also become consistent, which is why the default is 600″.

<h2 style="border-bottom: 3px solid white;">Notes</h2>

- **U is coarse.** The MPC derives U from the uncertainties of the perihelion time and the
  period. The check uses the upper end of U's runoff bin, exp(1.49 U)″ per decade, as the 1σ
  in-orbit longitude uncertainty. That uncertainty grows linearly from the orbit's epoch (or
  from the last observation, if that is longer ago), plus a tenth of the runoff (a year's
  growth). It is applied as a shift of the object along its orbit (r·σ), plus `cross_track` of
  that across.
  - Growing from the epoch also covers what the elements themselves lose over long
    propagations: rounding to 7 digits, and forces such as Yarkovsky that the MPC's fit
    included and this integration does not. For example, Bennu's 2005 detections fall ~28″
    from its 2025 MPCORB orbit.
  - For precise work, fit the orbits (`from_fits`) or use the MPC's covariances
    (`from_mpc_orb`).
- **Likelihood versus distance.** An orbit with a 1° ellipse is consistent (d < 3) with any
  detection along its line of variations. `max_uncertainty` keeps such orbits out of the
  consistency test. `log_likelihood` ranks the remaining matches.
- **Linear covariance.** The covariance is mapped linearly. For very short arcs (a week, say)
  the linearization is meaningless, and such orbits end up above `max_uncertainty`.
- **Comets.** `mpcorb_extended` holds asteroids only. Comets can be added with `from_fits` or
  `from_mpc_orb`, which keep their non-gravitational parameters.
- **Downloads.** The MPC was not reachable from the development sandbox, so the download itself
  (`download=True`) is untested. Parsing was tested on a real `mpcorb_extended.json.gz` and on
  `MPCORB.DAT` lines built to the MPC's column layout.
