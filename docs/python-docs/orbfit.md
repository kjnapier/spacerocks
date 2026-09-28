<h1 style="border-bottom: 5px solid white;">OrbitFit Module</h1>

### Table of Contents
1. [Overview](#overview)
2. [Functions](#functions)
3. [OrbitFit](#orbitfit)
4. [Examples](#examples)
5. [Notes](#notes)

<h2 style="border-bottom: 3px solid white;">Overview</h2>

`spacerocks.orbfit` determines orbits from optical astrometry, sky-motion rates (streaks) and radar
delay/Doppler, in any mix. It is a port of the orbit fitter in
[layup](https://github.com/Smithsonian/layup) and reproduces its results (see
`validation/layup/`):

1. The detections are split into arcs at gaps longer than 90 days; the longest arc is fitted first.
2. Gauss's method runs on a triplet from that arc spanning about 15° of mean anomaly for a
   main-belt orbit (or its first, middle and last detections when the arc is shorter). It gives
   up to eight candidate orbits.
3. Candidates that badly miss the detections (80th-percentile residual above 1000σ) are dropped.
   The rest are fitted to the longest arc, and the converged fit with the smallest χ² wins.
   Candidates passing within 0.1 AU of an observer are tried only if nothing else converges.
   If no Gauss candidate converges, or Gauss finds no root at all, a Bernstein–Khushalani linear
   fit to the longest arc (see `bk_iod()`) seeds one more fit. This is layup's `iod="auto"`.
   With `iod="herget"`, Herget's method on the longest arc (see `herget_iod()`) gives the
   candidate instead, and there is no fallback.
4. The winner is fitted to all detections, including rates and radar. If that fails, the arcs are
   added one at a time. Initial orbits come from the optical detections only, so radar needs
   either at least three optical detections or an `initial` orbit.
5. Optionally, non-gravitational parameters (A1, A2, A3) are fitted jointly with the state.

With `robust=True` (not in layup, and off by default), the fit also survives real survey data:

- **Outliers are rejected.** A detection whose normalized residual (the rms over its measurements
  of residual/σ) exceeds `outlier_sigma` (default 4) is left out. The fit is redone from the
  current orbit, and every detection is re-evaluated each round (rejected ones can come back),
  until the set stops changing.
- **Short windows seed the fit when layup's pipeline fails.** Layup fits each Gauss orbit to the
  whole longest arc at once. For a near-Earth object observed for months, that diverges even from
  a good start (Bennu's 313-day 2011–12 arc, for example). Robust mode then tries windows of at
  most 60 days, the most detections first. It runs the pipeline on each window and widens the
  first window that converges step by step to all the detections, rejecting outliers along the
  way.
- **Outliers in the IOD arc don't sink the fit.** Layup discards any initial orbit whose fit has
  χ²/ndof > 10, so one bad detection in the arc used for the initial orbit fails every candidate.
  If the steps above give nothing, they are repeated with that test lifted. The final fit still
  applies `chi2_threshold`.
- **Non-grav fits start from every detection.** The joint fit starts both from all detections
  and from those that survived the gravity-only fit, and keeps whichever settles on more.
  Otherwise, precise radar of a Yarkovsky drifter can look like an outlier to the gravity-only
  orbit.

`OrbitFit.used` marks the detections in the final fit. `OrbitFit.residuals` covers every
detection at the final orbit, the rejected ones included. If no stage settles, the result is
layup's.

Every fit is a Levenberg–Marquardt differential correction with ASSIST's force model (Sun, planets,
Moon, Pluto, 16 massive asteroids, GR, Earth J2–J4, solar J2, non-gravitational forces). The
partial derivatives come from the variational equations, and each detection is corrected for light
time. Rates follow layup's streak model, and radar follows its two-leg round-trip model (delay,
exact round-trip Doppler, Shapiro delay). The transmitter's state is computed exactly at the
transmit time where layup extrapolates it (see Notes).

All inputs and outputs are plain arrays. Angles are in radians (J2000/ICRF), and states are
barycentric J2000 in AU and AU/day.

<h2 style="border-bottom: 3px solid white;">Functions</h2>

**`fit()`**
```python
def fit(ra, dec, epoch, observer, kernel,
        sigma_ra=None, sigma_dec=None, timescale="utc",
        initial=None, prior=None, max_update_sigma=4.0, nongrav=None, nongrav_thresholds=None, gofr=None,
        max_iter=100, epsilon=1e-9, chi2_threshold=10.0, arc_gap=90.0,
        iod="auto", robust=False, outlier_sigma=4.0, conv_frac=0.0, per_arc=False, name="rock",
        ra_rate=None, dec_rate=None, sigma_ra_rate=None, sigma_dec_rate=None,
        delay=None, doppler=None, sigma_delay=None, sigma_doppler=None,
        frequency=None, transmitter=None) -> OrbitFit
```
- `ra`, `dec`: arrays of astrometry (radians). Use NaN for detections without optical astrometry
  (radar rows).
- `epoch`: a list of `Time`s, or an array of Julian dates in `timescale`. For radar this is the
  receive time.
- `observer`: where each detection was made. Any of:
  - an `Observatory`, or an MPC code (`"X05"`), for a single site;
  - a sequence of MPC codes, one per detection;
  - a list of `Observer`s;
  - an `(n, 3)`, `(n, 6)` or `(n, 9)` array of barycentric J2000 positions (AU), then velocities
    (AU/day), then accelerations (AU/day²). Rates and radar need velocities.
- `sigma_ra`, `sigma_dec`: on-sky 1σ uncertainties in radians (σ of RA·cos Dec, like ADES
  `rmsRA`). Each can be a scalar or one value per detection. The default is 1″ (1/206265 rad,
  layup's value).
- `initial`: a `SpaceRock` to start from, skipping steps 1–4.
- `prior`: an `OrbitFit` from an earlier call, to bring up to date with the detections passed
  now, as layup's incremental fitting does (its issue #419). The fit remembers a 64-bit key for
  each of its detections, a hash of every input the fit uses. `OrbitFit.route` then says what
  happened:
  - `"skip"`: the detections are the same (in any order), so the prior comes back unchanged.
  - `"sequential"`: detections were only added. Only the new ones are fitted, at the prior's
    epoch, with the prior's covariance as a Gaussian prior on the parameters it fitted. That
    integrates just the new detections, and equals a refit of all of them as long as the prior
    is close to Gaussian. `chi2` and `ndof` are then over the new detections alone (layup's
    convention), and residuals are NaN for the old ones.
  - `"sequential_fallback"`: that update moved the orbit more than `max_update_sigma` prior
    sigmas, or failed, so all detections were refitted from the prior.
  - `"full"`: detections were removed or changed, so all were refitted from the prior.
  - `"cold"`: the prior did not converge, so the fit starts from scratch.

  A per-arc prior is always refitted.
- `nongrav`: which non-gravitational parameters to fit: `"A2"`, `"A1A2A3"`, `["A1", "A2"]`, or
  `True` (A2). They are fitted after the gravity-only orbit converges. If the joint fit fails, the
  gravity-only orbit is returned.
  `nongrav="auto"` chooses the model as layup's `fit_nongrav="auto"` does. The gravity-only fit is
  kept if its χ²/ndof ≤ 1.5. Otherwise the ladder is A2, A1 and A3 alone, then A1+A2, then all
  three. It adopts the first tier with a warranted model: one that converges, lowers χ² by more
  than 9 per added parameter, and has every added parameter above 3σ. Within a tier, the model
  with the smallest χ² wins. `nongrav_thresholds=(accept_reduced_chi2, delta_chi2_per_param,
  nsigma)` changes those three numbers.

  The χ²/ndof gate assumes calibrated uncertainties. With conservative ones, the gravity-only
  χ²/ndof sits below 1 whether or not an acceleration is there, and nothing is tried. Apophis on
  its MPC-reported uncertainties has χ²/ndof 0.3 and keeps the gravity-only orbit, although its
  A2 is a 150σ detection. Pass `nongrav_thresholds=(0.0, 9.0, 3.0)` to always walk the ladder,
  at the cost of up to five extra fits. With `robust=True`, the model is chosen on the detections
  that survive the gravity-only fit.
- `gofr`: the Marsden g(r) law as `[alpha, m, n, k, r0]`. The default is the inverse-square law
  used for asteroids. For comets, pass the water-ice law `[0.1112620426, 2.15, 5.093, 4.6142, 2.808]`.
- `epsilon`: the IAS15 tolerance.
- `chi2_threshold`: converged fits with χ²/ndof above this get flag 2.
- `arc_gap`: the gap (days) that separates arcs.
- `iod`: `"auto"` (Gauss, then the Bernstein–Khushalani fallback), `"gauss"` (Gauss only), or
  `"herget"` (Herget's method, no fallback).
- `engine`: layup's fitting engine for the pipeline's gravity-only fits (screening, the fit to all
  detections, and the arc-by-arc build-up). `"cartesian"` (the default) fits the barycentric
  state. `"bk_native"` fits Bernstein–Khushalani parameters instead: the direction and inverse
  distance of the position, and the velocity scaled by the inverse distance, in a frame along
  the mean line of sight. It adds a fixed prior on the line-of-sight velocity from the
  bound-orbit condition |v|² < 2μ/r. That steadies short arcs of distant objects, where angles
  barely constrain the radial motion. Its χ² includes the prior's term. As in layup, fits with
  non-gravitational parameters stay Cartesian. It also applies to the fit from an `initial`
  orbit, which layup always fits Cartesian. Unlike layup's, it
  also fits rate and radar rows.
- `robust`: outlier rejection and the short-window fallback (see Overview). `outlier_sigma` is
  the rejection threshold.
- `per_arc`: piecewise-constant non-gravitational parameters, as layup's `per_arc=True` does for
  linking comet apparitions. The detections before the fit epoch (arc A) and after it (arc B) each
  get their own values of the parameters chosen with `nongrav`, with one shared state and g(r).
  Give an `initial` orbit whose epoch lies between the two apparitions. Arc A's values are in
  `nongrav` and `nongrav_sigma`, arc B's in `nongrav_arc2` and `nongrav_arc2_sigma`, and
  `fit_many` returns those arrays too. This needs an explicit `nongrav`, not `"auto"`.
- `conv_frac`: layup's scaled convergence test (its issue #477). By default a fit has converged
  when every component of the step is below 1e-12 (AU, AU/day, AU/day²). That demands more
  digits of a poorly determined parameter than the data support, and it runs into the
  integrator's noise floor. With `conv_frac > 0`, each parameter's bar becomes
  max(1e-12, `conv_frac`·σ), σ being its formal uncertainty, once a step has been accepted.
  0.1 stops within a tenth of a sigma of the minimum, usually in fewer iterations.
- `ra_rate`, `dec_rate`: sky-motion rates in radians/day, NaN where not measured. `ra_rate` is the
  on-sky rate cos(Dec)·dRA/dt, as in ADES `raRate`. It is **not** the `ra_rate` of
  `RockCollection.ephemeris`, which is dRA/dt. ADES rates in ″/hour convert by
  `π/(180·3600)·24`. `sigma_ra_rate`, `sigma_dec_rate` default to 24″/day.
- `delay`: round-trip radar delay in seconds. `doppler`: Doppler shift in Hz at the transmit
  `frequency` (Hz); approaching objects have positive shifts. Use NaN where not measured.
  `sigma_delay` is in seconds (default 1 µs) and `sigma_doppler` in Hz (default 1 Hz). JPL's
  radar astrometry is in µs and Hz, so divide the delays by 10⁶.
- `transmitter`: for bistatic radar, the transmitting antenna, as MPC code(s) or an `(n, 6)` array
  of its barycentric state at the transmit time. By default the receiving station transmitted.

The GIL is released while fitting.

**`fit_many()`**
```python
def fit_many(ids, ra, dec, epoch, observer, kernel, sigma_ra=None, sigma_dec=None,
             timescale="utc", prior=None, max_update_sigma=4.0, nongrav=None, nongrav_thresholds=None,
             gofr=None, max_iter=100, epsilon=1e-9,
             chi2_threshold=10.0, arc_gap=90.0, iod="auto", robust=False, outlier_sigma=4.0, conv_frac=0.0,
             per_arc=False,
             ra_rate=None, ...,
             transmitter=None) -> dict
```
Fits many objects in parallel. `ids` labels each detection with its object; the other arguments
are as for `fit()`. Returns a dict of arrays with one row per object, in order of first appearance:
`id`, `flag`, `chi2`, `ndof`, `niter`, `epoch` (TDB JD), `state` `(m, 6)`, `covariance`
`(m, 6, 6)`, `nongrav` `(m, 3)` and `nongrav_sigma` `(m, 3)`. It also returns `used`, one bool
per detection in input order: true where the detection is in its object's final fit.

For updating later, it also returns:
- `fit_nongrav` `(m, 3)` bool;
- `parameter_covariance`, a list of `(npar, npar)` arrays;
- `fingerprint`, 16 hex digits per object, the same for the same set of detections in any order;
- `keys`, a list of `uint64` arrays: the detections each object's fit covers;
- `route`.

Pass the whole dict back as `prior=` with the current detections. Each object is then routed as
in `fit(..., prior=...)`: layup's `incremental_orbitfit`. Objects without a converged prior go
`"cold"`. `route` is None without a prior.

**`sequential_update()`**
```python
def sequential_update(prior, ra, dec, epoch, observer, kernel, sigma_ra=None, sigma_dec=None,
                      timescale="utc", max_update_sigma=4.0, max_iter=100, epsilon=1e-9,
                      chi2_threshold=10.0, name="rock", ra_rate=None, ..., transmitter=None) -> OrbitFit
```
Updates the converged `prior` with *new* detections only (layup's `run_sequential_update` and
its gate), with no refit. The result keeps the prior's epoch and reports the posterior
covariance, and the χ² and ndof of the new detections. Flag 8 means the orbit moved more than
`max_update_sigma` prior sigmas, so the linearization behind the update can't be trusted: refit,
or use `fit(..., prior=...)`, which does. Flag 7 means the prior covariance is not positive
definite. Its `keys` are the prior's plus the new ones, so it can serve as a later `prior`.

**`residuals()`**
```python
def residuals(rock, ra, dec, epoch, observer, kernel, timescale="utc", epsilon=1e-9,
              ra_rate=None, ..., transmitter=None) -> np.ndarray
```
Residuals (observed minus computed) of a `SpaceRock`, with shape `(n, 6)`. The columns are RA
(times cos Dec) and Dec in radians, the RA and Dec rates in radians/day, radar delay in seconds
and Doppler in Hz, with NaN where a detection has no such measurement. They use the same force
model and light-time correction as the fit, including the rock's non-gravitational parameters.

**`bk_iod()`**
```python
def bk_iod(ra, dec, epoch, observer, kernel, sigma_ra=None, sigma_dec=None,
           timescale="utc", epoch0=None, name="rock") -> SpaceRock | None
```
Bernstein–Khushalani linear initial orbit (layup's `run_bk_iod`). Over a short arc and ignoring
gravity, each detection's gnomonic tangent-plane coordinates are linear in five of the six
Bernstein–Khushalani parameters: the direction (α, β), the inverse distance γ, and the sky-plane
motion. The line-of-sight velocity is fixed at zero. So the orbit comes from one weighted least
squares, with no iteration and no need for three well-spaced detections. That makes it a good IOD
for short arcs of distant objects, where Gauss's method is ill-conditioned or finds no root.

It returns a barycentric J2000 `SpaceRock` at `epoch0` (a TDB Julian date; by default the middle
detection in time), or `None` if there are fewer than three detections or the solution is
unphysical. Use the result as a seed (`fit(..., initial=rock)`), not as a final orbit: its radial
velocity is zero, and it assumes straight-line motion. On layup's own synthetic tests the distance
comes out to better than 0.01% for classical TNOs and within 4–8% for the main belt over a 10-day arc.

**`predict()`**
```python
def predict(fit, epoch, observer, kernel, timescale="utc", epsilon=1e-9) -> dict
```
Where a fitted orbit (an `OrbitFit`) will be seen from `observer` (as for `fit()`) at `epoch`, and
how well that is known (layup's `predict`). The orbit is integrated with its fitted
non-gravitational parameters and the variational equations, corrected for light time, and its
covariance is mapped linearly onto the sky. Returns a dict of arrays:
- `epoch` (TDB JD);
- `ra`, `dec` (radians, astrometric);
- `delta` (AU);
- `covariance` `(n, 2, 2)` in radians² along (RA·cos Dec, Dec);
- the 1σ error ellipse: `sigma_major` and `sigma_minor` (arcsec), and `pa` (degrees, North
  through East).

layup maps only the 6×6 state covariance. This maps the fit's full covariance, so fitted A1–A3
(and per-arc values) add their share.

**`comet_orbits()`**
```python
def comet_orbits(orbit, kernel, reference_distance=250.0, gofr=None, epsilon=1e-9) -> dict
```
The original and future orbits of a long-period comet (layup's `comet`): the barycentric
osculating orbit, with the mass of the Sun and planets, where the comet crosses
`reference_distance` AU inbound (`"original"`) and outbound (`"future"`). `orbit` is an `OrbitFit`
or a `SpaceRock`, and keeps its non-gravitational parameters (pass the water-ice `gofr` for
comets). Each entry is `None` if the orbit never gets that far, or else a dict:
- `epoch`;
- `distance`;
- `reached`;
- `inv_a` (1/AU; ×10⁶ for the CODE catalogue's units);
- `a`, `e`, `q`;
- `inc` (degrees, to the ecliptic).

If the planetary ephemeris ends before the comet gets there, the elements are taken where it
ended (`reached=False`). Out there the planets act as one mass at the barycenter, so this changes
1/a by at most ~1e-7/AU: on 25 CODE comets stopped at 92–168 AU by `de440s`'s 1849 start, the
median change was 0.01×10⁻⁶ and the largest 0.1×10⁻⁶ /AU.

**`herget_iod()`**
```python
def herget_iod(ra, dec, epoch, observer, kernel, sigma_ra=None, sigma_dec=None,
               timescale="utc", name="rock") -> SpaceRock | None
```
Herget's initial orbit (layup's `herget_iod`). The ranges to the first and last detections are
adjusted until the orbit through the two positions fits the detections in between: a two-body
shooting method gives the velocity, ASSIST integrations give how each detection moves with each
range, and the 2×2 normal equations correct both. It starts from ranges of 2 AU, then 5, then 40,
and stops once the ranges move by less than 0.003 AU. Every detection passed is used, so pass
one arc (the pipeline uses the longest). It returns a barycentric J2000 `SpaceRock` at the first
detection, or `None`. Like `bk_iod()`, it is a seed for `fit(..., initial=rock)`.

**`veres_sigma()`**
```python
def veres_sigma(station, epoch, catalog=None, program=None, timescale="utc") -> np.ndarray
```
One-sigma astrometric uncertainties in radians (per axis, on-sky), after Vereš et al. (2017,
Icarus 296, 139). They are assigned exactly as layup's `weight_data=True` does: by station, by
date (703, 691 and 644 improved at set dates), and for a few stations by star catalog and MPC
program code. Unlisted stations get 1″ with a catalog and 1.5″ without.

`station` is an MPC code or one per detection. `epoch` is as for `fit()`. `catalog` (the ADES
`astCat` name, e.g. `"Gaia2"`, or the MPC one-letter code, e.g. `"V"`) and `program` are None, a
string, or one per detection; None, NaN and `""` count as not given. Pass the result as both
`sigma_ra` and `sigma_dec`.

**`debias()`**
```python
def debias(ra, dec, epoch, catalog, table=None, timescale="utc", download=True) -> (np.ndarray, np.ndarray)
```
Star-catalog debiasing (Eggl, Farnocchia, Chamberlin & Chesley 2020), as layup applies it with
`debias_data=True`. It returns the corrected `(ra, dec)` in radians. Each detection's catalog
offset in its HEALPix pixel, plus the catalog's proper-motion error since J2000, is removed.
Apply it before fitting.

- `catalog`: the ADES `astCat` name (`"UCAC4"`) or the MPC one-letter code (`"q"`), as a string or
  one per detection. The table covers 26 catalogs. Anything else leaves the detection unchanged:
  Gaia DR2, DR3 and EDR3, UCAC-5, None, NaN or `""`. Gaia-referenced astrometry needs no
  correction.
- `table`: the path of JPL's `bias.dat` (from `debias_hires2018.tgz`, nside 256). By default it's
  in `$SPACEROCKS_DEBIAS_DIR` or `~/.spacerocks/debias`, and is downloaded from JPL if missing
  and `download` is true. The first load writes a binary copy next to it (`bias.bin`, 327 MB),
  and later loads memory-map that. The table is loaded once per process.

**`observers()`**
```python
def observers(station, epoch, kernel, sys=None, ctr=None, pos=None, vel=None, timescale="utc") -> np.ndarray
```
Barycentric J2000 observer states, shape `(n, 6)` (AU, AU/day), to pass as `fit(observer=...)`.
A detection that comes with its observer's position (ADES `sys`, `ctr`, `pos1`–`pos3`, optionally
`vel1`–`vel3`) is placed as layup places it:

- **`sys="ICRF_KM"` / `"ICRF_AU"`**: a geocentric ICRF vector in km or AU, added to the Earth's
  barycentric state. The optional velocity is in km/s or AU/day; without one, the observer moves
  with the Earth's center. Satellites such as TESS (C57) and WISE (C51) report this way.
- **`sys="WGS84"`**: east longitude and geodetic latitude in degrees and height in metres, for
  roving observers (247, 270). They are placed like fixed stations, so they turn with the Earth.

`ctr` must be 399 (the default). Detections without a position are placed from their MPC station
codes: fixed stations from their parallax constants, and spacecraft from a loaded SPK or from JPL
Horizons.

`pos` and `vel` are `(n, 3)` arrays with NaN rows where not given. `sys` and `ctr` are None, one
value, or one per detection.

**`occultation_radec()`**
```python
def occultation_radec(ra_star, dec_star, delta_ra, delta_dec) -> (np.ndarray, np.ndarray)
```
The object's RA and Dec (radians) from ADES occultation astrometry. The inputs are the occulted
star's `raStar`, `decStar` and the object's offset from it, `deltaRA` (which includes cos Dec) and
`deltaDec`, all in radians. ADES gives the star in degrees and the offsets in arcsec. The offset
is applied on the tangent plane at the star.

Occultation astrometry is precise to a few milliarcseconds. Use the reported `rmsRA`/`rmsDec`
as σ, not Vereš weights or a floor. The observer's position comes with the detection (station
275), so pass `observers(...)`. The correlation `rmsCorr` is not used, since weighting is
diagonal.

**`gauss()`**
```python
def gauss(o1: Observation, o2: Observation, o3: Observation, min_distance: float) -> list[SpaceRock] | None
```
Gauss's method on three `Observation`s. It returns barycentric J2000 candidates at the middle
observation, largest distance first, or `None` if there are no real roots beyond `min_distance`
(AU).

<h2 style="border-bottom: 3px solid white;">OrbitFit</h2>

| Attribute | Description |
|---|---|
| `rock` | The fitted orbit as a barycentric J2000 `SpaceRock`, with A1–A3 when they were fitted |
| `state` | `[x, y, z, vx, vy, vz]`, AU and AU/day |
| `epoch` | Epoch of `state`, a TDB `Time` |
| `covariance` | `(npar, npar)`: the state, then the fitted non-grav parameters |
| `state_covariance` | `(6, 6)` |
| `nongrav`, `nongrav_sigma` | A1, A2, A3 (AU/day²) and their 1σ uncertainties (NaN if not fitted); arc A's with `per_arc` |
| `nongrav_arc2`, `nongrav_arc2_sigma` | With `per_arc`, arc B's A1, A2, A3 and their uncertainties (NaN otherwise) |
| `chi2`, `ndof`, `niter` | χ², degrees of freedom (2n − npar), iterations of the last fit |
| `residuals` | `(n, 6)` residuals at the last iteration, as for `residuals()`. With `robust=True`, all detections at the final orbit. NaN rows if there is no orbit |
| `used` | `(n,)` bool: the detections in the final fit (all, unless `robust=True` rejected some); empty if there is no orbit |
| `route` | With `prior`: `"skip"`, `"sequential"`, `"sequential_fallback"`, `"full"` or `"cold"`; else None |
| `keys`, `fingerprint` | `uint64` keys of the detections the fit covers, and their order-independent fingerprint (16 hex digits) |
| `flag`, `status`, `converged` | Outcome (below) |

| `flag` | Meaning |
|---|---|
| −1 | Not attempted: fewer than three detections, a non-finite value, or a detection before 1801 |
| 0 | Converged |
| 1 | The differential correction did not converge |
| 2 | χ²/ndof above `chi2_threshold` |
| 3 | Candidate orbits were found, but none converged on the longest arc |
| 4 | The longest arc converged, but adding the other arcs failed |
| 5 | No candidate orbits: Gauss found none, and the Bernstein–Khushalani fallback was not possible or not enabled |
| 6 | Converged, but a fitted non-grav parameter has no usable variance |
| 7 | Sequential update: the prior covariance is not positive definite (without a refit) |
| 8 | Sequential update: the orbit moved more than `max_update_sigma` prior sigmas (without a refit) |
| 9 | Converged, but the hyperbolic excess speed exceeds 200 km/s |

The flags are layup's.

<h2 style="border-bottom: 3px solid white;">Examples</h2>

### Real MPC data, with outliers and long NEO arcs
```python
fit = orbfit.fit(ra, dec, epochs, stations, kernel, sigma_ra=sig_ra, sigma_dec=sig_dec, robust=True)
print(fit.used.sum(), "of", len(ra), "detections used")
rejected = fit.residuals[~fit.used, :2] * 206264.8   # their RA, Dec residuals in arcsec
```

### One object, from MPC-style data
```python
import numpy as np
from spacerocks import orbfit
from spacerocks.spice import SpiceKernel
from spacerocks.time import Time

kernel = SpiceKernel.defaults()   # DE440, sb441-n16, Earth orientation, leap seconds

ra = np.radians(data["ra"])       # degrees -> radians
dec = np.radians(data["dec"])
epochs = [Time.from_isot(t.rstrip("Z")) for t in data["obsTime"]]   # UTC

fit = orbfit.fit(ra, dec, epochs, data["stn"], kernel)
print(fit)                         # OrbitFit(flag=0 [converged], chi2=..., ...)
rock = fit.rock                    # a SpaceRock: propagate, observe, ...
sigma = np.sqrt(np.diag(fit.state_covariance))
```

### A survey's worth of objects
```python
res = orbfit.fit_many(data["provID"], ra, dec, epochs, data["stn"], kernel)
good = res["flag"] == 0
states, covs = res["state"][good], res["covariance"][good]
```

### Satellites and roving observers, from MPC ADES data
```python
pos = ades[["pos1", "pos2", "pos3"]].to_numpy(float)          # NaN rows for ground stations
obs = orbfit.observers(ades["stn"], epochs, kernel, sys=ades["sys"], pos=pos)
fit = orbfit.fit(ra, dec, epochs, obs, kernel)
```

### Stellar occultations
```python
occ = ades["ra"].isna() & ades["raStar"].notna()
ra[occ], dec[occ] = orbfit.occultation_radec(np.radians(ades["raStar"][occ]), np.radians(ades["decStar"][occ]),
                                             np.radians(ades["deltaRA"][occ] / 3600), np.radians(ades["deltaDec"][occ] / 3600))
sigma_ra[occ] = np.radians(ades["rmsRA"][occ] / 3600)          # milliarcsecond-level
sigma_dec[occ] = np.radians(ades["rmsDec"][occ] / 3600)
```

### Vereš et al. (2017) weights, from MPC ADES data
```python
sigma = orbfit.veres_sigma(ades["stn"], epochs, catalog=ades["astCat"], program=ades["prog"])
fit = orbfit.fit(ra, dec, epochs, ades["stn"], kernel, sigma_ra=sigma, sigma_dec=sigma)
```

### Debiased astrometry
```python
ra, dec = orbfit.debias(ra, dec, epochs, catalog=ades["astCat"])
fit = orbfit.fit(ra, dec, epochs, ades["stn"], kernel, sigma_ra=sigma, sigma_dec=sigma)
```

### Yarkovsky
```python
fit = orbfit.fit(ra, dec, epochs, data["stn"], kernel,
                 sigma_ra=data["rmsRA"] / 206265, sigma_dec=data["rmsDec"] / 206265,
                 nongrav="A2")
a2, a2_sigma = fit.nongrav[1], fit.nongrav_sigma[1]
```

### A short arc of a distant object
```python
seed = orbfit.bk_iod(ra, dec, epochs, "W84", kernel)        # e.g. two or three DECam nights
fit = orbfit.fit(ra, dec, epochs, "W84", kernel)            # does the same automatically if Gauss fails
```

### Streaks (astrometry plus rates)
```python
k = np.pi / (180 * 3600) * 24                                  # arcsec/hour -> rad/day
fit = orbfit.fit(ra, dec, epochs, data["stn"], kernel,
                 ra_rate=data["raRate"] * k, dec_rate=data["decRate"] * k,
                 sigma_ra_rate=data["rmsRArate"] * k, sigma_dec_rate=data["rmsDecrate"] * k)
```

### Radar, alone or with optical astrometry
```python
# JPL radar astrometry: delay in us, Doppler in Hz at freqTx; NaN in the column a row lacks
fit = orbfit.fit(np.full(n, np.nan), np.full(n, np.nan), epochs, stations, kernel,
                 delay=delay_us * 1e-6, sigma_delay=rms_delay_us * 1e-6,
                 doppler=doppler_hz, sigma_doppler=rms_doppler_hz, frequency=freq_tx,
                 initial=prior)                                 # radar alone needs a prior orbit
fit.residuals[:, 4] * 1e6, fit.residuals[:, 5]                  # delay (us) and Doppler (Hz) residuals
```
With optical astrometry in the same arrays (RA/Dec on those rows, NaN delay/Doppler), the full
pipeline runs without `initial`.

### Refining a known orbit, or checking residuals
```python
fit = orbfit.fit(ra, dec, epochs, "X05", kernel, initial=rock)
resid = orbfit.residuals(rock, ra, dec, epochs, "X05", kernel)[:, :2]   # RA, Dec (radians)
```

### Where to look, and how far to search
```python
fit = orbfit.fit(ra, dec, epochs, stations, kernel)
p = orbfit.predict(fit, [Time.from_isot("2027-03-01T05:00:00")], "W84", kernel)
print(np.degrees(p["ra"]), np.degrees(p["dec"]), p["sigma_major"], p["pa"])
```

### A comet's original orbit
```python
water = [0.1112620426, 2.15, 5.093, 4.6142, 2.808]
fit = orbfit.fit(ra, dec, epochs, stations, kernel, nongrav="A1A2", gofr=water)
o = orbfit.comet_orbits(fit, kernel, gofr=water)["original"]
print(o["inv_a"] * 1e6, "x 1e-6 / AU")     # below ~100: a first-time visitor from the Oort cloud
```

### Keeping a catalog up to date
```python
cat = orbfit.fit_many(ids, ra, dec, epochs, stations, kernel)            # first pass
# ... later, with every detection so far (old and new):
cat = orbfit.fit_many(ids2, ra2, dec2, epochs2, stations2, kernel, prior=cat)
print(collections.Counter(cat["route"]))   # e.g. {'skip': 9120, 'sequential': 804, 'cold': 31, ...}

fit = orbfit.fit(ra, dec, epochs, stations, kernel)                     # one object
fit = orbfit.fit(ra_all, dec_all, epochs_all, stations_all, kernel, prior=fit)
print(fit.route)                                                        # "sequential"
```

<h2 style="border-bottom: 3px solid white;">Notes</h2>

- Epochs given as Julian dates are interpreted in `timescale`, which defaults to UTC.
- Times before 1972 follow UTC's history: its 1960–1972 rate offsets, and before 1960, when
  times are UT, TT − UT = ΔT (about −2 s in 1900, 24 s in 1938, 29 s in 1950; see the time docs).
  SPICE, and so layup, take TT − UTC = 42.184 s for every date before 1972, so old observations
  sit 8–44 s later in layup than they did.
- Ground stations before 1962 Jan 20, where the binary Earth orientation starts, turn with the
  IAU rotation model of the Earth (from `pck00010.tpc`, which `SpiceKernel.defaults()` loads), as
  in sorcha and layup, but corrected for UT1 using ΔT. That matches ITRF93 to 0.3 km where both
  exist. Uncorrected, as in layup, it is off by up to 9 km.
- The fit's epoch is the middle Gauss observation minus the light time. As in layup, that light
  time is the barycentric distance divided by c. With `initial`, the epoch is the rock's.
- Weighting is diagonal, one σ per axis, and σ is whatever you pass: the reported `rmsRA`/`rmsDec`
  (layup's `weight_data="supplied"`), `veres_sigma()` (its `weight_data=True`), or your own. The
  default is 1″. `veres_sigma()` reproduces layup's table, which is a subset of the full Vereš
  et al. (2017) scheme.
- The radar transmitter differs from layup by design. For monostatic radar, layup extrapolates
  the receiving station back to the transmit time as x − vτ − ½aτ². The Taylor expansion is
  x − vτ + ½aτ², so that is a sign error of ~2 km over a 4-minute round trip. Even with the
  sign fixed, a second-order expansion leaves ~0.07 m/s of rotation error in the velocity, about
  1 Hz of Doppler. spacerocks instead computes the station's exact state at the model's
  transmit time when the station is given by MPC code or `Observatory`. On JPL's 2013 radar of
  Apophis, χ² at JPL's own orbit falls from 646 (layup) to 24.5 for 36 measurements, and the
  fitted orbit lands 16 km from JPL's instead of 59 km. With observer arrays, spacerocks
  extrapolates with the sign fixed. To reproduce layup exactly, pass its states as `transmitter`.
- The rate and radar partial derivatives are layup's approximations. The rate Jacobian omits the
  light-time-rate factor, and the radar Jacobian treats the two legs as parallel; the error is of
  order 1e-4 of the partials. The residuals are exact, so this affects only the path to the
  solution, not the solution.
- `iod="herget"`: layup lets an error from its two-body solver escape and end the whole run.
  Here it fails only that starting range, and the next is tried.
- Sequential updates: layup's update fits the state alone and drops a fitted non-grav model. Here
  the parameters the prior fitted are all updated, with its full covariance. For gravity-only
  priors the two are the same. The detection keys hash the fit's inputs (times, astrometry,
  uncertainties, observer states), where layup hashes the reported columns. So new Earth
  orientation data, or new weights, change the fingerprint and trigger a refit.
- `engine="bk_native"`: layup reports a Bernstein–Khushalani fit that does not converge with flag
  2 (its χ² stays infinite, so the χ²/ndof test fires). Here it is flag 1. Inside the pipeline
  both are just "not converged", and the pipeline's flags agree.
- Speed: the full pipeline on layup's 99-object MPC test set takes 6.7 s on 2 cores, against 378 s
  for layup's `orbitfit()`. The Levenberg–Marquardt iterations alone are 2–3× faster, since the
  forward and backward integrations run in parallel. The rest of layup's time is per-object setup.
