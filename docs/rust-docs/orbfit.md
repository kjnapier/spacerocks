# Orbit fitting (`spacerocks::orbfit`)

A port of [layup](https://github.com/Smithsonian/layup)'s orbit fitter. The inputs are plain
arrays and the result is a small struct. The Python bindings (`spacerocks.orbfit`) are thin
wrappers; see `docs/python-docs/orbfit.md` for the pipeline and the flags.

```rust
use spacerocks::orbfit::{determine_orbit, observer_positions, Astrometry, FitOptions};

// epochs: TDB Julian dates; ra, dec: radians; codes: MPC observatory codes
let observers = observer_positions(&codes, &epochs, &kernel)?;          // barycentric J2000, AU
let astrometry = Astrometry::new(epochs, ra, dec, observers)?          // default sigma 1"
    .with_sigma(&sigma_ra, &sigma_dec)?;                               // optional, radians
let fit = determine_orbit(&astrometry, None, &kernel, &FitOptions::default());
if fit.converged() {
    println!("{:?} at {} TDB, chi2 {} / {}", fit.state, fit.epoch, fit.chi2, fit.ndof);
}
```

## Data

- **`Astrometry`**: parallel `Vec`s with one entry per detection, laid out like an ADES table with
  nullable columns (NaN where not measured):
  - always: `epoch` (TDB JD; radar: receive time) and `observer` (`[f64; 3]`, barycentric J2000
    AU);
  - optical: `ra`, `dec`, `sigma_ra`, `sigma_dec` (radians; `sigma_ra` on-sky);
  - rates: `ra_rate` (on-sky, cos Dec dRA/dt), `dec_rate`, and their sigmas (radians/day);
  - radar: `delay` (round trip, days), `doppler` (round-trip range rate, AU/day), and their
    sigmas;
  - observer motion, for rates and radar: `observer_velocity`, `observer_acceleration`;
  - the radar transmitter: `transmitter` (its state at the transmit time) or `transmitter_site`
    (its Earth-fixed position).

  The optional columns are empty or full length. Builders: `new`, `with_sigma`, `with_rates`,
  `with_radar`, `with_observer_velocity`, `with_observer_acceleration`, `with_transmitter`,
  `with_transmitter_site`, and `from_observations`. Helpers: `rows(i)` (the `RowKind`s a detection
  contributes), `n_rows()`, `sigma(i, kind)`, `optical()`, `subset(&idx)`, `time_order()`,
  `rho_hat(i)` and `tangent_basis(i)`.
- **Orbits** are `(epoch, [f64; 6])`: a TDB Julian date and a barycentric J2000 state in AU and
  AU/day. Non-gravitational parameters are `[f64; 3]` (A1, A2, A3 in AU/day²).
- **`OrbitFit`**: `epoch`, `state`, `nongrav`, `fit_nongrav`, `covariance` (row-major
  `npar × npar`: the state, then the fitted non-grav parameters), `npar`, `chi2`, `ndof`, `niter`,
  `flag: FitFlag`, `residuals` (`[ra, dec, ra_rate, dec_rate, delay, doppler]` per detection, NaN
  where not measured), and `used` (the detections in the final fit; empty if there is no orbit).
  Methods: `converged()`, `state_covariance()`, `nongrav_sigma()`.
- **`FitOptions`**: layup's settings by default. `max_iter` 100, `screen_iter` 80, `tolerance`
  1e-12, `conv_frac` 0 (layup's scaled convergence test: a step converges below
  `max(tolerance, conv_frac * sigma_i)` once one has been accepted), `chi2_threshold` 10, IAS15 `epsilon` 1e-9, `fit_nongrav`, `nongrav_auto`
  (`Option<NongravAuto>`: layup's `fit_nongrav="auto"`, with its thresholds 1.5, 9 and 3), `gofr`, `arc_gap` 90 d,
  `min_distance` 0.3 AU, `prefilter_sigma` 1000, `bk_fallback` (true), `herget` (false: seed from Herget's method
  instead of Gauss, with no Bernstein–Khushalani fallback, as layup's `iod="herget"`), and `parallel` (run the
  forward and backward integrations concurrently). Robust mode (not in layup): `robust` (false),
  `outlier_sigma` 4, `seed_window` 60 d, `seed_tries` 10. Sequential updates: `max_update_sigma` 4.
  `engine`: `Engine::Cartesian` (default) or `Engine::BkNative` for the pipeline's gravity-only
  fits, including the fit from an `initial` orbit (layup's `engine`; non-grav fits stay Cartesian,
  and layup's `initial_guess` path always is).

## Functions

| Function | What it does |
|---|---|
| `determine_orbit(&a, initial, &kernel, &opts)` | The whole pipeline: IOD, candidate screening, the fit, the build-up fallback, and optional non-grav fitting. With `initial = Some((epoch, state, nongrav))`, only the fits run. With `opts.robust`, outlier rejection and the short-window fallback (below). |
| `select_nongrav(&a, &gravity_fit, &kernel, &opts, &thresholds)` | layup's `_select_nongrav_auto`: the non-grav model ladder from a converged gravity-only fit. `determine_orbit` runs it when `opts.nongrav_auto` is set. |
| `fit_orbit_bk(&a, epoch, &state, &kernel, &opts, max_iter)` | layup's `run_bk_native_fit`: LM in Bernstein–Khushalani parameters (`bk_fit::cartesian_to_bk`, `bk_to_cartesian`, `dcart_dbk`, `choose_fiducial`), with the Cartesian partials chained through `dcart_dbk` and a fixed bound-orbit prior on gdot (variance `2 mu gamma^3 - adot^2 - bdot^2`, mu = layup's MU_SUN). χ² includes the prior; the covariance is carried back to Cartesian. Rate and radar rows are fitted too (layup: RA/Dec only). |
| `predict(&fit, &epochs, &observers, &kernel, &opts)`, `predict_astrometry(...)` | layup's `predict_sequence`: light-time-corrected RA/Dec and distance with the fit's covariance mapped to the sky (`Prediction::cov`, `ellipse()`), forward and backward passes from the epoch. Uses the fit's full covariance and non-grav (layup: the state's). |
| `comet_orbit(epoch, &state, &nongrav, future, reference, &kernel, &opts)`, `original_and_future(...)`, `comet::barycentric_elements` | layup's `comet`: barycentric osculating elements (GM of the Sun and planets, `MU_TOTAL`) at the inbound or outbound crossing of `reference` AU (250), found by marching and bisection. If the ephemeris ends first, the elements where it stopped (`reached = false`). |
| `sequential_update(&new, &prior, all, &kernel, &opts)` | layup's `sequential_update`: fits `new` alone at the prior's epoch with the prior's covariance as a Gaussian prior on the parameters it fitted. If that fails or moves more than `opts.max_update_sigma` prior sigmas, it refits `all` (when given) from the prior, or flags 7/8. Returns `(fit, accepted)`. |
| `update_orbit(&current, prior, &kernel, &opts)` | layup's `incremental_orbitfit` for one object, with `prior: Option<&PriorFit>` (a fit plus the detection keys it covers). Returns `(fit, UpdateRoute)`: `Skip`, `Sequential`, `SequentialFallback`, `Full` or `Cold`. |
| `fit_orbit_with_prior(&a, epoch, &state, &nongrav, &info, &kernel, &opts, max_iter)`, `update_mahalanobis(&prior, &fit)` | The LM fit with a Gaussian prior given by its information matrix (`L^T` rows added to the QR, steps judged on χ² plus the prior term, ndof = rows), and the move in prior sigmas. |
| `Astrometry::detection_key(i)` / `detection_keys()` / `fingerprint()` | FNV-1a hash of every input of a detection, and an order-independent fingerprint of a set (`orbfit::fingerprint(&keys)`). |
| `fit_orbit(&a, epoch, &state, &nongrav, &kernel, &opts, max_iter)` | One Levenberg–Marquardt differential correction (layup's `run_from_vector_with_initial_guess`). |
| `fit_orbit_per_arc(&a, epoch, &state, &nongrav, &nongrav_arc2, &kernel, &opts, max_iter)`, `residuals_per_arc(...)` | Per-arc non-grav (layup's `per_arc`): separate parameters for the detections before (arc A) and after (arc B) the epoch, one state. `FitOptions::nongrav_per_arc` makes `fit_orbit` and `determine_orbit` do this. `OrbitFit::nongrav_arc2`, `nongrav_arc2_sigma()`. |
| `residuals(&a, epoch, &state, &nongrav, &kernel, &opts, partials)` | Residuals and, optionally, the `rows × npar` Jacobian from the variational equations. Rows are packed per detection, with `row_detection` and `row_kind` labelling each; `per_detection(n)` unpacks them. |
| `veres_sigma(station, jd_tdb, catalog, program)` | One-sigma astrometric uncertainty (arcsec) after Vereš et al. (2017), exactly as layup's `astrometric_uncertainty_Veres2017` assigns it. |
| `BiasTable::load(path)` / `BiasTable::load_default(download)`, `table.debias(ra, dec, jd_tdb, catalog)` | Star-catalog debiasing (Eggl et al. 2020) exactly as layup's `debias`: nested HEALPix lookup as healpy does it, the offset plus proper motion since J2000, and layup's unit-vector round trip. The table is JPL's `bias.dat`. It is parsed once into a binary `bias.bin` and memory-mapped. |
| `bk_iod(&a, epoch)` | Bernstein–Khushalani linear initial orbit (layup's `run_bk_iod`): barycentric state at `epoch`, or `None`. It is the pipeline's fallback when no Gauss candidate converges (`FitOptions::bk_fallback`, default on, as in layup's `iod="auto"`). |
| `herget_iod(&a, &idx, &kernel, &opts)`, `herget(&a, &idx, rho0, &kernel, &opts)` | Herget's method (layup's `herget_iod` and `herget_with_assist`) on the detections `idx` in time order: returns the epoch of the first and the barycentric state there, or `None`. `herget_iod` starts from ranges of 2, then 5, then 40 AU. `FitOptions::herget` makes the pipeline use it on the longest arc. |
| `gauss_states(rho, observer, t, mu, min_distance)` | Gauss's method on three unit vectors, observer positions and epochs. `gauss_astrometry` picks the three from an `Astrometry`, and `gauss` takes `Observation`s. |
| `build_sequence(&epochs, gap)` | Arcs in fitting order. |
| `gauss::select_triplet(&epochs, &idx, target_days)` | layup's span-targeted Gauss triplet. |
| `observer_positions` / `observer_states` / `observer_accelerations` (`&codes, &epochs, &kernel`) | Barycentric positions, velocities and accelerations of MPC observatories, computed in parallel. |
| `occultation_radec(ra_star, dec_star, delta_ra, delta_dec)` | RA/Dec of an occulting object from ADES occultation astrometry (radians), applying the offset on the tangent plane at the star. |
| `ades_observer_state(sys, ctr, pos, vel, jd_tdb, &kernel)` | Barycentric state of an observer given with its detection (ADES `sys` ICRF_KM, ICRF_AU or WGS84; `ctr` 399), as layup places it. `geodetic_to_earth_fixed(lon, lat, h)` converts WGS84 to ITRF93. |
| `observing::earth_fixed_state(&p_itrf, jd, &kernel)` | Barycentric state of an Earth-fixed point: the radar transmitter's exact state. |

## How a residual evaluation works

For a trial orbit, `residuals` builds a `SpiceSimulation` at the orbit's epoch. The simulation
holds the object, its six state variational particles, and one particle per fitted non-grav
parameter. It is integrated forward through the later detections and, in a second simulation,
backward through the earlier ones. States are read from IAS15's dense output (`integrate_rel`),
so each pass is a single integration.

Light time comes from four fixed-point iterations starting from zero delay. The residual is the
observed minus computed direction projected on the observed tangent plane. The partials include
the light-time term, as in layup's `compute_optical_residuals`.

The integrator matches layup's: IAS15 with the PRS23 step criterion (layup sets REBOUND's
`adaptive_mode = 2`), ε = 1e-9, a first step of 0.001 d, and ASSIST's default forces.

## Robust mode

`FitOptions::robust` wraps layup's pipeline in stages. The first stage that *settles* wins.
Settling means the final fit's detections are exactly those within `outlier_sigma` of its orbit.
If the set still flips after ten rounds, an average χ² per degree of freedom within
`outlier_sigma²` also counts.

1. Layup's pipeline, or the fit from `initial`, then outlier rejection. Each rejection round
   refits the detections within `outlier_sigma` from the current orbit, re-evaluating every
   detection, rejected ones included.
2. Orbits grown from short windows. The detections are split into windows of at most
   `seed_window` days, and up to `seed_tries` windows are tried, the most detections first.
   Layup's pipeline runs on each window. The first that converges is widened by its own span on
   each side, repeatedly, with a refit and rejection at each step, until it covers everything.
3. Stages 1 and 2 again with layup's χ² test on initial orbits lifted. That test fails every
   candidate when the arc used for the initial orbit holds an outlier.

If nothing settles, layup's own result is returned. Non-gravitational parameters are fitted after
the gravity-only stage settles. The joint fit starts both from every detection and from the
surviving ones, and keeps whichever settles on more detections. On Apophis, the gravity-only
orbit otherwise rejects the precise radar that the Yarkovsky term explains.

On real MPC data (`notebooks/examples/orbit_fitting.ipynb`, and `tests/data/bennu_mpc_2010_2025.csv`):

- Bennu (2010–2025, 310 detections): layup's pipeline gives flag 3; robust mode converges with
  308 detections, a = 1.1260 AU.
- Apophis, optical and radar with A2 (8212 detections): A2 = −2.913e-14 ± 2.0e-16 AU/day². JPL
  has −2.902e-14 ± 1.9e-16.
- With data from 2005 on, `fit_many` of Holman, Eros, Bennu, Sedna and Apophis: layup's pipeline
  loses Bennu and Apophis (flag 3); robust mode fits all five in 6 s.

With robust mode off, nothing changes: the 99-object layup comparison is identical.

## Differences from layup

- Detections are ordered by epoch. layup sorts them by the `obsTime` *string*, which misorders
  times whose seconds are not zero-padded, and its own Obs80 reader writes those
  (`09:09:7.2576`).
- Near the minimum, both codes can stop at slightly different points, because layup's
  acceptance rule (kept here) is noise-limited. It compares χ² at the current point with χ² at the
  previous point, and after a rejected step it never re-evaluates. So once integration noise makes
  χ² tick up, λ doubles until the step drops below the tolerance. On layup's test sets this moves
  states by up to 5e-4σ and never changes a flag.
- `Astrometry::default` uncertainties are 1/206265 rad, like layup's. `ARCSEC` is the exact value.
- Rates and radar use layup's models (`compute_streak_residuals`, `compute_radar_residuals`)
  and its approximate partials, except for where the radar transmitter is at the transmit time.
  layup extrapolates the receiving station as x − vτ − ½aτ², which has the wrong sign on the
  acceleration term, and a second-order expansion misses ~1 Hz of Doppler even with the sign
  fixed. Given `transmitter_site`, which the Python layer fills for stations given by code,
  spacerocks evaluates the station's state exactly. Without a site it extrapolates with the sign
  fixed. See `validation/layup/README.md`.
- IOD uses optical detections only. layup would pass radar rows to Gauss's method if they fell
  in the triplet.
- Before 1972, times follow UTC's actual history (see `docs/rust-docs/time.md`): the 1960–1972
  rate offsets, and before 1960 TT − UT = ΔT. SPICE, and so layup, take TT − UTC = 42.184 s for
  every date before 1972, which puts observations 8–44 s off in the 20th century. Ground
  stations before 1962 Jan 20 use the IAU rotation model of the Earth (`pck00010.tpc`), as
  sorcha and layup do, but corrected for the Earth's actual rotation
  (`spice::iau_earth_rotation_delay`, from ΔT). That matches ITRF93 to 0.3 km over 1962–2019,
  where the uncorrected model (layup's) is off by up to 9 km.
- Occultations: layup means to take occultation positions from the star plus offset
  (`_use_star_astrometry`), but that path never runs. It tests `ra is None`, where a missing
  value is NaN, and it reads `deltRA` instead of ADES `deltaRA`. `occultation_radec` does it,
  exactly on the tangent plane.
- Herget: layup lets a `KeplerConvergenceError` from its two-body solver escape `herget_iod`,
  which ends the whole `orbitfit` run (and leaves the detection times light-time shifted).
  Here the error fails only that starting range, and the next one is tried.
- Sequential updates cover every parameter the prior fitted (layup's fit the state alone). The
  detection keys hash the fit's inputs, not the reported ADES columns.
- `Engine::BkNative`: a Bernstein–Khushalani fit that doesn't converge is flag 1 here; layup
  reports 2, since its χ² stays infinite and trips the χ²/ndof test.

## Herget's method

`herget` guesses the ranges `rho_1`, `rho_n` to the first and last detections, which fixes the
positions `r_1`, `r_n`. A shooting method on two-body motion (universal variables, Danby's
formulation as layup ports it) finds the velocity that carries `r_1` to `r_n`. The two-body state
transition matrix gives how that velocity changes with each range. Two ASSIST integrations,
one from each end with a variational particle along that change, give each detection's
predicted direction and its derivative with respect to each range. The 2×2 normal equations then
correct both ranges (each capped at half its value), with the detection times corrected by the
mean light time. The iteration stops when the mean correction is below 0.003 AU, and fails
after 100 iterations or if the normal equations become degenerate (median normalized
determinant ≤ 1e-3). The state returned is at the first detection. It is an initial orbit
(good to about 0.01 AU in position), so the pipeline refines it as it would a Gauss candidate.

## The Bernstein–Khushalani fallback

`bk_iod` picks a right-handed frame `(a, b, n0)` with `n0` along the mean line of sight. Each
detection's gnomonic coordinates `(x, y)` in that frame then satisfy

```text
x = alpha + adot t - gamma (X - x z),   y = beta + bdot t - gamma (Y - y z)
```

Here `(X, Y, z)` is the observer's position in the frame, and `t` is the time from the epoch
corrected by `z / c`. This is the exact gnomonic relation multiplied through by its perspective
denominator, with the line-of-sight velocity `gdot` fixed at 0. The weighted least-squares
solution is refined once with weights divided by the squared denominator, the factor
`s = sqrt(1 + alpha² + beta²)` is divided out of `gamma` and the dots, and the result is converted
to a Cartesian state.

In the pipeline, the fallback runs on the longest arc at its middle detection, after every Gauss
candidate has failed, and it also runs when Gauss finds no root at all. On the real test sets it
agrees with layup's `run_bk_iod` to ~1e-8 relative. See `validation/layup/README.md` for the
end-to-end comparison.
