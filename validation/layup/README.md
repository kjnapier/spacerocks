# Cross-check of `spacerocks.orbfit` against layup

These scripts run [layup](https://github.com/Smithsonian/layup) and spacerocks on the same
astrometry and compare the fitted orbits.

## Setup

- Build layup (`pip install -e .`, with the `eigen` and `autodiff` submodules). It needs `rebound`,
  `assist` and `sorcha`.
- layup's cache directory must hold `linux_p1550p2650.440`, `sb441-n16.bsp`, `de440s.bsp`, Earth
  orientation `.bpc` files, `naif0012.tls`, `ObsCodes.json` (layup bundles a gzipped copy) and a
  `meta_kernel.txt`. Without network access, symlink existing files under the names layup expects
  (`earth_620120_260806.bpc`, `earth_2026_260806_2126_predict.bpc`, `pck00010.pck`). layup also
  reads `naif0012.tls` from its *default* cache (`~/.cache/layup`) regardless of the directory
  passed in.
- `SPACEROCKS_KERNELS` must hold `de440.bsp`, `sb441-n16.bsp`, `latest_leapseconds.tls` and the
  Earth orientation files.

## Scripts

| Script | Purpose |
|---|---|
| `run_layup.py in.csv out.npz` | layup's `orbitfit()` per object. Also saves layup's TDB epochs and observer positions. |
| `export_inputs.py out.npz dir` | Writes those inputs (`detections.csv`) and layup's results (`layup_fits.csv`) for the Rust driver. |
| `cargo run --release --example orbfit_layup -- dir/detections.csv sr.csv [spacerocks]` | Fits the same inputs with `orbfit::determine_orbit`. The last argument makes it use spacerocks' own observatory positions. |
| `compare_fits.py layup_fits.csv sr.csv` | Per object: flags, Mahalanobis distance of the states (with layup's covariance), χ², epoch and iterations. |
| `run_spacerocks.py in.csv out.npz` | End to end in Python: spacerocks parses the times, places the observatories and runs `orbfit.fit_many`. |
| `compare_nongrav.py` | A1, A2, A3 and A1A2A3 fits on layup's synthetic Apophis-like arcs (`tests/layup/test_nongrav_a2.py`). |
| `compare_nongrav_auto.py` | `fit_nongrav="auto"`: model choice on layup's synthetic arcs with A1, A2, A3 at 30–1000× its test values, plus 0.1″ noise, with the default thresholds and with the χ²/ndof gate off. |
| `compare_veres.py` | `orbfit.veres_sigma` vs layup's `astrometric_uncertainty_Veres2017` over every station, catalog (names and MPC codes), program and date split the model knows, plus others. |
| `compare_debias.py` | `orbfit.debias` vs layup's `debias` on synthetic bias tables: nside 256 through the binary cache, and nside 64 through the text parser. |
| `compare_observers.py` | `orbfit.observers` vs layup's `obscodes_to_barycentric`: satellites with ICRF_KM/ICRF_AU positions (with and without velocities), roving observers (WGS84), ground stations, and ground stations before 1962 (IAU_EARTH). |
| `compare_per_arc.py` | Per-arc non-grav: layup's piecewise-amplitude arcs (A2 equal and distinct, A1, A3, A1A2A3; noise-free and 0.1″ noise), gravity fit then per-arc fit in both codes. |
| `compare_sequential.py run.npz` | `orbfit.sequential_update` vs layup's `run_sequential_update`: each object's last 10% of detections as the update to a fit of the rest, plus each code's full refit. |
| `compare_incremental.py in.csv [n]` | `fit_many(..., prior=...)` vs layup's `incremental_orbitfit` on a previous/current pair of data sets (a quarter of the objects unchanged, appended to, with a detection removed, and new). |
| `compare_bk_engine.py run.npz ...` | The Bernstein–Khushalani engine alone: `fit(initial=seed, engine="bk_native")` vs `run_bk_native_fit`, from layup's BK IOD seed over all detections. |
| `compare_predict.py run.npz` | `orbfit.predict` vs `predict_sequence`, same orbit and covariance, at every detection and 12 future geocentric epochs. |
| `compare_comet.py [n]` | Original 1/a of the CODE catalogue sample in layup's tests, from spacerocks and layup's `comet`, against the catalogue. |
| `showcase.py mpc_cache out.json [jpl.json]` | Five real objects (Holman, Eros, Bennu, Apophis, Sedna) from the notebook's MPC cache, fitted by both codes and refitted at one epoch, against the MPC's orbit (and JPL's states when given). |
| `compare_herget.py run.npz ...` | `orbfit.herget_iod` vs layup's `herget_iod` on every object's primary arc. For the whole pipeline, run `run_layup.py` and `run_spacerocks.py` with `IOD=herget`. |
| `compare_bk.py run.npz ...` | `orbfit.bk_iod` vs layup's `run_bk_iod` on every object's primary arc. |
| `make_short_arcs.py out.csv` | 90 synthetic short arcs (TNOs, Centaurs, main belt; 1–4 nights over up to 14 days), which exercise the Bernstein–Khushalani fallback. |
| `compare_rates_radar.py` | Streaks and radar: layup's synthetic streak and radar arcs, real JPL radar of Apophis (2013), layup's streak ADES file, and both full pipelines on streaks and on Apophis optical + radar. |
| `bench.py out.json run.npz ...` | Per-object timing of both pipelines, and of the LM alone, with setup excluded (see Benchmark). |
| `layup_self_consistency.py run.npz out.npz` | Runs layup's `do_fit` twice per object, with RA/Dec differing in the last bit (`x*pi/180` vs `np.radians(x)`), to show which outcomes are reproducible even within layup. |

All scripts convert degrees to radians as `x * pi / 180`, which is what layup's `_orbitfit` does.
`np.radians` can differ in the last bit, and on degenerate arcs that is enough to change layup's
own result (see below).

## Results (2026-09-27, layup `578e51c`)

Data sets: layup's demo (1979 HP, 4135 detections); `tests/data/100_random_mpc_ADES_provIDs_no_sats.csv`
(99 MPC objects, 50,714 detections); and a mixed set of 58 objects. The mixed set is 55 Rubin short
arcs from `x05_short_arcs.obs80`, 3I/ATLAS, and the Bernstein et al. KBOs minus their HST rows.

| | Flags equal | Epoch equal | Iterations equal | State difference (σ) |
|---|---|---|---|---|
| layup's inputs, 99 MPC objects | 99 / 99 | 98 / 98 | 98 / 98 | ≤ 3.4e-7 |
| layup's inputs, mixed set | 58 / 58 | 57 / 57 | 51 / 57 | ≤ 4.5e-4 (one object; the rest ≤ 6e-7) |
| end to end, 99 MPC objects | 99 / 99 | | | ≤ 3.7e-5 |
| end to end, mixed set | 58 / 58 | | | ≤ 4.5e-4 |
| end to end, 1979 HP | 1 / 1 | | | 2.8e-5 |

On layup's inputs, χ² agrees to 1e-12. The end-to-end differences come from observatory positions
(median 0.5 m, at most 0.1 km) and from rounding the times to a single f64 Julian date (≤ 80 µs).

The iteration differences on short arcs come from layup's LM acceptance rule, which spacerocks
keeps. Near the minimum, whichever code's integration noise first raises χ² stops there, so the
two can stop at slightly different points. On the one object where this shows (K21RO0K), the
spacerocks χ² is the lower of the two.

### Bernstein–Khushalani fallback

`bk_iod` matches layup's `run_bk_iod` to ≤2e-8 relative on all 156 primary arcs of the real data
sets, and both reject the same one. None of those objects needs the fallback. On the synthetic
short arcs, layup takes the BK path for 8 of 90 objects.

| Synthetic short arcs (90) | |
|---|---|
| Flags equal | 88 / 90 |
| Converged in both | 71, of which 67 agree to ≤5e-6σ |
| Differences | Only on 1- and 2-night arcs: 4 states 0.1–0.72σ apart and 2 flag mismatches |

The differences are not a porting error:

- **tno-2n-5.** layup converges from the BK seed, spacerocks does not. With `np.radians`
  inputs, layup itself also does not converge: its fit ends at the same χ² (0.0469593) as
  spacerocks'.
- **The four states.** On single- and two-night arcs the χ² valley is essentially flat; the "orbits"
  have main-belt test objects at 4.4 AU moving at 100 km/s. layup's LM stops at the first
  noise-induced χ² increase (see the note on its acceptance rule), so where it stops depends on
  integrator round-off. The two codes stop at different points along the same valley, within
  0.72σ of each other. spacerocks ends with the lower χ² in 3 of the 4.
- **cen-1n-3.** layup's end point is flagged implausible (>200 km/s excess); spacerocks' has
  a lower χ² and passes.

A noise-free case where Gauss finds no root at all (a 40 AU TNO over two nights) is in
`tests/orbfit_tests.rs`. layup and spacerocks both take the BK path there and converge to the
same orbit (3e-11σ apart, distance correct to 1e-5).

Non-grav fits (`compare_nongrav.py`) give the same A1, A2, A3 to six digits, the same formal
uncertainties, the same iteration counts, and states within 6e-7σ.

Scaled convergence (`CONV_FRAC=0.1` for both `run_layup.py` and `run_spacerocks.py`, 99 MPC
objects): 0 flag mismatches, states ≤3.8e-5σ apart, and the same iteration count for every object.
Summed over the objects, the final fits take 145 iterations instead of 259, in both codes.

Vereš et al. (2017) weights: `compare_veres.py` finds no difference in 74,088 combinations of
station, catalog, program and date, including dates on the split boundaries. End to end
(`WEIGHT_DATA=veres` for both runners, 99 MPC objects): 0 flag mismatches, states ≤6.6e-5σ
apart, χ² within 2.5e-7.

Debiasing (`compare_debias.py`): JPL's table isn't reachable from the sandbox, so both codes read
the same synthetic tables, with random offsets and proper motions per catalog and pixel. On 200,084
positions at nside 256, including the poles, the polar-cap boundary and RA 0/360, the results
agree to 2.3e-6 mas; with the nside-64 text file, 3.0e-5 mas. Every HEALPix pixel lookup therefore
agrees with healpy's, since a wrong pixel would be off by ~0.5″.

Observer positions (`compare_observers.py`, the same TDB instants for both) agree to ≤0.6 m and
≤5e-5 m/s for satellites (ICRF_KM, ICRF_AU, with and without velocities), roving observers
(WGS84), and ground stations. The 0.6 m floor is the precision of a Julian date. For ground
stations between 1901 and 1962, both codes use the IAU rotation model from `pck00010.tpc`, and
the layup cache's meta-kernel has to include it. spacerocks corrects that model for UT1
(ΔT), and layup does not. The stations therefore differ by up to 7.7 km, and match to 0.58 m
once layup's are turned by the same correction. Given UTC instead of TDB, times before 1972 also
differ. SPICE takes TT − UTC = 42.184 s there. spacerocks follows the 1960–1972 UTC table
and, before 1960, ΔT (see `validation/time/compare_erfa.py`).

Per-arc non-grav (`compare_per_arc.py`, 10 cases): flags and χ² equal. Arc A and arc B
amplitudes agree to ≤1.3e-7σ and states to ≤4.3e-7σ, with identical formal uncertainties.

Herget IOD (`compare_herget.py`, the 99 primary arcs): both codes find a seed for the same 97
arcs. The states agree to 8e-14 relative (median). The largest difference, 7e-4 on 244185, is
within the method's 0.003 AU stopping tolerance. layup's two-body solver raised
`KeplerConvergenceError` four times (48550 and 9946 from 2 AU, 415974 from 2 and 5 AU). layup
lets that escape and end the run; the scripts, like spacerocks, move on to the next starting
range. The whole pipeline with `IOD=herget` gives 0 flag mismatches, 97 converged in both,
states to ≤2.9e-5σ (median 2.6e-6σ), and χ² to ≤3.6e-7 relative. It runs in 3.1 s against 322 s
for layup. On the primary arcs alone, `herget_iod` takes 1.1 s against layup's 15 s.

Sequential updates (`compare_sequential.py`, 98 objects with a converged prior): 0 flag
mismatches. The updated states agree to ≤5.9e-7σ (median 4e-8σ) and χ² to 2.5e-8 relative, with
equal ndof. The gate's measure of the move agrees to 4e-8 prior σ. Both codes' updates come
within ≤6.0e-4σ of their own full refits (median 7e-7σ); `fit(..., prior=)` routes all 98
through the sequential update. Incremental fitting (`compare_incremental.py`, 40 objects): the
same routes in both codes (10 skip, 10 sequential, 10 full, 10 cold), 0 flag mismatches, and
final states within ≤2.4e-5σ (median 5e-6σ).

Bernstein–Khushalani engine:
- End to end on the 99 MPC objects (`ENGINE=bk_native` for both runners): 0 flag mismatches,
  and 80 converge in both (the BK engine converges on fewer of these long arcs than the
  Cartesian one, in both codes). States agree to ≤3.0e-5σ and χ² to ≤3.6e-7 relative. The run
  takes 42 s here and 423 s in layup.
- End to end on the synthetic short arcs, layup stalls: IAS15 ground on one main-belt arc for
  over half an hour.
- The engine alone (`compare_bk_engine.py scratch_short.npz`, 80 short arcs from layup's BK
  seed):
  - Converged in both: 75, with states within a median of 2.4e-9σ (max 9e-4σ) and χ² within
    4.9e-5.
  - Three more converge in both codes, but spacerocks' pipeline check flags them 9: a
    hyperbolic excess speed above 200 km/s, a check layup's engine doesn't make.
  - One (cen-2n-2) converges in layup at 74 iterations and not in spacerocks within 100.
  - Iteration counts are equal in 51 of 76. These one- and two-night arcs are nearly
    degenerate, so the two integrators' round-off steers the last steps.

Predictions (`compare_predict.py`, 98 objects, 51,721 predictions): directions agree to
≤2.6e-5 mas (median 1.5e-6), and the 2×2 sky covariances to ≤2e-11 relative. The 1σ major axes
range from 0.014″ to 1″.

Comet original orbits (`compare_comet.py`, 334 CODE catalogue comets with an uncertainty):
- |spacerocks − CODE| / σ has median 0.0041, 90th percentile 0.070 and maximum 5.5. layup's is
  the same: 0.0042, 0.070 and 5.5.
- spacerocks and layup agree to ≤0.001×10⁻⁶/AU (median 0.0002).
- The two comets that layup's test sets aside for non-gravitational forces are 22.5σ and 7.3σ
  off in both codes.
- The run takes 3.2 s here and 316 s in layup.
- With `de440s`, which starts in 1849, every comet stops short of 250 AU, at 92–168 AU. The
  elements taken there change 1/a by a median of 0.01×10⁻⁶/AU and at most 0.1×10⁻⁶/AU.

Automatic model selection (`compare_nongrav_auto.py`) makes the same choice as layup's
`_select_nongrav_auto` in 42 of 42 cases. The choices cover gravity only, A1, A2 and A3, and each
case runs with the default gate and with the gate off. The adopted parameters agree to ≤1.2e-7σ
and χ² to the digits printed.

### Rates and radar (`compare_rates_radar.py`)

The model is compared by χ² at fixed states: a 1-iteration fit returns χ² at its starting state.
The fits are compared from the same perturbed start.

| Case | χ² at the perturbed start | Fit, layup vs spacerocks |
|---|---|---|
| Synthetic streaks (7 detections, 28 rows) | agree to 1.3e-10 | both converge, 1.4e-8σ apart |
| ... one rate corrupted by 50σ | agree to 1.3e-10 | both flag 2, χ² 2500.0 |
| Synthetic radar, delay + Doppler | agree to 4.0e-10 | both converge, 1.2e-4σ apart |
| ... delay only / Doppler only | agree to 4.0e-10 / 1.5e-11 | 1.8e-6σ / 8.7e-9σ apart |
| Apophis 2013, real radar, layup's stations and transmitter model | agree to 6e-5\* | 3.5e-3σ apart |
| Full pipelines: synthetic streaks | | same epoch, flag and ndof; 8.6e-9σ |
| Full pipelines: Apophis synthetic optical + real radar (126 rows) | | same epoch, flag and ndof; 1.5e-3σ |
| layup's streak ADES file (2 objects, 4 detections each) | | both flag 5 on both |

\*The transmitter states were passed explicitly, with τ from the JPL orbit, so the model at the
perturbed start differs slightly.

**The radar transmitter.** layup places the monostatic transmitter at x − vτ − ½aτ², but the
Taylor expansion of x(t − τ) is x − vτ + ½aτ². On Apophis (τ ≈ 4 min) the sign error puts the
transmitter ~1.7 km from the station's true position. Even with the sign fixed, the expansion
leaves ~0.07 m/s of Earth-rotation error in the transmitter's velocity, about 1 Hz of Doppler
against uncertainties of 0.1–0.2 Hz. spacerocks computes the station's exact state at the model's
transmit time (`transmitter_site`) when the station is known by code, which is always the case
for real radar.

| Apophis 2013, 36 radar measurements | χ² at JPL's orbit | Fitted χ² (ndof 30) | Fit vs JPL Horizons |
|---|---|---|---|
| layup | 645.8 | 5.39 | 59.5 km, 3.2 cm/s |
| sign fixed, second-order extrapolation | 37.9 | 5.70 | 59.5 km, 3.2 cm/s |
| spacerocks (exact transmitter state) | 24.5 | 3.87 | 16.1 km, 0.59 cm/s |

JPL's orbit is fitted to these same radar data together with the optical astrometry, so the
right model should find it consistent. spacerocks' residuals there are 0.23 µs rms in delay and
0.06 Hz in Doppler.

The rate and radar partials are layup's approximations (see `docs/rust-docs/orbfit.md`) and are
tested to that accuracy. The residuals are exact.

### Benchmark (`bench.py`, 2026-09-27)

Both codes are called from Python per object, on identical inputs: TDB epochs, barycentric
observer positions, and RA/Dec as layup converts them. Setup is excluded on both sides; layup's
`get_ephem` is memoized, though it costs only ~1 ms anyway. Each time is the best of 3 on a 2-core
sandbox (x86-64).

- **Pipeline** is layup's `do_fit(iod="auto")` against `orbfit.fit`: IOD, candidate screening,
  fit and build-up.
- **LM** is `run_from_vector_with_initial_guess` against `fit(initial=...)`, starting from layup's
  answer perturbed by 1e-6 AU and 1e-8 AU/d.

| Set (objects) | Pipeline, layup | spacerocks, 1 thread | Ratio | spacerocks, 2 threads | Ratio |
|---|---|---|---|---|---|
| 99 MPC objects | 41.1 s | 13.0 s | 3.2× | 10.6 s | 3.9× |
| 57 mixed (Rubin short arcs, 3I/ATLAS, KBO) | 7.4 s | 3.3 s | 2.2× | 1.9 s | 3.8× |
| 90 synthetic short arcs | 19.2 s | 10.2 s | 1.9× | 6.3 s | 3.1× |
| 1979 HP (4135 detections) | 0.55 s | 0.14 s | 3.8× | 0.13 s | 4.2× |
| **All 247** | **68.2 s** | **26.7 s** | **2.6×** | **19.0 s** | **3.6×** |

| LM alone (228 fits) | layup 14.2 s | spacerocks 6.7 s (1 thread) | 2.1× | 5.0 s (2 threads) | 2.8× |
|---|---|---|---|---|---|

- On the MPC set the per-object single-thread ratio is tight: median 3.3×, 10th–90th percentile
  2.6–3.7×. It narrows to about 2× on short arcs, where fixed per-evaluation costs such as building
  the simulation weigh more.
- The single-thread pipeline gain exceeds the LM gain because spacerocks' candidate prefilter
  integrates without variational particles. They aren't needed there, and they don't change the
  step sizes, so the result is identical.
- The second thread helps least on long MPC arcs (1.2×), because most detections fall on one side
  of the epoch and the forward and backward passes are unbalanced. `fit_many` also parallelizes
  across objects, which scales with cores.
- The slowest per-object ratios (0.3–0.9×, i.e. spacerocks slower) are all 1-night arcs. There
  layup's LM stops early at a noise-induced χ² uptick, while spacerocks keeps accepting steps down
  the flat valley (up to 98 iterations vs 14) and ends at lower χ².
- Iteration counts in the LM test agree for all 98 MPC objects. They differ on many short arcs, for
  the stopping-rule reason above.
- Flags agree except for tno-2n-5 (the one-ulp case). layup's `do_fit` omits the >200 km/s check,
  which lives in its `_orbitfit` wrapper, so that check is left out of this comparison.

For scale, layup's full `orbitfit()` takes about 3 s per object (378 s for the 99 MPC objects),
almost all of it outside `do_fit`.

The `x05` rows had to be rewritten with zero-padded seconds before the comparison. layup's Obs80
reader writes times like `09:09:7.2576`, and `_orbitfit` sorts by the `obsTime` string. That puts
detections out of time order and changes the Gauss triplet.
