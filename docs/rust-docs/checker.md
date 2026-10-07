# `spacerocks::checker`

Identifies detections with known objects, like the MPC's MPChecker. The Python interface and
the validation results are in `docs/python-docs/checker.md`.

## Types

- **`Catalog`**: orbits stored as columns.
  - `orbits`: a J2000 / SSB `Population` with the names, epochs (TDB JD) and barycentric
    states, the same type `batch` works on.
  - `nongrav`, `fit_nongrav`.
  - `covariance`: row-major npar×npar (the state, then the fitted A1–A3), empty if unknown.
  - `h`, `g`, `u` (the MPC's uncertainty parameter), `last_obs`.
  - An optional `snapshot`: every orbit integrated to `snapshot_epoch`.
- **`CheckOptions`**: `nsigma`, `radius`, `max_age`, `margin`, `margin_growth`,
  `max_uncertainty`, `floor`, `cross_track`, `window`, `epsilon`, `block`, `batch`,
  `parallel`. Angles are in radians and times in days.
- **`Match`**: one (detection, object) pair.
  - The prediction: `ra`, `dec`, and its `covariance` along (RA·cos Dec, Dec), `floor`
    included.
  - The comparison: `offset` (observed − predicted, tangent plane), `separation`, `distance`
    (Mahalanobis), `consistent`, `log_likelihood`.
  - Also `ra_rate`, `dec_rate`, `delta`, `r_helio`, `mag` and `from_covariance`.
  - `ellipse()` gives the 1σ ellipse.
- **`Checked`**: `matches`, sorted by detection, then with consistent pairs first, then by
  likelihood. Also `failed` (object, reason) and `candidates`, the number of pairs refined.
- **`mpc::MpcElements`, `mpc::MpcOrb`**: records read from the MPC's files.

## Functions

```rust
Catalog::mpcorb(path: Option<&Path>, download: bool, kernel) -> Result<Catalog>
Catalog::from_mpc_elements(&[MpcElements], kernel) / from_mpc_orbs(&[MpcOrb], kernel)
Catalog::from_fits(names, fits, h) -> Result<Catalog>
catalog.make_snapshot(epoch, kernel, &BatchOptions, block)
catalog.save(path) / Catalog::load(path)
catalog.select(&indices) / extend(&other) / states_at(&indices, &targets, kernel, &opts)

checker::check(&catalog, &Astrometry, correlation: &[f64], kernel, &CheckOptions) -> Result<Checked>

mpc::read_mpcorb(path)            // mpcorb_extended.json(.gz) or MPCORB.DAT(.gz)
mpc::parse_mpc_orb(&serde_json::Value)
mpc::mpcorb_path(download, update), mpc::download(url, path), mpc::mpc_dir()
mpc::unpack_epoch, mpc::unpack_designation, mpc::elements_to_state, mpc::runoff_from_u
```

The detections are an `orbfit::Astrometry`: TDB epochs, RA/Dec, σ (with RA's along
RA·cos Dec), and barycentric observer positions. Observer velocities are optional; without them
the rates come out NaN.

## How `check` works

1. **Grouping.** Detections are grouped by epoch and observer. The groups are gathered into
   windows (`window`, 1 day) and then into clusters that span at most 2 × `max_age` days (30 by default).
2. **Coarse search**, over the catalog in blocks of `block` objects, in parallel.
   - Each object is moved by two-body motion about the Sun from its reference state: the
     orbit's epoch or the snapshot, whichever is closer. If neither is within `max_age` days of
     the cluster, the object is first integrated to the cluster's middle, from the closer of the
     two (`batch`, N-body by default).
   - At each cluster, then each window, a bound on how far the object can appear to move is
     checked against a sky grid of that cluster's or window's detections. Only then is each
     detection tested exactly.
   - A detection is a candidate if the object falls within `radius`, or within the `nsigma`
     region, plus `margin + growth × dt²`. Here `growth` is `margin_growth` (1″/day²) for
     perihelia inside 1.3 AU and `margin_growth_distant` (0.02″/day²) beyond.
   - For orbits without a covariance, that region is U's ellipse along the motion. For orbits
     with a covariance, it is a circle of the covariance's rough positional σ (position plus
     4 × velocity × time).
3. **Refinement**, per candidate object, in parallel.
   - The object is predicted at its candidate detections by `orbfit::predict`: ASSIST's force
     model, light time, and the variational equations when there is a covariance. The
     integration starts from the orbit's epoch, or from the snapshot when there is no covariance
     and the snapshot is closer.
   - For orbits without a covariance, U's along-track uncertainty is added, with
     `cross_track` × that across. Then `floor` is added.
   - The Mahalanobis distance uses the prediction's covariance plus the detection's.

## Performance

These figures are for 2 cores, on the full MPCORB (1,439,233 orbits).

| step | time |
|---|---|
| first `Catalog::mpcorb` (parse the JSON, integrate the 13,074 older-epoch orbits, write the cache) | ~40 s |
| later `Catalog::mpcorb` (the binary cache) | 0.4 s |
| `check` of 36 detections over 20 days, near the catalog epoch | 2.9 s |
| `make_snapshot` 16 months from the catalog epoch | 197 s |
| `check` of 67 detections 16 months out, without a snapshot | 189 s |
| the same with a snapshot nearby | 3.5 s |

The coarse search costs about 0.2 µs per object per window (a Kepler step and a grid lookup).
Refinement costs a few milliseconds per candidate object.
