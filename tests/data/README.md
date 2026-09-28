# Test data

- `layup_streak_synthetic.json`, `layup_radar_synthetic.json`: noise-free synthetic streak
  (astrometry + sky-motion rates) and radar (delay + Doppler) arcs of a main-belt orbit, generated
  with REBOUND/ASSIST by the layup project (https://github.com/Smithsonian/layup,
  `tests/data/`, commit 578e51c). Copyright (c) 2024 Matthew Holman, MIT License.
  Used by `tests/orbfit_tests.rs` as independent references for spacerocks' rate and radar models.

- `bennu_mpc_2010_2025.csv`: the ground-based optical astrometry of (101955) Bennu from 2010 on,
  from the Minor Planet Center's observations API (https://data.minorplanetcenter.net/api/get-obs,
  retrieved 2026-09-27). Columns: TDB Julian date, RA and Dec (radians), on-sky 1-sigma
  uncertainties (radians; the reported rmsRA/rmsDec, 1" where none, floored at 0.1"), MPC station.
  Used by `tests/orbfit_tests.rs` as a real long NEO arc on which layup's pipeline finds no orbit.
