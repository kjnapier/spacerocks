# Work in progress (not compiled)

Source files that were in `src/` but not part of the module tree, moved here during the
2026 cleanup so the crate only contains code that builds. Nothing here is compiled or tested.

- `orbfit/` — orbit-fitting work: Levenberg–Marquardt core (`lm.rs`), `Model` trait
  (`model.rs`), analytic and numerical fitters, the unit-vector-from-invariants solver, and
  `plans.md`, and `fitter.rs`, the earlier finite-difference LM fitter on the `nbody` engine
  (superseded by `spacerocks::orbfit`, the layup port).
- `heliostack/initial_condition.rs` — heliostack initial conditions.
- `mpc/detections.rs` — MPC detection parsing.
- `nbody_forces/` — commented-out drag and radiation-pressure forces for the n-body code.

To revive a file, move it back under `src/` and add it to the parent `mod.rs`.
