# spacerocks vs REBOUND: Wisdom–Holman and TRACE

`compare.py` integrates the same initial conditions with spacerocks (`WisdomHolman`, `Trace`,
`IAS15`, through `examples/nbody_compare.rs`) and REBOUND (WHFast in democratic heliocentric
coordinates, TRACE, IAS15), and writes final-position errors against REBOUND IAS15, energy
errors and per-step wall times to `results.json`.

    pip install rebound numpy
    cargo build --release --example nbody_compare
    python3 validation/rebound/compare.py > validation/rebound/results.json

`results.json` here is from REBOUND 5.2.1 on a 4-core cloud container.

The spacerocks timings go through `Simulation::steps`, so `WisdomHolman` and `Trace` split test
particles between threads (rayon). Set `RAYON_NUM_THREADS=1` for single-threaded numbers: on
one thread the blocked steps are as fast as or a little faster than stepping one at a time.
REBOUND runs on one thread.
