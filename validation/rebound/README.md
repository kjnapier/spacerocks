# spacerocks vs REBOUND: Wisdom–Holman and TRACE

`compare.py` integrates the same initial conditions with spacerocks (`WisdomHolman`, `Trace`,
`IAS15`, through `examples/nbody_compare.rs`) and REBOUND (WHFast in democratic heliocentric
coordinates, TRACE, IAS15), and writes final-position errors against REBOUND IAS15, energy
errors and per-step wall times to `results.json`.

    pip install rebound numpy
    cargo build --release --example nbody_compare
    python3 validation/rebound/compare.py > validation/rebound/results.json

`results.json` here is from REBOUND 5.2.1 on a 4-core cloud container.
