//! Per-step cost of IAS15 and of its compensated summation (ASSIST's benchmark: N copies of
//! Holman, 10 years).
//!
//!     SPACEROCKS_KERNELS=/path/to/kernels cargo run --release --example bench_ias15

use std::path::PathBuf;
use std::time::Instant;

use spacerocks::assist::{SpiceSimulation, Summation};
use spacerocks::{SpaceRock, SpiceKernel, Time};

fn main() {
    let dir = PathBuf::from(std::env::var("SPACEROCKS_KERNELS").expect("set SPACEROCKS_KERNELS"));
    let mut k = SpiceKernel::new();
    k.load(dir.join("de440s.bsp")).or_else(|_| k.load(dir.join("de440.bsp"))).unwrap();
    k.load(dir.join("sb441-n16.bsp")).unwrap();
    let t0v = 2451545.0 + 8416.5;
    let t0 = Time::new(t0v, "tdb", "jd").unwrap();
    let x0 = [-2.724183384883979E+00, -3.523994546329214E-02, 9.036596202793466E-02,
              -1.374545432301129E-04, -1.027075301472321E-02, -4.195690627695180E-03];
    let mut modes: Vec<(&str, Summation)> = vec![("none", Summation::None), ("kahan", Summation::Kahan), ("full", Summation::Full)];
    if let Ok(m) = std::env::var("BENCH_SUM") {
        modes.retain(|(name, _)| *name == m);
    }
    // BENCH_N=100 BENCH_MODE=newtonian restricts the run (for profiling).
    let ns: Vec<usize> = std::env::var("BENCH_N").map(|v| vec![v.parse().unwrap()]).unwrap_or(vec![1, 100, 1000]);
    let only = std::env::var("BENCH_MODE").ok();
    for n in ns {
        for newtonian in [false, true] {
            if let Some(o) = &only { if (o == "newtonian") != newtonian { continue; } }
            let mut line = format!("N={:5} {:9}", n, if newtonian { "newtonian" } else { "full" });
            for (name, m) in &modes {
                let mut best = f64::INFINITY;
                let mut steps = 0;
                for _ in 0..(if n < 1000 { 5 } else { 3 }) {
                    let mut sim = SpiceSimulation::horizons(&t0, &k).unwrap();
                    if newtonian { sim.newtonian_only(); }
                    sim.integrator.set_summation(*m);
                    for i in 0..n {
                        let r = SpaceRock::from_xyz(&format!("p{i}"), x0[0] + i as f64 * 1e-10, x0[1], x0[2], x0[3], x0[4], x0[5], t0.clone(), "J2000", "SSB").unwrap();
                        sim.add(r).unwrap();
                    }
                    let t = Instant::now();
                    steps = 0;
                    // Step manually (epochs are days since jd_ref) so the count is known.
                    loop {
                        if sim.state.particles_1[0].epoch >= 3652.5 { break; }
                        sim.step(&k).unwrap();
                        steps += 1;
                    }
                    best = best.min(t.elapsed().as_secs_f64());
                }
                line += &format!("  {name} {:8.2} ms ({steps} steps, {:.2} us/step/particle)", best * 1e3, best / steps as f64 / n as f64 * 1e6);
            }
            println!("{line}");
        }
    }
}
