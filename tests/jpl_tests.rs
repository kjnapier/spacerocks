//! Integrations against JPL Horizons (ASSIST's `holman_spk` unit test). Needs `de440s.bsp` and
//! `sb441-n16.bsp` in `SPACEROCKS_KERNELS`; skipped otherwise.

use std::path::PathBuf;

use spacerocks::assist::SpiceSimulation;
use spacerocks::{SpaceRock, SpiceKernel, Time};

fn kernel() -> Option<SpiceKernel> {
    let dir = PathBuf::from(std::env::var("SPACEROCKS_KERNELS").ok()?);
    let mut k = SpiceKernel::new();
    k.load(dir.join("de440s.bsp")).ok()?;
    k.load(dir.join("sb441-n16.bsp")).ok()?;
    Some(k)
}

#[test]
fn holman_30_days_matches_horizons_to_10_cm() {
    let Some(k) = kernel() else {
        eprintln!("kernels not available; skipping");
        return;
    };
    let t0 = Time::new(2451545.0 + 8416.5, "tdb", "jd").unwrap();
    let rock = SpaceRock::from_xyz("holman", -2.724183384883979E+00, -3.523994546329214E-02, 9.036596202793466E-02,
        -1.374545432301129E-04, -1.027075301472321E-02, -4.195690627695180E-03, t0.clone(), "J2000", "SSB").unwrap();
    let mut sim = SpiceSimulation::horizons(&t0, &k).unwrap();
    sim.add(rock).unwrap();
    sim.integrate_jd(2451545.0 + 8446.5, &k).unwrap();
    let p = &sim.state.particles[0];
    // JPL Horizons, DES=2003666, 2023-02-16 TDB, barycentric ICRF.
    let jpl = [-2.710320457933958E+00, -3.424507930535848E-01, -3.582442972611413E-02,
               1.059255302926290E-03, -1.018748422976772E-02, -4.207712906489264E-03];
    let au_m = 149597870700.0;
    for i in 0..3 {
        let dx = (p.position[i] - jpl[i]).abs() * au_m;
        let dv = (p.velocity[i] - jpl[3 + i]).abs() * au_m;
        assert!(dx < 0.1, "position component {} off by {} m", i, dx);
        assert!(dv < 5e-3, "velocity component {} off by {} m/day", i, dv);
    }
}

/// Round trips in the round-off-limited regime (tight epsilon, many steps): compensated
/// summation of positions and velocities should bring the clones back closer than plain sums.
#[test]
fn compensated_summation_reduces_round_trip_error() {
    use spacerocks::assist::{AdaptiveMode, Summation};
    let Some(k) = kernel() else {
        eprintln!("kernels not available; skipping");
        return;
    };
    let t0v = 2451545.0 + 8416.5;
    let t0 = Time::new(t0v, "tdb", "jd").unwrap();
    let x0 = [-2.724183384883979E+00, -3.523994546329214E-02, 9.036596202793466E-02,
              -1.374545432301129E-04, -1.027075301472321E-02, -4.195690627695180E-03];
    let n = 32;
    let median_error = |summation: Summation| -> f64 {
        let mut sim = SpiceSimulation::horizons(&t0, &k).unwrap();
        sim.newtonian_only();
        sim.integrator.set_adaptive_mode(AdaptiveMode::Global);
        sim.integrator.set_epsilon(1e-11);
        sim.integrator.set_summation(summation);
        for i in 0..n {
            let d = (i as f64 - n as f64 / 2.0) * 1e-9;
            sim.add(SpaceRock::from_xyz(&format!("c{i}"), x0[0] + d, x0[1] - d, x0[2], x0[3], x0[4], x0[5], t0.clone(), "J2000", "SSB").unwrap()).unwrap();
        }
        sim.integrate_jd(t0v + 1000.0, &k).unwrap();
        sim.integrate_jd(t0v, &k).unwrap();
        let mut e: Vec<f64> = (0..n)
            .map(|i| {
                let d = (i as f64 - n as f64 / 2.0) * 1e-9;
                let p = &sim.state.particles[i];
                ((p.position[0] - x0[0] - d).powi(2) + (p.position[1] - x0[1] + d).powi(2) + (p.position[2] - x0[2]).powi(2)).sqrt()
            })
            .collect();
        e.sort_by(|a, b| a.partial_cmp(b).unwrap());
        e[n / 2] * 149597870700.0
    };
    let plain = median_error(Summation::None);
    let kahan = median_error(Summation::Kahan);
    eprintln!("median round-trip error: plain {plain:.3e} m, kahan {kahan:.3e} m");
    assert!(kahan < plain, "compensated {kahan} m vs plain {plain} m");
}
