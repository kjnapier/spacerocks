//! N-body integrators on synthetic systems (no kernels needed).

use spacerocks::constants::GRAVITATIONAL_CONSTANT;
use spacerocks::nbody::{Integrator, Leapfrog, Simulation, Trace, WisdomHolman, IAS15};
use spacerocks::{SpaceRock, Time};

const T0: f64 = 2460000.5;

fn rock(name: &str, mass: f64, x: [f64; 3], v: [f64; 3]) -> SpaceRock {
    let mut r = SpaceRock::from_xyz(name, x[0], x[1], x[2], v[0], v[1], v[2], Time::new(T0, "tdb", "jd").unwrap(), "J2000", "SSB").unwrap();
    r.set_mass(mass);
    r
}

/// Sun, Jupiter and Saturn on roughly circular, mutually inclined orbits, plus a test particle.
fn outer_system(integrator: Box<dyn Integrator + Send + Sync>, test_particle: bool) -> Simulation {
    let mut sim = Simulation::new(&Time::new(T0, "tdb", "jd").unwrap(), "J2000", "SSB").unwrap();
    sim.integrator = integrator;
    let vc = |a: f64| (GRAVITATIONAL_CONSTANT / a).sqrt();
    sim.add(rock("sun", 1.0, [0.0; 3], [0.0; 3])).unwrap();
    sim.add(rock("jupiter", 9.547919e-4, [5.2, 0.0, 0.0], [0.0, vc(5.2) * 1.02, 2e-4])).unwrap();
    sim.add(rock("saturn", 2.858860e-4, [0.0, 9.55, 0.3], [-vc(9.55) * 0.97, 0.0, 0.0])).unwrap();
    if test_particle {
        sim.add(rock("tp", 0.0, [-2.7, 0.3, 0.1], [1e-3, -vc(2.7) * 1.1, 5e-4])).unwrap();
    }
    sim
}

fn max_offset(a: &Simulation, b: &Simulation) -> (f64, f64) {
    a.particles.iter().zip(&b.particles).fold((0.0, 0.0), |(dx, dv), (p, q)| {
        (f64::max(dx, (p.position - q.position).norm()), f64::max(dv, (p.velocity - q.velocity).norm()))
    })
}

#[test]
fn wisdom_holman_two_body_is_exact() {
    // A test particle around the Sun alone: the Kepler drift is the whole step, so even a
    // seventh of an orbit per step returns the particle to its start after 100 periods.
    let mut sim = Simulation::new(&Time::new(T0, "tdb", "jd").unwrap(), "J2000", "SSB").unwrap();
    let a: f64 = 1.3;
    let e: f64 = 0.4;
    let mu = GRAVITATIONAL_CONSTANT;
    let period = std::f64::consts::TAU * (a.powi(3) / mu).sqrt();
    let rp = a * (1.0 - e);
    let vp = (mu * (1.0 + e) / rp).sqrt();
    sim.add(rock("sun", 1.0, [0.0; 3], [0.0; 3])).unwrap();
    sim.add(rock("tp", 0.0, [rp, 0.0, 0.0], [0.0, vp, 0.0])).unwrap();
    sim.integrator = Box::new(WisdomHolman::new(period / 7.0));
    let start = sim.clone();
    for _ in 0..700 {
        sim.step();
    }
    let (dx, dv) = max_offset(&sim, &start);
    assert!(dx < 1e-10 && dv < 1e-12, "dx = {dx:e} AU, dv = {dv:e} AU/day");
    assert!((sim.epoch.tdb().jd() - (T0 + 100.0 * period)).abs() < 1e-6);
}

#[test]
fn wisdom_holman_is_time_reversible() {
    let mut sim = outer_system(Box::new(WisdomHolman::new(20.0)), true);
    let start = sim.clone();
    for _ in 0..2000 {
        sim.step();
    }
    sim.integrator.set_timestep(-20.0);
    for _ in 0..2000 {
        sim.step();
    }
    let (dx, dv) = max_offset(&sim, &start);
    assert!(dx < 1e-10 && dv < 1e-12, "dx = {dx:e} AU, dv = {dv:e} AU/day");
}

#[test]
fn wisdom_holman_energy_is_bounded() {
    // ~1100 years with a 100-day step (about 43 steps per Jupiter orbit).
    let dt = 100.0;
    let mut wh = outer_system(Box::new(WisdomHolman::new(dt)), false);
    let mut lf = outer_system(Box::new(Leapfrog::new(dt)), false);
    let e0 = wh.energy();
    let (mut wh_early, mut wh_late, mut lf_max) = (0.0f64, 0.0f64, 0.0f64);
    let n = 4000;
    for i in 0..n {
        wh.step();
        lf.step();
        let err = ((wh.energy() - e0) / e0).abs();
        if i < n / 2 {
            wh_early = wh_early.max(err);
        } else {
            wh_late = wh_late.max(err);
        }
        lf_max = lf_max.max(((lf.energy() - e0) / e0).abs());
    }
    // Bounded: the second half is no worse than the first, and far better than leapfrog.
    let wh_max = wh_early.max(wh_late);
    assert!(wh_max < 1e-5, "Wisdom-Holman max |dE/E| = {wh_max:e}");
    assert!(wh_late <= 2.0 * wh_early, "energy error drifting: {wh_early:e} then {wh_late:e}");
    assert!(wh_max < 1e-2 * lf_max, "WH {wh_max:e} vs leapfrog {lf_max:e}");
}

#[test]
fn wisdom_holman_matches_ias15() {
    // 50 years with a 5-day step against IAS15.
    let target = Time::new(T0 + 50.0 * 365.25, "tdb", "jd").unwrap();
    let mut wh = outer_system(Box::new(WisdomHolman::new(5.0)), true);
    let mut ias = outer_system(Box::new(IAS15::new(1.0)), true);
    wh.integrate(&target);
    ias.integrate(&target);
    assert!((wh.epoch.tdb().jd() - target.tdb().jd()).abs() < 1e-9);
    let (dx, _) = max_offset(&wh, &ias);
    assert!(dx < 1e-5, "max position difference = {dx:e} AU");
}


/// Sun, Jupiter, and a test particle 0.08 AU from Jupiter (inside its Hill sphere).
fn encounter_system(integrator: Box<dyn Integrator + Send + Sync>) -> Simulation {
    let mut sim = Simulation::new(&Time::new(T0, "tdb", "jd").unwrap(), "J2000", "SSB").unwrap();
    sim.integrator = integrator;
    let vj = (GRAVITATIONAL_CONSTANT / 5.2).sqrt();
    sim.add(rock("sun", 1.0, [0.0; 3], [0.0; 3])).unwrap();
    sim.add(rock("jupiter", 9.547919e-4, [5.2, 0.0, 0.0], [0.0, vj, 0.0])).unwrap();
    sim.add(rock("tp", 0.0, [5.2 - 0.3, -0.25, 0.01], [0.0, vj + 1.2e-3, 0.0])).unwrap();
    sim
}

/// A sungrazing comet (q = 0.05 AU, e = 0.98) and Jupiter.
fn comet_system(integrator: Box<dyn Integrator + Send + Sync>) -> Simulation {
    let mut sim = Simulation::new(&Time::new(T0, "tdb", "jd").unwrap(), "J2000", "SSB").unwrap();
    sim.integrator = integrator;
    let mu = GRAVITATIONAL_CONSTANT;
    let vj = (mu / 5.2).sqrt();
    let (q, e) = (0.05, 0.98);
    let a = q / (1.0 - e);
    let ra = a * (1.0 + e);
    let va = (mu * (1.0 - e) / ra).sqrt();
    sim.add(rock("sun", 1.0, [0.0; 3], [0.0; 3])).unwrap();
    sim.add(rock("jupiter", 9.547919e-4, [0.0, 5.2, 0.0], [-vj, 0.0, 0.0])).unwrap();
    sim.add(rock("comet", 0.0, [-ra, 0.0, 0.0], [0.0, -va, 0.0])).unwrap();
    sim
}

/// Position error of particle `idx` after integrating to `days` with each integrator, against IAS15.
fn errors_vs_ias15(system: fn(Box<dyn Integrator + Send + Sync>) -> Simulation, dt: f64, days: f64, idx: usize) -> (f64, f64) {
    let target = Time::new(T0 + days, "tdb", "jd").unwrap();
    let mut ias = system(Box::new(IAS15::new(1.0)));
    let mut wh = system(Box::new(WisdomHolman::new(dt)));
    let mut tr = system(Box::new(Trace::new(dt)));
    for sim in [&mut ias, &mut wh, &mut tr] {
        sim.integrate(&target);
    }
    let err = |s: &Simulation| (s.particles[idx].position - ias.particles[idx].position).norm();
    (err(&wh), err(&tr))
}

#[test]
fn trace_is_wisdom_holman_without_encounters() {
    let mut wh = outer_system(Box::new(WisdomHolman::new(20.0)), true);
    let mut tr = outer_system(Box::new(Trace::new(20.0)), true);
    for _ in 0..2000 {
        wh.step();
        tr.step();
    }
    let (dx, dv) = max_offset(&wh, &tr);
    assert!(dx < 1e-14 && dv < 1e-16, "dx = {dx:e} AU, dv = {dv:e} AU/day");
}

#[test]
fn trace_handles_close_encounters() {
    // A test particle passing deep inside Jupiter's Hill sphere.
    let (wh, tr) = errors_vs_ias15(encounter_system, 20.0, 3000.0, 2);
    assert!(tr < 1e-4 && tr < 1e-2 * wh, "WH {wh:e} AU, TRACE {tr:e} AU");
}

#[test]
fn trace_handles_pericenter_passages() {
    // Five passages of a sungrazer at q = 0.05 AU.
    let (wh, tr) = errors_vs_ias15(comet_system, 5.0, 8000.0, 2);
    assert!(tr < 0.05 && tr < 0.1 * wh, "WH {wh:e} AU, TRACE {tr:e} AU");
}

#[test]
fn trace_is_time_reversible() {
    for system in [encounter_system as fn(Box<dyn Integrator + Send + Sync>) -> Simulation, comet_system] {
        let mut sim = system(Box::new(Trace::new(20.0)));
        let start = sim.clone();
        for _ in 0..200 {
            sim.step();
        }
        sim.integrator.set_timestep(-20.0);
        for _ in 0..200 {
            sim.step();
        }
        let (dx, dv) = max_offset(&sim, &start);
        assert!(dx < 1e-8 && dv < 1e-10, "dx = {dx:e} AU, dv = {dv:e} AU/day");
    }
}
