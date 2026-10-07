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
    a.particles().iter().zip(&b.particles()).fold((0.0, 0.0), |(dx, dv), (p, q)| {
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
    let err = |s: &Simulation| (s.particle(idx).position - ias.particle(idx).position).norm();
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

#[test]
fn ias15_does_not_depend_on_output_cadence() {
    // Frequent outputs force short final steps, which used to make rejected steps apply twice.
    let mut once = comet_system(Box::new(IAS15::new(1.0)));
    let mut often = comet_system(Box::new(IAS15::new(1.0)));
    once.integrate(&Time::new(T0 + 8000.0, "tdb", "jd").unwrap());
    for k in 1..=1600 {
        often.integrate(&Time::new(T0 + 5.0 * k as f64, "tdb", "jd").unwrap());
    }
    let (dx, dv) = max_offset(&once, &often);
    assert!(dx < 1e-7 && dv < 1e-9, "dx = {dx:e} AU, dv = {dv:e} AU/day");
}

#[test]
fn integrate_or_interpolate_matches_integrate() {
    // Interpolated states agree with exact integration, and the integrator's own path doesn't
    // depend on how many epochs are requested.
    let mut exact = comet_system(Box::new(IAS15::new(1.0)));
    let mut interp = comet_system(Box::new(IAS15::new(1.0)));
    let mut once = comet_system(Box::new(IAS15::new(1.0)));
    let mut worst: f64 = 0.0;
    for k in 1..=1600 {
        let epoch = Time::new(T0 + 5.0 * k as f64 + 0.37, "tdb", "jd").unwrap();
        exact.integrate(&epoch);
        interp.integrate_or_interpolate(&epoch);
        assert!((interp.epoch.jd() - epoch.jd()).abs() < 1e-9);
        worst = worst.max(max_offset(&exact, &interp).0);
    }
    assert!(worst < 1e-7, "worst offset {worst:e} AU");

    // Going back inside the last step, then on to the end, lands where a single call does.
    let end = Time::new(T0 + 8100.0, "tdb", "jd").unwrap();
    interp.integrate_or_interpolate(&Time::new(T0 + 8000.0, "tdb", "jd").unwrap());
    interp.integrate_or_interpolate(&end);
    once.integrate_or_interpolate(&end);
    let (dx, dv) = max_offset(&once, &interp);
    assert!(dx == 0.0 && dv == 0.0, "dx = {dx:e} AU, dv = {dv:e} AU/day");
}

#[test]
fn integrate_or_interpolate_runs_backwards() {
    let mut fwd = outer_system(Box::new(IAS15::new(1.0)), true);
    let start = fwd.clone();
    fwd.integrate_or_interpolate(&Time::new(T0 + 3000.0, "tdb", "jd").unwrap());
    fwd.integrate_or_interpolate(&Time::new(T0, "tdb", "jd").unwrap());
    let (dx, dv) = max_offset(&fwd, &start);
    assert!(dx < 1e-9 && dv < 1e-11, "dx = {dx:e} AU, dv = {dv:e} AU/day");
}

#[test]
fn integrate_or_interpolate_falls_back_without_dense_output() {
    let mut a = outer_system(Box::new(WisdomHolman::new(20.0)), true);
    let mut b = outer_system(Box::new(WisdomHolman::new(20.0)), true);
    let epoch = Time::new(T0 + 1234.5, "tdb", "jd").unwrap();
    a.integrate(&epoch);
    b.integrate_or_interpolate(&epoch);
    let (dx, dv) = max_offset(&a, &b);
    assert!(dx == 0.0 && dv == 0.0);
}

/// Sun, Jupiter and Saturn with a population of test particles on quiet orbits, and with
/// `encounters` also some passing close to Jupiter and a sungrazer, so TRACE flags both kinds
/// of encounter now and then.
fn population(integrator: Box<dyn Integrator + Send + Sync>, encounters: bool) -> Simulation {
    let mut sim = outer_system(integrator, false);
    let vc = |a: f64| (GRAVITATIONAL_CONSTANT / a).sqrt();
    let mut seed: u64 = 12345;
    let mut uniform = || {
        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (seed >> 11) as f64 / (1u64 << 53) as f64
    };
    for k in 0..48 {
        let a = 2.0 + 30.0 * uniform();
        let phase = std::f64::consts::TAU * uniform();
        let (s, c) = phase.sin_cos();
        let v = vc(a) * (0.9 + 0.2 * uniform());
        sim.add(rock(&format!("tp{k}"), 0.0, [a * c, a * s, 0.1 * a * (uniform() - 0.5)], [-v * s, v * c, 0.02 * v * (uniform() - 0.5)])).unwrap();
    }
    if !encounters {
        return sim;
    }
    // Near Jupiter's orbit, close behind it.
    let vj = vc(5.2);
    for k in 0..8 {
        let dy = -0.2 - 0.05 * k as f64;
        sim.add(rock(&format!("near{k}"), 0.0, [5.2 - 0.1, dy, 0.01], [0.0, vj * 1.02 + 2e-4 * k as f64, 0.0])).unwrap();
    }
    sim.add(rock("comet", 0.0, [-9.5, 0.5, 0.0], [0.0, -1.2e-3, 0.0])).unwrap();
    sim
}

fn assert_same(a: &Simulation, b: &Simulation) {
    assert_eq!(a.epoch.tdb().jd(), b.epoch.tdb().jd());
    for (p, q) in a.particles().iter().zip(&b.particles()) {
        assert_eq!(p.name, q.name);
        assert_eq!((p.position, p.velocity), (q.position, q.velocity), "{}", p.name);
        assert_eq!(p.epoch.tdb().jd(), q.epoch.tdb().jd());
    }
}

#[test]
fn steps_match_step_by_step() {
    // `steps` takes test particles through blocks of steps (in parallel); it must do exactly
    // the arithmetic of single steps, forwards and backwards, with and without encounters.
    for dt in [20.0, -20.0, 3.0] {
        let integrators: [fn(f64) -> Box<dyn Integrator + Send + Sync>; 2] = [|dt| Box::new(WisdomHolman::new(dt)), |dt| Box::new(Trace::new(dt))];
        for (make, encounters) in integrators.iter().flat_map(|m| [(m, false), (m, true)]) {
            let mut one = population(make(dt), encounters);
            let mut many = one.clone();
            for _ in 0..800 {
                one.step();
            }
            many.steps(500);
            many.steps(300);
            assert_same(&one, &many);
        }
    }
}

#[test]
fn steps_match_step_by_step_with_other_forces() {
    let integrators: [Box<dyn Integrator + Send + Sync>; 2] = [Box::new(WisdomHolman::new(20.0)), Box::new(Trace::new(20.0))];
    for integrator in integrators {
        let mut one = population(integrator, true);
        one.add_force(Box::new(spacerocks::nbody::forces::SolarGR));
        let mut many = one.clone();
        for _ in 0..300 {
            one.step();
        }
        many.steps(300);
        assert_same(&one, &many);
    }
}

#[test]
fn integrate_matches_stepping() {
    // `integrate` hands a fixed-step integrator all its full steps at once.
    for dt in [20.0, -20.0] {
        let mut a = population(Box::new(Trace::new(dt)), true);
        let mut b = a.clone();
        let target = T0 + 7013.25 * dt.signum();
        a.integrate(&Time::new(target, "tdb", "jd").unwrap());
        while (target - b.epoch.tdb().jd()).abs() >= dt.abs() {
            b.step();
        }
        b.integrator.set_timestep(target - b.epoch.tdb().jd());
        b.step();
        b.integrator.set_timestep(dt);
        assert_same(&a, &b);
        assert_eq!(a.integrator.timestep(), dt);
    }
}
