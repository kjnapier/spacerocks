//! Numerical propagation tests. Need `de440s.bsp` and `sb441-n16.bsp` in the directory named by
//! `SPACEROCKS_KERNELS`; skipped otherwise.

use spacerocks::{SpaceRock, SpiceKernel, Time};
use std::path::PathBuf;

fn kernel() -> Option<SpiceKernel> {
    let dir = PathBuf::from(std::env::var("SPACEROCKS_KERNELS").ok()?);
    let mut k = SpiceKernel::new();
    k.load(dir.join("de440s.bsp")).ok()?;
    k.load(dir.join("sb441-n16.bsp")).ok()?;
    Some(k)
}

#[test]
fn propagate_is_independent_of_input_frame_and_origin() {
    let Some(k) = kernel() else {
        eprintln!("kernels not available; skipping");
        return;
    };
    let t0 = Time::new(2460762.549988426, "TDB", "JD").unwrap();
    let t1 = Time::new(2460762.549988426 + 300.0, "TDB", "JD").unwrap();
    let a = SpaceRock::from_xyz("x", 2.963305899720348, -1.627586306680811, -0.7799786968810375,
        0.004951381894813546, 0.006677060249157604, 0.002540378471598749, t0, "J2000", "SSB").unwrap();

    let mut b = a.clone();
    b.to_helio(&k).unwrap();
    b.change_reference_plane("ECLIPJ2000").unwrap();

    let mut a1 = a.clone();
    a1.propagate(&t1, &k).unwrap();
    b.propagate(&t1, &k).unwrap();

    // Output keeps the input's frame and origin.
    assert_eq!(b.reference_plane, spacerocks::ReferencePlane::ECLIPJ2000);
    assert_eq!(b.origin, spacerocks::Origin::SUN);

    b.to_ssb(&k).unwrap();
    b.change_reference_plane("J2000").unwrap();
    assert!((a1.position - b.position).norm() < 1e-10);
    assert!((a1.velocity - b.velocity).norm() < 1e-12);
}

#[test]
fn ias15_rejections_with_variational_particles_and_coverage_errors() {
    use spacerocks::assist::SpiceSimulation;
    let Some(k) = kernel() else {
        eprintln!("kernels not available; skipping");
        return;
    };
    let t0 = Time::new(2460762.549988426, "TDB", "JD").unwrap();
    let t1 = Time::new(2460762.549988426 + 200.0, "TDB", "JD").unwrap();
    let rock = SpaceRock::from_xyz("holman", 2.963305899720348, -1.627586306680811, -0.7799786968810375,
        0.004951381894813546, 0.006677060249157604, 0.002540378471598749, t0.clone(), "J2000", "SSB").unwrap();

    // A deliberately huge first step forces rejections while variational particles exist.
    let mut sim = SpiceSimulation::horizons(&t0, &k).unwrap();
    sim.add(rock.clone()).unwrap();
    sim.add_full_variation("holman").unwrap();
    sim.integrator.set_timestep(80.0);
    sim.integrate(&t1, &k).unwrap();

    let mut plain = SpiceSimulation::horizons(&t0, &k).unwrap();
    plain.add(rock.clone()).unwrap();
    plain.integrate(&t1, &k).unwrap();
    assert!((sim.state.particles[0].position - plain.state.particles[0].position).norm() < 1e-10);

    // Integrating past the ephemeris coverage is an error, not a panic.
    let mut r = rock.clone();
    let far = Time::new(2460762.549988426 + 1.0e6, "TDB", "JD").unwrap();
    assert!(r.propagate(&far, &k).is_err());
}
