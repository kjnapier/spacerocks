//! Force models: analytic partial derivatives against finite differences, and the variational
//! equations against finite-difference propagation. Need `de440s.bsp` and `sb441-n16.bsp` in
//! `SPACEROCKS_KERNELS`; skipped otherwise.

use std::path::PathBuf;

use nalgebra::Vector3;
use spacerocks::assist::forces::{
    assist_default_forces, EarthHarmonics, Force, GrEih, GrPotential, GrSimple, NewtonianGravity, NonGravitational, SolarJ2,
};
use spacerocks::assist::{EphemerisConstants, SimulationState, SpiceSimulation};
use spacerocks::{SpaceRock, SpiceKernel, Time};

const T0: f64 = 2460300.5;

fn kernel() -> Option<SpiceKernel> {
    let dir = PathBuf::from(std::env::var("SPACEROCKS_KERNELS").ok()?);
    let mut k = SpiceKernel::new();
    k.load(dir.join("de440s.bsp")).ok()?;
    k.load(dir.join("sb441-n16.bsp")).ok()?;
    Some(k)
}

fn tdb(jd: f64) -> Time {
    Time::new(jd, "tdb", "jd").unwrap()
}

/// A simulation state holding one particle near the Earth (so every force matters).
fn state_near_earth(k: &SpiceKernel, offset: Vector3<f64>, vel: Vector3<f64>, nongrav: [f64; 3]) -> SimulationState {
    let e = k.state_au(399, 0, T0).unwrap();
    let mut rock = SpaceRock::from_xyz("p", e[0] + offset.x, e[1] + offset.y, e[2] + offset.z, e[3] + vel.x, e[4] + vel.y, e[5] + vel.z, tdb(T0), "J2000", "SSB").unwrap();
    rock.set_nongrav(nongrav[0], nongrav[1], nongrav[2]);
    let mut sim = SpiceSimulation::horizons(&tdb(T0), k).unwrap();
    sim.add(rock).unwrap();
    sim.state
}

fn accel(force: &dyn Force, state: &mut SimulationState) -> Vector3<f64> {
    force.calculate_acceleration(state)[0]
}

/// Compare a force's analytic Jacobian with central differences.
fn check_jacobian(name: &str, force: &dyn Force, state: &SimulationState) {
    let mut s = state.clone();
    for p in s.particles_1.iter_mut() {
        p.acceleration = Vector3::zeros();
        p.stm = [[0.0; 6]; 6];
        p.nongrav_partials = [[0.0; 3]; 3];
    }
    force.apply_acceleration_and_stm(&mut s);
    let a0 = s.particles_1[0].acceleration;
    let jac = s.particles_1[0].stm;
    let dadk = s.particles_1[0].nongrav_partials;
    assert!(a0.norm() > 0.0, "{}: no acceleration", name);
    assert_eq!(a0, accel(force, &mut state.clone()), "{}: apply_acceleration and apply_acceleration_and_stm differ", name);

    let scale = [1e-7, 1e-7, 1e-7, 1e-9, 1e-9, 1e-9];
    let mut worst = 0.0f64;
    let mut norm = 0.0f64;
    for j in 0..6 {
        let h = scale[j];
        let mut plus = state.clone();
        let mut minus = state.clone();
        if j < 3 {
            plus.particles_1[0].position[j] += h;
            minus.particles_1[0].position[j] -= h;
        } else {
            plus.particles_1[0].velocity[j - 3] += h;
            minus.particles_1[0].velocity[j - 3] -= h;
        }
        let fd = (accel(force, &mut plus) - accel(force, &mut minus)) / (2.0 * h);
        for i in 0..3 {
            worst = worst.max((fd[i] - jac[3 + i][j]).abs());
            norm = norm.max(fd[i].abs()).max(jac[3 + i][j].abs());
        }
    }
    let rel_r_v = worst / norm.max(1e-300);
    let jac_norm = norm;
    let (mut worst, mut norm) = (0.0f64, 0.0f64);
    for kk in 0..3 {
        let h = 1e-3 * state.particles_1[0].nongrav[kk].abs().max(1e-12);
        let mut plus = state.clone();
        let mut minus = state.clone();
        plus.particles_1[0].nongrav[kk] += h;
        minus.particles_1[0].nongrav[kk] -= h;
        let fd = (accel(force, &mut plus) - accel(force, &mut minus)) / (2.0 * h);
        for i in 0..3 {
            worst = worst.max((fd[i] - dadk[i][kk]).abs());
            norm = norm.max(fd[i].abs());
        }
    }
    let rel = rel_r_v.max(if norm > 0.0 { worst / norm } else { 0.0 });
    eprintln!("{:<18} |da/dx| {:.2e}, max relative error {:.1e}", name, jac_norm, rel);
    assert!(rel < 1e-6, "{}: Jacobian disagrees with finite differences (rel {:e})", name, rel);
}

#[test]
fn force_jacobians_match_finite_differences() {
    let Some(k) = kernel() else {
        eprintln!("kernels not available; skipping");
        return;
    };
    let c = EphemerisConstants::from_kernel(&k);
    assert_eq!(c, EphemerisConstants::DE440, "constants read from de440s.bsp");
    // 0.002 AU (~300,000 km, ~47 Earth radii) from the Earth, off the equator.
    let st = state_near_earth(&k, Vector3::new(0.0012, -0.0009, 0.0011), Vector3::new(0.004, 0.002, -0.003), [3e-9, -2e-10, 1e-10]);
    let forces: Vec<(&str, Box<dyn Force + Send + Sync>)> = vec![
        ("newtonian", Box::new(NewtonianGravity)),
        ("earth J2-J4", Box::new(EarthHarmonics::new(&c))),
        ("solar J2", Box::new(SolarJ2::new(&c))),
        ("GR simple", Box::new(GrSimple::new(&c))),
        ("GR potential", Box::new(GrPotential::new(&c))),
        ("GR EIH (Sun)", Box::new(GrEih::new(&c))),
        ("GR EIH (all 11)", Box::new(GrEih::new(&c).with_sources(11))),
        ("nongrav", Box::new(NonGravitational::default())),
        ("nongrav (comet)", Box::new(NonGravitational::comet())),
    ];
    for (name, f) in &forces {
        check_jacobian(name, f.as_ref(), &st);
    }
}

/// Integrate the state transition matrix (and d/dA) with variational particles and compare it
/// with finite differences of whole propagations.
#[test]
fn variational_particles_match_finite_differences() {
    let Some(k) = kernel() else {
        eprintln!("kernels not available; skipping");
        return;
    };
    let c = EphemerisConstants::from_kernel(&k);
    // An Apophis-like orbit with non-gravitational parameters, over 60 days.
    let mut base = SpaceRock::from_kepler("x", 0.75, 0.19, 0.06, 2.2, 3.5, 1.0, tdb(T0), "J2000", "SSB").unwrap();
    base.set_nongrav(2e-10, -5e-12, 1e-12);
    let t1 = T0 + 60.0;

    let run = |rock: &SpaceRock, variations: bool| {
        let mut sim = SpiceSimulation::horizons(&tdb(T0), &k).unwrap();
        sim.set_forces(assist_default_forces(&c));
        sim.add(rock.clone()).unwrap();
        if variations {
            sim.add_full_variation("x").unwrap();
            for a in ["A1", "A2", "A3"] {
                sim.add_variation(a, "x").unwrap();
            }
        }
        sim.integrate_jd(t1, &k).unwrap();
        let p = &sim.state.particles[0];
        let state = [p.position.x, p.position.y, p.position.z, p.velocity.x, p.velocity.y, p.velocity.z];
        let vars: Vec<[f64; 6]> = sim
            .state
            .variational_particles
            .iter()
            .map(|v| [v.position.x, v.position.y, v.position.z, v.velocity.x, v.velocity.y, v.velocity.z])
            .collect();
        (state, vars)
    };
    let (_, vars) = run(&base, true);
    assert_eq!(vars.len(), 9);

    // Steps large enough that the integrator's own noise (different step sequences) is
    // negligible, small enough that the response is linear.
    let steps = [1e-5, 1e-5, 1e-5, 1e-7, 1e-7, 1e-7];
    let mut worst = 0.0f64;
    for j in 0..9 {
        let (plus, minus) = if j < 6 {
            let h = steps[j];
            let mut p = base.clone();
            let mut m = base.clone();
            if j < 3 {
                p.position[j] += h;
                m.position[j] -= h;
            } else {
                p.velocity[j - 3] += h;
                m.velocity[j - 3] -= h;
            }
            ((p, h), (m, h))
        } else {
            let h = 1e-9;
            let mut ng_p = base.nongrav().unwrap();
            let mut ng_m = ng_p;
            ng_p[j - 6] += h;
            ng_m[j - 6] -= h;
            let mut p = base.clone();
            let mut m = base.clone();
            p.set_nongrav(ng_p[0], ng_p[1], ng_p[2]);
            m.set_nongrav(ng_m[0], ng_m[1], ng_m[2]);
            ((p, h), (m, h))
        };
        let (sp, _) = run(&plus.0, false);
        let (sm, _) = run(&minus.0, false);
        let mut wj = 0.0f64;
        for i in 0..6 {
            let fd = (sp[i] - sm[i]) / (2.0 * plus.1);
            let rel = (fd - vars[j][i]).abs() / vars[j].iter().fold(0.0f64, |a, b| a.max(b.abs()));
            wj = wj.max(rel);
        }
        eprintln!("variation {}: max relative difference {:.2e}", j, wj);
        worst = worst.max(wj);
    }
    eprintln!("variational vs finite-difference propagation: max relative difference {:.2e}", worst);
    assert!(worst < 2e-5, "{:e}", worst);
}
