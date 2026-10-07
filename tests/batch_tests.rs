//! Tests for `spacerocks::batch`, `SpaceRock::apparent` and `PerturberCache`. The kernel tests
//! need `de440s.bsp` and `sb441-n16.bsp` in the directory named by `SPACEROCKS_KERNELS` and are
//! skipped otherwise. Epochs are inside the coverage of the 16-asteroid kernel.

use std::path::PathBuf;

use nalgebra::Vector3;
use spacerocks::assist::{PerturberCache, SpiceSimulation};
use spacerocks::batch::{self, BatchOptions, Method};
use spacerocks::{Observer, Origin, Population, ReferencePlane, SpaceRock, SpiceKernel, Time};

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

/// A few main-belt objects and TNOs, in a mix of frames and origins.
fn rocks(epoch: f64) -> Vec<SpaceRock> {
    let mut v = Vec::new();
    for i in 0..24 {
        let (a, e) = if i % 2 == 0 { (2.3 + 0.04 * i as f64, 0.12) } else { (38.0 + 0.3 * i as f64, 0.15) };
        let (plane, origin) = match i % 4 {
            0 => ("ECLIPJ2000", "SUN"),
            1 => ("J2000", "SSB"),
            2 => ("INVARIABLE", "SUN"),
            _ => ("J2000", "SUN"),
        };
        let mut r = SpaceRock::from_kepler(&format!("r{}", i), a * (1.0 - e), e, 0.05 * i as f64, 0.3 * i as f64, 0.7 * i as f64, 0.25 * i as f64, tdb(epoch), plane, origin).unwrap();
        r.set_absolute_magnitude(10.0 + 0.1 * i as f64);
        if i % 5 == 0 {
            r.set_nongrav(1e-10, -3e-12, 5e-13);
        }
        v.push(r);
    }
    v
}

fn observer_at(jd: f64, k: &SpiceKernel) -> Observer {
    // Geocentric observer (no Earth orientation kernel needed).
    let s = k.state_au(399, 0, jd).unwrap();
    let sun = k.state_au(10, 0, jd).unwrap();
    Observer::from_xyz(
        Vector3::new(s[0], s[1], s[2]),
        Some(Vector3::new(s[3], s[4], s[5])),
        tdb(jd),
        ReferencePlane::J2000,
        Origin::SSB,
        None,
    )
    .with_sun_position(Vector3::new(sun[0], sun[1], sun[2]))
}

#[test]
fn apparent_matches_observe() {
    let t = tdb(T0);
    let mut r = SpaceRock::from_xyz("x", 2.1, -1.4, 0.3, 0.004, 0.009, -0.001, t.clone(), "J2000", "SSB").unwrap();
    r.set_absolute_magnitude(12.3);
    let o = Observer::from_xyz(Vector3::new(-0.2, 0.97, 0.01), Some(Vector3::new(-0.017, -0.003, 0.0)), t, ReferencePlane::J2000, Origin::SSB, None)
        .with_sun_position(Vector3::new(0.004, -0.002, 0.0001));
    let a = r.apparent(&o).unwrap();
    let ob = r.observe(&o).unwrap();
    assert_eq!(a.ra, ob.ra());
    assert_eq!(a.dec, ob.dec());
    assert_eq!(Some(a.ra_rate), ob.ra_rate());
    assert_eq!(Some(a.dec_rate), ob.dec_rate());
    assert_eq!(Some(a.range), ob.range());
    assert_eq!(Some(a.range_rate), ob.range_rate());
    assert_eq!(Some(a.mag), ob.mag());
    assert!(a.phase > 0.0 && a.phase < std::f64::consts::PI);
    assert!(a.elong > 0.0 && a.elong < std::f64::consts::PI);

    // Wrong epoch or plane is an error.
    let mut o2 = o.clone();
    o2.epoch = tdb(T0 + 1.0);
    assert!(r.apparent(&o2).is_err());
    let mut o3 = o.clone();
    o3.change_reference_plane("ECLIPJ2000").unwrap();
    assert!(r.apparent(&o3).is_err());
}

#[test]
fn batch_propagate_with_chunk_size_one_is_identical_to_propagate() {
    let Some(k) = kernel() else {
        eprintln!("kernels not available; skipping");
        return;
    };
    let t1 = tdb(T0 + 150.0);
    let mut batch_rocks = rocks(T0);
    let opts = BatchOptions { chunk_size: 1, ..Default::default() };
    batch::propagate_batch(&mut batch_rocks, &t1, &k, &opts).unwrap();
    for (r, b) in rocks(T0).into_iter().zip(&batch_rocks) {
        let mut r = r;
        r.propagate(&t1, &k).unwrap();
        assert_eq!(r.reference_plane, b.reference_plane);
        assert_eq!(r.origin, b.origin);
        assert_eq!(r.epoch, b.epoch);
        assert!((r.position - b.position).norm() < 1e-14, "{}: {:e}", r.name, (r.position - b.position).norm());
        assert!((r.velocity - b.velocity).norm() < 1e-16);
    }
}

#[test]
fn batch_propagate_matches_propagate_and_handles_mixed_epochs() {
    let Some(k) = kernel() else {
        eprintln!("kernels not available; skipping");
        return;
    };
    let t1 = tdb(T0 + 120.0);
    // Two starting epochs (two groups) and one rock already at the target epoch.
    let mut input = rocks(T0);
    input.extend(rocks(T0 - 40.0));
    let at_target = SpaceRock::from_xyz("already", 3.0, 0.5, 0.1, -0.002, 0.009, 0.0, t1.clone(), "ECLIPJ2000", "SUN").unwrap();
    input.push(at_target.clone());

    let mut out = input.clone();
    batch::propagate_batch(&mut out, &t1, &k, &BatchOptions::default()).unwrap();
    assert_eq!(out.last().unwrap(), &at_target);

    for (r, b) in input.iter().zip(&out) {
        let mut r = r.clone();
        r.propagate(&t1, &k).unwrap();
        assert_eq!(r.reference_plane, b.reference_plane);
        assert_eq!(r.origin, b.origin);
        // Different step sequences; agreement at the integrator's tolerance.
        assert!((r.position - b.position).norm() < 1e-9, "{}: {:e}", r.name, (r.position - b.position).norm());
    }

    // Two-body batch equals analytic_propagate.
    let mut two = input.clone();
    batch::propagate_batch(&mut two, &t1, &k, &BatchOptions { method: Method::TwoBody, ..Default::default() }).unwrap();
    for (r, b) in input.iter().zip(&two) {
        let mut r = r.clone();
        if r.epoch != t1 {
            r.analytic_propagate(&t1).unwrap();
        }
        assert_eq!(&r, b);
    }
}

#[test]
fn ephemeris_matches_propagate_and_observe() {
    let Some(k) = kernel() else {
        eprintln!("kernels not available; skipping");
        return;
    };
    let input = rocks(T0);
    // Unsorted, before and after the rocks' epoch, with a duplicate and the epoch itself.
    let jds = [T0 + 30.0, T0 - 25.0, T0, T0 + 30.0, T0 + 3.5, T0 - 60.0];
    let observers: Vec<Observer> = jds.iter().map(|&t| observer_at(t, &k)).collect();
    let opts = BatchOptions { with_states: true, ..Default::default() };
    let eph = batch::ephemeris_rocks(&input, &observers, &k, &opts).unwrap();
    assert_eq!((eph.n_rocks, eph.n_epochs), (input.len(), jds.len()));
    assert_eq!(eph.epochs, jds.to_vec());
    let states = eph.states.as_ref().unwrap();

    let mut worst = 0.0f64;
    for (i, r) in input.iter().enumerate() {
        for (j, o) in observers.iter().enumerate() {
            let mut x = r.clone();
            x.propagate(&o.epoch, &k).unwrap();
            x.change_reference_plane("J2000").unwrap();
            x.to_ssb(&k).unwrap();
            let ob = x.observe(o).unwrap();
            let a = eph.get(i, j);
            let d = ((a.ra - ob.ra()) * ob.dec().cos()).hypot(a.dec - ob.dec());
            worst = worst.max(d);
            assert!((a.mag - ob.mag().unwrap()).abs() < 1e-9);
            let s = states[i * jds.len() + j];
            assert!((Vector3::new(s[0], s[1], s[2]) - x.position).norm() < 1e-9);
        }
        // The duplicate epoch gives identical results.
        assert_eq!(eph.get(i, 0), eph.get(i, 3));
    }
    // Well below a milliarcsecond.
    assert!(worst.to_degrees() * 3.6e6 < 0.1, "worst offset {} mas", worst.to_degrees() * 3.6e6);

    // Two-body ephemeris agrees with analytic propagation + observe exactly in structure.
    let eph2 = batch::ephemeris_rocks(&input, &observers, &k, &BatchOptions { method: Method::TwoBody, ..Default::default() }).unwrap();
    for (i, r) in input.iter().enumerate() {
        for (j, o) in observers.iter().enumerate() {
            let mut x = r.clone();
            x.analytic_propagate(&o.epoch).unwrap();
            x.change_reference_plane("J2000").unwrap();
            x.to_ssb(&k).unwrap();
            let ob = x.observe(o).unwrap();
            let a = eph2.get(i, j);
            assert!((a.ra - ob.ra()).abs() < 1e-12 && (a.dec - ob.dec()).abs() < 1e-12);
        }
    }
}

#[test]
fn perturber_cache_reproduces_kernel_and_falls_back() {
    let Some(k) = kernel() else {
        eprintln!("kernels not available; skipping");
        return;
    };
    let ids = SpiceSimulation::horizons_body_ids().unwrap();
    let cache = PerturberCache::build(&k, &ids, T0 - 20.0, T0 + 80.0).unwrap();
    let mut a = vec![[0.0; 6]; ids.len()];
    let mut b = vec![[0.0; 6]; ids.len()];
    for i in 0..2000 {
        let t = T0 - 20.0 + 100.0 * (i as f64 + 0.37) / 2000.0;
        assert!(cache.states(&ids, t, &mut a));
        k.barycentric_states_au(&ids, t, &mut b).unwrap();
        for (x, y) in a.iter().zip(&b) {
            let dp = ((x[0] - y[0]).powi(2) + (x[1] - y[1]).powi(2) + (x[2] - y[2]).powi(2)).sqrt();
            assert!(dp < 2e-13, "position error {:e} AU", dp);
        }
    }
    // Outside the span, or for a different body list, the cache declines.
    assert!(!cache.states(&ids, T0 + 81.0, &mut a));
    assert!(!cache.states(&ids[..3], T0, &mut a));

    // A simulation using the cache agrees with one reading the kernel.
    let rock = SpaceRock::from_kepler("x", 1.1, 0.3, 0.2, 1.0, 2.0, 0.5, tdb(T0), "J2000", "SSB").unwrap();
    let mut with = SpiceSimulation::horizons(&tdb(T0), &k).unwrap();
    with.set_perturber_cache(std::sync::Arc::new(cache));
    with.add(rock.clone()).unwrap();
    let mut without = SpiceSimulation::horizons(&tdb(T0), &k).unwrap();
    without.add(rock).unwrap();
    with.integrate_jd(T0 + 60.0, &k).unwrap();
    without.integrate_jd(T0 + 60.0, &k).unwrap();
    let d = (with.state.particles[0].position - without.state.particles[0].position).norm();
    assert!(d < 1e-10, "{:e}", d);
}

/// The Population forms run the same code as the rock forms: identical results for a population
/// in one plane and origin, both methods.
#[test]
fn population_forms_match_rock_forms() {
    let Some(k) = kernel() else {
        eprintln!("kernels not available; skipping");
        return;
    };
    for (plane, origin) in [("ECLIPJ2000", "SUN"), ("J2000", "SSB")] {
        let mut input = rocks(T0);
        for r in input.iter_mut() {
            r.change_reference_plane(plane).unwrap();
            if origin == "SSB" {
                r.to_ssb(&k).unwrap();
            } else {
                r.to_helio(&k).unwrap();
            }
        }
        let pop = Population::from_rocks(input.clone()).unwrap();
        let jds = [T0 + 30.0, T0 - 25.0, T0 + 3.5];
        let observers: Vec<Observer> = jds.iter().map(|&t| observer_at(t, &k)).collect();
        for method in [Method::NBody, Method::TwoBody] {
            let opts = BatchOptions { method, with_states: true, ..Default::default() };
            let a = batch::ephemeris(&pop, &observers, &k, &opts).unwrap();
            let b = batch::ephemeris_rocks(&input, &observers, &k, &opts).unwrap();
            assert_eq!(a.states, b.states);
            for i in 0..input.len() {
                for j in 0..jds.len() {
                    // Debug strings, so NaN fields compare equal.
                    assert_eq!(format!("{:?}", a.get(i, j)), format!("{:?}", b.get(i, j)));
                }
            }

            let t1 = tdb(T0 + 40.0);
            let mut p = pop.clone();
            batch::propagate(&mut p, &t1, &k, &opts).unwrap();
            let mut rs = input.clone();
            batch::propagate_batch(&mut rs, &t1, &k, &opts).unwrap();
            for (i, r) in rs.iter().enumerate() {
                assert_eq!(p.states[i], r.state(), "{plane} {origin} {method:?} {i}");
                assert_eq!(p.epochs[i], r.epoch.tdb().jd());
            }
        }
    }
}
