//! `state` kernels and `Population` agree with the per-rock `SpaceRock` methods.

use nalgebra::Vector3;
use spacerocks::observing::Observer;
use spacerocks::{state, Origin, Population, ReferencePlane, SpaceRock, Time};

fn rocks(n: usize, epoch: &Time) -> Vec<SpaceRock> {
    (0..n)
        .map(|i| {
            let x = i as f64;
            let e = if i % 5 == 0 { 0.0 } else { (0.13 * x).sin().abs() * 0.95 };
            let inc = if i % 4 == 0 { 0.0 } else { 0.3 * x % 3.0 };
            SpaceRock::from_kepler(&format!("r{}", i), 1.0 + x, e, inc, 0.7 * x % 6.2, 1.1 * x % 6.2, 0.9 * x % 6.2 - 3.1, epoch.clone(), "ECLIPJ2000", "SUN").unwrap()
        })
        .collect()
}

#[test]
fn elements_match_individual_methods() {
    let t = Time::new(2460000.5, "TDB", "JD").unwrap();
    for r in rocks(50, &t) {
        let el = r.elements();
        assert_eq!(el.a, r.a());
        assert_eq!(el.e, r.e());
        assert_eq!(el.q, r.q());
        assert_eq!(el.inc, r.inc());
        assert_eq!(el.node, r.node());
        assert_eq!(el.arg, r.arg());
        assert_eq!(el.true_anomaly, r.true_anomaly());
        assert_eq!(el.conic_anomaly, r.conic_anomaly());
        assert_eq!(el.mean_anomaly, r.mean_anomaly());
    }
}

#[test]
fn from_kepler_round_trips() {
    let mu = Origin::SUN.mu();
    let s = state::from_kepler(2.5, 0.3, 0.4, 1.0, 2.0, 0.5, mu).unwrap();
    let el = state::elements(&s, mu);
    assert!((el.q - 2.5).abs() < 1e-12);
    assert!((el.e - 0.3).abs() < 1e-12);
    assert!((el.inc - 0.4).abs() < 1e-12);
    assert!((el.arg - 1.0).abs() < 1e-12);
    assert!((el.node - 2.0).abs() < 1e-12);
    assert!((el.true_anomaly - 0.5).abs() < 1e-12);
}

#[test]
fn push_rotates_into_the_population_plane_and_checks_origin() {
    let t = Time::new(2460000.5, "TDB", "JD").unwrap();
    let mut pop = Population::new(ReferencePlane::J2000, Origin::SUN);
    let r = rocks(1, &t).remove(0);
    pop.push(r.clone()).unwrap();
    let mut expected = r.clone();
    expected.change_reference_plane("J2000").unwrap();
    assert_eq!(pop.get(0).unwrap(), expected);

    let ssb = SpaceRock::from_xyz("b", 1.0, 0.0, 0.0, 0.0, 0.017, 0.0, t, "J2000", "SSB").unwrap();
    assert!(pop.push(ssb).is_err());
    assert_eq!(pop.len(), 1);
}

#[test]
fn empty_population_takes_the_first_rocks_metadata() {
    let t = Time::new(2460000.5, "TDB", "JD").unwrap();
    let rs = rocks(10, &t);
    let pop = Population::from_rocks(rs.clone()).unwrap();
    assert_eq!(pop.reference_plane, ReferencePlane::ECLIPJ2000);
    assert_eq!(pop.origin, Origin::SUN);
    assert_eq!(pop.to_rocks(), rs);
    assert_eq!(pop.index_of("r3"), Some(3));
    let f = pop.filter(&(0..10).map(|i| i % 2 == 0).collect::<Vec<_>>()).unwrap();
    assert_eq!(f.names, vec!["r0", "r2", "r4", "r6", "r8"]);
}

#[test]
fn population_operations_match_spacerock() {
    let t0 = Time::new(2460000.5, "TDB", "JD").unwrap();
    let t1 = Time::new(2460321.25, "TDB", "JD").unwrap();
    // mixed epochs
    let mut rs = rocks(40, &t0);
    for (i, r) in rs.iter_mut().enumerate() {
        r.epoch = Time::new(2460000.5 + i as f64, "TDB", "JD").unwrap();
    }
    let mut pop = Population::from_rocks(rs.clone()).unwrap();

    pop.analytic_propagate(&t1).unwrap();
    pop.change_reference_plane(&ReferencePlane::J2000).unwrap();
    for r in rs.iter_mut() {
        r.analytic_propagate(&t1).unwrap();
        r.change_reference_plane("J2000").unwrap();
    }
    assert_eq!(pop.to_rocks(), rs);

    let obs = Observer::from_xyz(Vector3::new(0.3, -0.9, 0.1), Some(Vector3::new(0.015, 0.005, 0.0)), t1, ReferencePlane::J2000, Origin::SUN, None);
    let apps = pop.apparent(&obs).unwrap();
    for (a, r) in apps.iter().zip(&rs) {
        let b = r.apparent(&obs).unwrap();
        assert_eq!(a.ra, b.ra);
        assert_eq!(a.dec, b.dec);
        assert_eq!(a.range_rate, b.range_rate);
    }

    let wrong_epoch = Observer::from_xyz(Vector3::new(0.3, -0.9, 0.1), Some(Vector3::zeros()), t0, ReferencePlane::J2000, Origin::SUN, None);
    assert!(pop.apparent(&wrong_epoch).is_err());
}
