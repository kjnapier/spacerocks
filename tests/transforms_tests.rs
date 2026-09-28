use spacerocks::OrbitType;
use spacerocks::SpaceRock;
use spacerocks::Observer;
use spacerocks::Time;
use spacerocks::transforms::calc_conic_anomaly_from_mean_anomaly;
use spacerocks::transforms::calc_conic_anomaly_from_true_anomaly;
use spacerocks::transforms::calc_mean_anomaly_from_conic_anomaly;
use spacerocks::transforms::calc_true_anomaly_from_conic_anomaly;
use spacerocks::transforms::calc_true_anomaly_from_mean_anomaly;
use spacerocks::transforms::correct_for_ltt;
use spacerocks::constants::SPEED_OF_LIGHT;
use spacerocks::errors::OrbitError;
use spacerocks::constants::MU_BARY;

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;
    use nalgebra::Vector3;
    

    const EPSILON: f64 = 1e-10;

    /// Orbit Classification Tests
    #[test]
    fn test_orbit_type_classification() {
        let threshold = 1e-10;
        
        // Test circular orbit
        let e = 0.0;
        match OrbitType::from_eccentricity(e, EPSILON) {
            Ok(result) => assert_eq!(result, OrbitType::Circular),
            Err(_) => assert!(false, "Failed to classify circular orbit"),
        }

        // Test elliptical orbit
        let e = 0.5;
        match OrbitType::from_eccentricity(e, EPSILON) {
            Ok(result) => assert_eq!(result, OrbitType::Elliptical),
            Err(_) => assert!(false, "Failed to classify elliptical orbit"),
        }

        // Test parabolic orbit
        let e = 1.0;
        match OrbitType::from_eccentricity(e, EPSILON) {
            Ok(result) => assert_eq!(result, OrbitType::Parabolic),
            Err(_) => assert!(false, "Failed to classify parabolic orbit"),
        }

        // Test hyperbolic orbit
        let e = 1.5;
        match OrbitType::from_eccentricity(e, EPSILON) {
            Ok(result) => assert_eq!(result, OrbitType::Hyperbolic),
            Err(_) => assert!(false, "Failed to classify hyperbolic orbit"),
        }

        // Test invalid eccentricity
        let e = -0.1;
        match OrbitType::from_eccentricity(e, EPSILON) {
            Ok(_) => assert!(false, "Should reject negative eccentricity"),
            Err(err) => assert_eq!(err, OrbitError::NegativeEccentricity(e)),
        }
    }

    /// Testing conic anomaly from mean anomaly

    #[test]
    fn test_calc_conic_anomaly_from_mean_anomaly() {

        // Test circular orbits (e = 0)

        let e = 0.0;
        let test_anomalies = vec![
            0.0, 
            std::f64::consts::PI / 4.0,
            std::f64::consts::PI / 2.0,
            std::f64::consts::PI,
            3.0 * std::f64::consts::PI / 2.0,
            2.0 * std::f64::consts::PI
        ];
        
        for &M in &test_anomalies {
            match calc_conic_anomaly_from_mean_anomaly(e, M) {
                Ok(result) => assert!((result - M).abs() < EPSILON, 
                    "Circular orbit: for M = {}, expected {}, got {}", M, M, result),
                Err(_) => assert!(false, "Circular orbit should not fail for e = {}, M = {}", e, M),
            }
        }
    
        // Test elliptical orbits (0 < e < 1)

        let elliptical_es = vec![0.1, 0.5, 0.7, 0.9];
        for &e in &elliptical_es {
            for &M in &test_anomalies {
                match calc_conic_anomaly_from_mean_anomaly(e, M) {
                    Ok(E) => {
                        // Check Kepler's equation: M = E - e*sin(E)
                        let computed_M = E - e * E.sin();
                        assert!((computed_M - M).abs() < EPSILON, 
                            "Elliptical orbit: for e = {}, M = {}, got E = {}, which gives M = {}", 
                            e, M, E, computed_M);
                    },
                    Err(_) => assert!(false, "Elliptical orbit should not fail for e = {}, M = {}", e, M),
                }
            }
        }
    
        // Test parabolic orbit (e = 1)

        let e = 1.0;
        let parabolic_Ms = vec![-1.0, -0.5, 0.0, 0.5, 1.0];
        for &M in &parabolic_Ms {
            match calc_conic_anomaly_from_mean_anomaly(e,M) {
                Ok(B) => {
                    // Verify Barker's equation: M = tan(f/2)/2 + tan³(f/2)/6
                    // where tan(f/2) = B
                    let computed_M = B/2.0 + B.powi(3)/6.0;
                    assert!((computed_M - M).abs() < EPSILON,
                        "Parabolic orbit: for M = {}, got B = {}, which gives M = {}",
                        M, B, computed_M);
                },
                Err(_) => assert!(false, "Parabolic orbit should not fail for M = {}", M),
        }
    }
    
        // Test hyperbolic orbits (e > 1)

        let hyperbolic_es = vec![1.1, 1.5, 2.0];
        let hyperbolic_Ms = vec![-2.0, -1.0, 0.0, 1.0, 2.0];
        for &e in &hyperbolic_es {
            for &M in &hyperbolic_Ms {
                match calc_conic_anomaly_from_mean_anomaly(e, M) {
                    Ok(H) => {
                        // Check Kepler's equation for hyperbolic orbits: M = e*sinh(H) - H
                        let computed_M = e * H.sinh() - H;
                        assert!((computed_M - M).abs() < EPSILON, 
                            "Hyperbolic orbit: for e = {}, M = {}, got H = {}, which gives M = {}", 
                            e, M, H, computed_M);
                    },
                    Err(_) => assert!(false, "Hyperbolic orbit should not fail for e = {}, M = {}", e, M),
                }
            }
        }
    
        // Test invalid eccentricities

        let e = -0.1;
        match calc_conic_anomaly_from_mean_anomaly(e, 0.0) {
            Ok(_) => assert!(false, "Should reject negative eccentricity"),
            Err(err) => assert_eq!(err, OrbitError::NegativeEccentricity(e)),
        }
    
    }

    /// Testing conic anomaly from true anomaly
    
    #[test]
    fn test_calc_conic_anomaly_from_true_anomaly() {
    
        // Test circular orbits (e = 0)
        
        let e = 0.0;
        let test_anomalies = vec![
            0.0, 
            std::f64::consts::PI / 4.0,
            std::f64::consts::PI / 2.0,
            std::f64::consts::PI,
            3.0 * std::f64::consts::PI / 2.0,
            2.0 * std::f64::consts::PI
        ];
        
        for &f in &test_anomalies {
            match calc_conic_anomaly_from_true_anomaly(e, f) {
                Ok(result) => assert!((result - f).abs() < EPSILON, 
                    "Circular orbit: for f = {}, expected {}, got {}", f, f, result),
                Err(_) => assert!(false, "Circular orbit should not fail for e = {}, f = {}", e, f),
            }
        }
    
        // Test elliptical orbits (0 < e < 1)
        
        let elliptical_es = vec![0.1, 0.5, 0.7, 0.9];
        for &e in &elliptical_es {
            for &f in &test_anomalies {
                match calc_conic_anomaly_from_true_anomaly(e, f) {
                    Ok(E) => {
                        // Verify with the relation: f = 2*atan(sqrt((1+e)/(1-e)) * tan(E/2))
                        let computed_f = 2.0 * ((1.0 + e).sqrt() * (E/2.0).sin()).atan2((1.0 - e).sqrt() * (E/2.0).cos());
                        assert!((computed_f - f).abs() < EPSILON, 
                            "Elliptical orbit: for e = {}, f = {}, got E = {}", e, f, E);
                    },
                    Err(_) => assert!(false, "Elliptical orbit should not fail for e = {}, f = {}", e, f),
                }
            }
        }
    
        // Test parabolic orbit (e = 1)
        
        let e = 1.0;
        let parabolic_fs = vec![-1.0, -0.5, 0.0, 0.5, 1.0];
        for &f in &parabolic_fs {
            match calc_conic_anomaly_from_true_anomaly(e, f) {
                Ok(B) => {
                    // For parabolic orbits, B = tan(f/2)
                    assert!((B - (f/2.0).tan()).abs() < EPSILON,
                        "Parabolic orbit: for f = {}, got B = {}, expected tan(f/2) = {}",
                        f, B, (f/2.0).tan());
                },
                Err(_) => assert!(false, "Parabolic orbit should not fail for f = {}", f),
            }
        }
    
        // Test hyperbolic orbits (e > 1)
        
        let hyperbolic_es = vec![1.1, 1.5, 2.0];
        let hyperbolic_fs = vec![-1.0, -0.5, 0.0, 0.5, 1.0];
        for &e in &hyperbolic_es {
            for &f in &hyperbolic_fs {
                match calc_conic_anomaly_from_true_anomaly(e, f) {
                    Ok(H) => {
                        // Verify using the inverse relation: tan(f/2) = sqrt((e+1)/(e-1)) * tanh(H/2)
                        let tan_f_half = ((e + 1.0)/(e - 1.0)).sqrt() * (H/2.0).tanh();
                        assert!((tan_f_half - (f/2.0).tan()).abs() < EPSILON, 
                            "Hyperbolic orbit: for e = {}, f = {}, got H = {}, which gives tan(f/2) = {}", 
                            e, f, H, tan_f_half);
                    },
                    Err(_) => assert!(false, "Hyperbolic orbit should not fail for e = {}, f = {}", e, f),
                }
            }
        }
    }

    #[test]
    fn test_state_from_kepler_elements() {
        let epoch = Time::new(2451545.0, "tdb", "jd").unwrap();
        let mu = spacerocks::Origin::SSB.mu();

        // Circular orbit at periapsis: position (a, 0, 0), velocity (0, sqrt(mu/a), 0)
        let r = SpaceRock::from_kepler("c", 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, epoch.clone(), "J2000", "SSB").unwrap();
        assert!((r.position - Vector3::new(1.0, 0.0, 0.0)).norm() < EPSILON);
        assert!((r.velocity - Vector3::new(0.0, (mu / 1.0f64).sqrt(), 0.0)).norm() < EPSILON);

        // Elliptical and hyperbolic orbits at periapsis: r = q, v from vis-viva
        for (q, e) in [(1.0f64, 0.5f64), (2.0, 2.0)] {
            let a = q / (1.0 - e);
            let r = SpaceRock::from_kepler("e", q, e, 0.0, 0.0, 0.0, 0.0, epoch.clone(), "J2000", "SSB").unwrap();
            let v = (mu * (2.0 / q - 1.0 / a)).sqrt();
            assert!((r.position.norm() - q).abs() < EPSILON);
            assert!((r.velocity.norm() - v).abs() < EPSILON);
        }

        // 90 degree inclination: orbit in the x-z plane
        let r = SpaceRock::from_kepler("i", 1.0, 0.0, PI / 2.0, 0.0, 0.0, 0.0, epoch.clone(), "J2000", "SSB").unwrap();
        assert!(r.position.y.abs() < EPSILON && r.velocity.y.abs() < EPSILON);

        // Round trip through the element accessors
        let (q, e, inc, arg, node, nu) = (2.5, 0.3, 0.4, 1.1, 0.7, 0.9);
        let r = SpaceRock::from_kepler("rt", q, e, inc, arg, node, nu, epoch, "J2000", "SSB").unwrap();
        assert!((r.q() - q).abs() < 1e-12);
        assert!((r.e() - e).abs() < 1e-12);
        assert!((r.inc() - inc).abs() < 1e-12);
        assert!((r.arg() - arg).abs() < 1e-12);
        assert!((r.node() - node).abs() < 1e-12);
        assert!((r.true_anomaly() - nu).abs() < 1e-12);
    }

    #[test]
    fn test_light_time_correction_matches_two_body() {
        let epoch = Time::new(2460000.5, "tdb", "jd").unwrap();
        let rock = SpaceRock::from_kepler("neo", 0.9, 0.3, 0.2, 1.0, 0.5, 0.3, epoch.clone(), "J2000", "SSB").unwrap();
        let obs = Observer::from_xyz(rock.position + Vector3::new(0.15, -0.1, 0.05), Some(Vector3::new(0.0, 0.017, 0.0)),
            epoch.clone(), spacerocks::ReferencePlane::J2000, spacerocks::Origin::SSB, None);
        let cr = correct_for_ltt(&rock, &obs);
        let tau = cr.position.norm() / SPEED_OF_LIGHT;
        let exact = rock.analytic_at(&(epoch - tau)).unwrap();
        assert!((cr.position - (exact.position - obs.position)).norm() < 1e-10);
        assert!((cr.velocity - (exact.velocity - obs.velocity.unwrap())).norm() < 1e-10);
    }

    #[test]
    fn test_elements_round_trip_including_degenerate_orbits() {
        let epoch = Time::new(2451545.0, "tdb", "jd").unwrap();
        let two_pi = 2.0 * PI;
        let wrap = |x: f64| x.rem_euclid(two_pi);
        let close = |a: f64, b: f64| {
            let d = (wrap(a) - wrap(b)).abs();
            d.min(two_pi - d) < 1e-9
        };
        // (q, e, inc, arg, node, nu)
        let cases = [
            (2.5, 0.3, 0.4, 1.1, 0.7, 0.9),
            (2.5, 0.3, 0.4, 5.9, 4.0, 4.5),   // quadrants > pi
            (1.0, 1.7, 2.8, 0.5, 3.3, -0.6),  // retrograde hyperbolic
            (40.0, 0.05, 1.0e-3, 2.0, 1.0, 3.5),
        ];
        for &(q, e, inc, arg, node, nu) in &cases {
            let r = SpaceRock::from_kepler("x", q, e, inc, arg, node, nu, epoch.clone(), "J2000", "SSB").unwrap();
            assert!((r.q() - q).abs() < 1e-9 * q && (r.e() - e).abs() < 1e-11 && (r.inc() - inc).abs() < 1e-11);
            assert!(close(r.arg(), arg) && close(r.node(), node) && close(r.true_anomaly(), nu), "{:?}", (q, e, inc, arg, node, nu));
        }
        // Equatorial: node undefined (0); arg is then the longitude of perihelion.
        let r = SpaceRock::from_kepler("eq", 2.0, 0.2, 0.0, 1.0, 0.5, 2.0, epoch.clone(), "J2000", "SSB").unwrap();
        assert_eq!(r.node(), 0.0);
        assert!(close(r.arg(), 1.5) && close(r.true_anomaly(), 2.0));
        // Circular: perihelion undefined (arg 0); true anomaly is the argument of latitude.
        let r = SpaceRock::from_kepler("c", 2.0, 0.0, 0.3, 0.0, 0.5, 4.0, epoch.clone(), "J2000", "SSB").unwrap();
        assert!(close(r.arg(), 0.0) && close(r.node(), 0.5) && close(r.true_anomaly(), 4.0));
        // Circular and equatorial: true longitude.
        let r = SpaceRock::from_kepler("ce", 2.0, 0.0, 0.0, 0.0, 0.0, 4.0, epoch, "J2000", "SSB").unwrap();
        assert!(!r.true_anomaly().is_nan() && close(r.true_anomaly(), 4.0));
    }

    #[test]
    fn test_apparent_magnitude_geometry() {
        use spacerocks::spacerock::hg_magnitude;
        let epoch = Time::new(2460000.5, "tdb", "jd").unwrap();
        let obs = |x: f64, y: f64| Observer::from_xyz(Vector3::new(x, y, 0.0), Some(Vector3::zeros()),
            epoch.clone(), spacerocks::ReferencePlane::J2000, spacerocks::Origin::SUN, None);
        let rock = |x: f64, y: f64| {
            let mut r = SpaceRock::from_xyz("m", x, y, 0.0, 0.0, 0.0, 0.0, epoch.clone(), "J2000", "SUN").unwrap();
            r.set_absolute_magnitude(7.0);
            r
        };
        // Opposition: phase 0, r = 2, delta = 1  ->  H + 5 log10(2)
        let o = rock(2.0, 0.0).observe(&obs(1.0, 0.0)).unwrap();
        assert!((o.mag.unwrap() - (7.0 + 5.0 * 2f64.log10())).abs() < 1e-6);
        // Quadrature-like geometry: phase = atan(1/2) (26.57 deg), r = 2, delta = sqrt(5)
        let o = rock(0.0, 2.0).observe(&obs(1.0, 0.0)).unwrap();
        let expected = hg_magnitude(7.0, 0.15, 2.0, 5f64.sqrt(), (0.5f64).atan());
        assert!((o.mag.unwrap() - expected).abs() < 1e-6);
        // With the Sun offset from the origin, the observer's sun_position is used.
        let shift = Vector3::new(0.01, -0.02, 0.0);
        let mut r = rock(2.01, -0.02);
        r.origin = spacerocks::Origin::SSB;
        let mut ob = obs(1.01, -0.02).with_sun_position(shift);
        ob.origin = spacerocks::Origin::SSB;
        let o = r.observe(&ob).unwrap();
        assert!((o.mag.unwrap() - (7.0 + 5.0 * 2f64.log10())).abs() < 1e-6);
    }
}
