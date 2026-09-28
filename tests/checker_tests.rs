//! Tests of `spacerocks::checker`. The end-to-end tests need `de440s.bsp` (or `de440.bsp`),
//! `sb441-n16.bsp` and an Earth orientation file (`earth_*.bpc`) in `SPACEROCKS_KERNELS`;
//! skipped otherwise. Checks against real MPC data and the full MPCORB are in
//! `validation/checker/`.

use std::path::PathBuf;

use spacerocks::checker::{self, mpc, Catalog, CheckOptions, MpcElements};
use spacerocks::orbfit::{self, Astrometry, FitFlag, FitOptions, OrbitFit};
use spacerocks::SpiceKernel;

const ARCSEC: f64 = std::f64::consts::PI / (180.0 * 3600.0);

fn kernel() -> Option<SpiceKernel> {
    let dir = PathBuf::from(std::env::var("SPACEROCKS_KERNELS").ok()?);
    let mut k = SpiceKernel::new();
    let planets = if dir.join("de440.bsp").exists() { "de440.bsp" } else { "de440s.bsp" };
    k.load(dir.join(planets)).ok()?;
    k.load(dir.join("sb441-n16.bsp")).ok()?;
    for e in std::fs::read_dir(&dir).ok()? {
        let p = e.ok()?.path();
        if p.extension().map(|x| x == "bpc").unwrap_or(false) {
            k.load(&p).ok()?;
        }
    }
    Some(k)
}

#[test]
fn packed_epochs_and_designations() {
    // 2025 May 5.0 and 1999 Dec 31.0 (day 'V' = 31).
    assert_eq!(mpc::unpack_epoch("K2555"), Some(2460800.5));
    assert_eq!(mpc::unpack_epoch("J99CV"), Some(2451543.5));
    assert_eq!(mpc::unpack_epoch("X2555"), None);
    assert_eq!(mpc::unpack_designation("00001"), "1");
    assert_eq!(mpc::unpack_designation("03666"), "3666");
    assert_eq!(mpc::unpack_designation("A0345"), "100345");
    assert_eq!(mpc::unpack_designation("~0000"), "620000");
    assert_eq!(mpc::unpack_designation("K07Tf8A"), "2007 TA418");
    assert_eq!(mpc::unpack_designation("J95X00A"), "1995 XA");
}

/// A line of MPCORB.DAT, laid out column by column as the MPC documents it.
fn mpcorb_dat_line() -> String {
    let mut s = vec![b' '; 202];
    let mut put = |col: usize, text: &str| {
        // 1-based column where `text` starts.
        s[col - 1..col - 1 + text.len()].copy_from_slice(text.as_bytes());
    };
    put(1, "00001");
    put(9, " 3.34");
    put(15, " 0.15");
    put(21, "K2555");
    put(27, "188.70269");
    put(38, " 73.27343");
    put(49, " 80.25221");
    put(60, " 10.58780");
    put(71, "0.0794013");
    put(81, " 0.21424651");
    put(93, "  2.7660512");
    put(106, "0");
    put(108, "E2024-V47");
    put(167, "(1) Ceres");
    put(195, "20241101");
    String::from_utf8(s).unwrap()
}

#[test]
fn reads_mpcorb_dat_lines() {
    let dir = std::env::temp_dir().join("spacerocks_checker_dat");
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("MPCORB.DAT");
    std::fs::write(&path, format!("MPCORB header text\n-----------\n{}\n\n", mpcorb_dat_line())).unwrap();
    let els = mpc::read_mpcorb(&path).unwrap();
    assert_eq!(els.len(), 1);
    let e = &els[0];
    assert_eq!(e.name, "1");
    assert_eq!(e.epoch_tt, 2460800.5);
    assert_eq!((e.a, e.e, e.inc, e.node, e.peri, e.m), (2.7660512, 0.0794013, 10.5878, 80.25221, 73.27343, 188.70269));
    assert_eq!((e.h, e.g, e.u), (3.34, 0.15, 0.0));
    assert_eq!(e.last_obs, 2460615.5);
}

#[test]
fn elements_to_state_matches_the_orbit() {
    // Ceres' elements: the state has the right semimajor axis, and a quarter period later the
    // mean anomaly has advanced by 90 degrees.
    let (a, e) = (2.7660512, 0.0794013);
    let d = std::f64::consts::PI / 180.0;
    let (p, v) = mpc::elements_to_state(a, e, 10.5878 * d, 80.25221 * d, 73.27343 * d, 188.70269 * d).unwrap();
    let mu = spacerocks::Origin::SUN.mu();
    let a2 = 1.0 / (2.0 / p.norm() - v.norm_squared() / mu);
    assert!((a2 - a).abs() < 1e-12, "{} vs {}", a2, a);
    let h = p.cross(&v);
    let inc = (h.z / h.norm()).acos() / d;
    assert!((inc - 10.5878).abs() < 1e-10);
}

#[test]
fn parses_mpc_orb_json() {
    let text = std::fs::read_to_string("tests/data/mpc_orb_3666.json").unwrap();
    let o = mpc::parse_mpc_orb(&serde_json::from_str(&text).unwrap()).unwrap();
    assert_eq!(o.name, "3666");
    // MJD 61200 TT.
    assert!((o.epoch_tdb - 2461200.5).abs() < 1e-6);
    // Heliocentric distance is preserved by the rotation to the equator.
    let r = (o.helio_j2000[0].powi(2) + o.helio_j2000[1].powi(2) + o.helio_j2000[2].powi(2)).sqrt();
    let r_ecl: f64 = [2.93294851396238f64, 1.74610708363794, -0.140989125026942].iter().map(|x| x * x).sum::<f64>().sqrt();
    assert!((r - r_ecl).abs() < 1e-14);
    // The covariance is symmetric, and its trace (rotation invariant) is the ecliptic one's.
    let c = &o.covariance;
    for i in 0..6 {
        for j in 0..6 {
            assert!((c[i * 6 + j] - c[j * 6 + i]).abs() <= 1e-12 * c[i * 6 + i].abs().max(c[j * 6 + j].abs()));
        }
    }
    let tr_pos = c[0] + c[7] + c[14];
    let tr_ecl = 2.552022069421863e-15 + 3.253676125691538e-15 + 2.339998991818776e-15;
    assert!((tr_pos - tr_ecl).abs() < 1e-6 * tr_ecl);
    assert_eq!(o.u, 0.0);
}

#[test]
fn runoff_brackets_the_mpc_definition() {
    // U = int(ln(runoff) / 1.49) + 1: the runoff returned for U maps back to U.
    for u in 1..=9 {
        let r = mpc::runoff_from_u(u as f64) / ARCSEC;
        assert_eq!(((r * 0.999).ln() / 1.49).floor() as i32 + 1, u);
    }
    assert!((mpc::runoff_from_u(0.0) / ARCSEC - 1.0).abs() < 1e-12);
}

fn elements(name: &str, a: f64, e: f64, inc: f64, node: f64, peri: f64, m: f64, u: f64) -> MpcElements {
    MpcElements { name: name.into(), epoch_tt: 2460800.5, a, e, inc, node, peri, m, h: 15.0, g: 0.15, u, last_obs: 2460700.5 }
}

/// Noise-free detections of catalog object `i` from Rubin (X05) at TDB `times`, with the
/// catalog's own (N-body) orbit.
fn detect(cat: &Catalog, i: usize, times: &[f64], k: &SpiceKernel) -> Astrometry {
    let mut fit = OrbitFit::failed(FitFlag::Converged);
    fit.epoch = cat.epoch[i];
    fit.state = cat.state[i];
    let codes = vec!["X05".to_string(); times.len()];
    let obs = orbfit::observer_positions(&codes, times, k).unwrap();
    let p = orbfit::predict(&fit, times, &obs, k, &FitOptions::default()).unwrap();
    let mut a = Astrometry::new(times.to_vec(), p.iter().map(|x| x.ra).collect(), p.iter().map(|x| x.dec).collect(), obs).unwrap();
    a = a.with_sigma(&vec![0.2 * ARCSEC; times.len()], &vec![0.2 * ARCSEC; times.len()]).unwrap();
    a
}

fn append(a: &mut Astrometry, b: &Astrometry) {
    a.epoch.extend(&b.epoch);
    a.ra.extend(&b.ra);
    a.dec.extend(&b.dec);
    a.sigma_ra.extend(&b.sigma_ra);
    a.sigma_dec.extend(&b.sigma_dec);
    a.observer.extend(&b.observer);
}

/// A small MPCORB-like catalog: three named objects among 2,000 main-belt orbits spread in
/// mean anomaly and node.
fn test_catalog(k: &SpiceKernel) -> Catalog {
    let mut els = vec![
        elements("mba", 2.6, 0.12, 8.0, 40.0, 100.0, 20.0, 0.0),
        elements("neo", 1.3, 0.35, 12.0, 200.0, 30.0, 60.0, 3.0),
        elements("tno", 44.0, 0.05, 3.0, 120.0, 10.0, 250.0, 1.0),
    ];
    for n in 0..2000 {
        let x = n as f64;
        els.push(elements(&format!("bg{}", n), 2.2 + (x * 0.37) % 1.1, (x * 0.013) % 0.25, (x * 0.7) % 20.0, (x * 7.3) % 360.0, (x * 11.9) % 360.0, (x * 3.7) % 360.0, (x % 10.0).min(9.0)));
    }
    Catalog::from_mpc_elements(&els, k).unwrap()
}

#[test]
fn identifies_detections_of_catalog_objects() {
    let Some(k) = kernel() else { return };
    let cat = test_catalog(&k);
    let t0 = cat.epoch[0];
    // Near the epoch (two-body from the orbits), and 200 days later (integrated first).
    let times = [t0 + 3.0, t0 + 3.02, t0 + 5.0, t0 + 200.0, t0 + 200.02];
    let mut det = Astrometry::default();
    let mut truth = Vec::new();
    for i in 0..3 {
        append(&mut det, &detect(&cat, i, &times, &k));
        truth.extend(std::iter::repeat(i).take(times.len()));
    }
    let found = checker::check(&cat, &det, &[], &k, &CheckOptions::default()).unwrap();
    assert!(found.failed.is_empty());
    for (j, &i) in truth.iter().enumerate() {
        let first = found.matches.iter().find(|m| m.detection == j).unwrap_or_else(|| panic!("detection {} not matched", j));
        assert_eq!(first.object, i, "detection {}", j);
        assert!(first.consistent && first.distance < 0.1, "detection {}: distance {}", j, first.distance);
        assert!(first.separation < 0.01 * ARCSEC, "detection {}: {}\"", j, first.separation / ARCSEC);
        assert!(first.mag.is_finite() && first.delta > 0.0);
    }

    // Moved 20" (100 sigma for the U = 0 main-belt object): reported (within the 60" radius),
    // but not consistent.
    let mut moved = detect(&cat, 0, &times, &k);
    for d in moved.dec.iter_mut() {
        *d += 20.0 * ARCSEC;
    }
    let found = checker::check(&cat, &moved, &[], &k, &CheckOptions::default()).unwrap();
    for j in 0..times.len() {
        let m = found.matches.iter().find(|m| m.detection == j && m.object == 0).expect("reported within the radius");
        assert!(!m.consistent && m.distance > 50.0, "distance {}", m.distance);
        assert!((m.offset[1] / ARCSEC - 20.0).abs() < 0.1);
    }
    // With radius 0 only consistent pairs are reported.
    let opts = CheckOptions { radius: 0.0, ..CheckOptions::default() };
    let found = checker::check(&cat, &moved, &[], &k, &opts).unwrap();
    assert!(found.matches.iter().all(|m| m.consistent && m.object != 0));
}

#[test]
fn covariance_orbits_and_snapshots() {
    let Some(k) = kernel() else { return };
    let base = test_catalog(&k);
    let t0 = base.epoch[0];
    // Object 0 as a "fit" with a covariance: 1000 km in position, 1 cm/s in velocity.
    let mut fit = OrbitFit::failed(FitFlag::Converged);
    fit.epoch = t0;
    fit.state = base.state[0];
    fit.covariance = vec![0.0; 36];
    for i in 0..3 {
        fit.covariance[i * 7] = (1000.0 / 1.495978707e8f64).powi(2);
        fit.covariance[(i + 3) * 7] = (1e-5 / 1.495978707e8 * 86400.0f64).powi(2);
    }
    let mut cat = Catalog::from_fits(&["fitted".to_string()], &[fit], None).unwrap();
    cat.extend(&base.select(&(1..base.len()).collect::<Vec<_>>()));
    let times = [t0 + 150.0, t0 + 150.03];
    let mut det = detect(&base, 0, &times, &k);
    let found = checker::check(&cat, &det, &[], &k, &CheckOptions::default()).unwrap();
    let m = found.matches.iter().find(|m| m.object == 0).unwrap();
    assert!(m.from_covariance && m.consistent);
    let (smaj, _, pa) = m.ellipse();
    assert!(smaj > 0.3 * ARCSEC, "{}", smaj / ARCSEC);
    // Displaced by 1.5 sigma along the ellipse's major axis: still consistent.
    let pa = pa.to_radians();
    for j in 0..times.len() {
        det.ra[j] += 1.5 * smaj * pa.sin() / det.dec[j].cos();
        det.dec[j] += 1.5 * smaj * pa.cos();
    }
    let found = checker::check(&cat, &det, &[], &k, &CheckOptions::default()).unwrap();
    let m = found.matches.iter().find(|m| m.object == 0 && m.detection == 0).unwrap();
    assert!(m.consistent && m.distance > 0.3, "distance {}", m.distance);

    // A snapshot near the detections gives the same answers, and survives save/load.
    let det = detect(&base, 1, &times, &k);
    let plain = checker::check(&base, &det, &[], &k, &CheckOptions::default()).unwrap();
    let mut snap = base.clone();
    snap.make_snapshot(t0 + 148.0, &k, &Default::default(), 500);
    let path = std::env::temp_dir().join("spacerocks_checker_test.srcat");
    snap.save(&path).unwrap();
    let snap = Catalog::load(&path).unwrap();
    assert_eq!(snap.len(), base.len());
    assert_eq!(snap.snapshot_epoch, t0 + 148.0);
    assert_eq!(snap.last_obs, base.last_obs);
    let with = checker::check(&snap, &det, &[], &k, &CheckOptions::default()).unwrap();
    let a: Vec<(usize, usize, bool)> = plain.matches.iter().map(|m| (m.detection, m.object, m.consistent)).collect();
    let b: Vec<(usize, usize, bool)> = with.matches.iter().map(|m| (m.detection, m.object, m.consistent)).collect();
    assert_eq!(a, b);
    for (x, y) in plain.matches.iter().zip(&with.matches) {
        assert!((x.distance - y.distance).abs() < 1e-3 * (1.0 + x.distance), "{} vs {}", x.distance, y.distance);
    }
}
