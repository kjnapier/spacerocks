//! Tests of `spacerocks::orbfit`. The fits need `de440s.bsp` (or `de440.bsp`), `sb441-n16.bsp`
//! and an Earth orientation file (`earth_*.bpc`) in `SPACEROCKS_KERNELS`; skipped otherwise.
//! Agreement with layup itself is checked by the scripts in `validation/layup/`.

use std::path::PathBuf;

use nalgebra::Vector3;

use spacerocks::assist::SpiceSimulation;
use spacerocks::orbfit::{self, build_sequence, determine_orbit, fit_orbit, gauss::select_triplet, Astrometry, FitFlag, FitOptions};
use spacerocks::{Observatory, SpaceRock, SpiceKernel, Time};

const C: f64 = 173.1446326742403;

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

/// Noise-free astrometry of `state` (barycentric J2000 at TDB `epoch`) from Rubin (X05) at
/// the TDB Julian dates `times`, with ASSIST's force model and light-time correction.
fn synthesize(k: &SpiceKernel, epoch: f64, state: [f64; 6], nongrav: [f64; 3], times: &[f64]) -> Astrometry {
    let site = Observatory::from_obscode("X05").unwrap();
    let t0 = Time::new(epoch, "tdb", "jd").unwrap();
    let mut a = Astrometry::default();
    for &t in times {
        let obs = site.at(&Time::new(t, "tdb", "jd").unwrap(), "J2000", "SSB", k).unwrap();
        let mut sim = SpiceSimulation::horizons(&t0, k).unwrap();
        sim.integrator.set_epsilon(1e-12);
        let mut rock = SpaceRock::from_xyz("rock", state[0], state[1], state[2], state[3], state[4], state[5], t0.clone(), "J2000", "SSB").unwrap();
        if nongrav != [0.0; 3] {
            rock.set_nongrav(nongrav[0], nongrav[1], nongrav[2]);
        }
        sim.add(rock).unwrap();
        let mut lt = 0.0;
        let mut rho = Vector3::zeros();
        for _ in 0..6 {
            sim.integrate_rel(t - lt - epoch, k).unwrap();
            rho = sim.state.particles[0].position - obs.position;
            lt = rho.norm() / C;
        }
        let rho = rho / rho.norm();
        a.epoch.push(t);
        a.ra.push(rho.y.atan2(rho.x).rem_euclid(std::f64::consts::TAU));
        a.dec.push(rho.z.asin());
        a.observer.push([obs.position.x, obs.position.y, obs.position.z]);
        a.sigma_ra.push(orbfit::DEFAULT_SIGMA);
        a.sigma_dec.push(orbfit::DEFAULT_SIGMA);
    }
    a
}

/// A main-belt object observed on pairs of nights over one opposition and again a year later.
fn main_belt_times(start: f64) -> Vec<f64> {
    let mut t = Vec::new();
    for night in [0.0, 1.0, 4.0, 9.0, 16.0, 30.0, 45.0, 60.0, 400.0, 401.0, 410.0] {
        t.push(start + night);
        t.push(start + night + 0.03);
    }
    t
}

const EPOCH: f64 = 2460400.5;
const MBA: [f64; 6] = [-2.0443, 1.3172, 0.6531, -0.006381, -0.008037, -0.003028];

#[test]
fn sequence_puts_longest_arc_first_then_nearest() {
    // arcs: [0, 10], [200, 300] (longest), [400, 405], [700, 701]
    let t = [0.0, 10.0, 200.0, 250.0, 300.0, 400.0, 405.0, 700.0, 701.0];
    let seq = build_sequence(&t, 90.0);
    assert_eq!(seq, vec![vec![2, 3, 4], vec![5, 6], vec![0, 1], vec![7, 8]]);
    // input order does not matter
    let shuffled = [300.0, 0.0, 701.0, 250.0, 10.0, 200.0];
    assert_eq!(build_sequence(&shuffled, 90.0), vec![vec![5, 3, 0], vec![1, 4], vec![2]]);
}

#[test]
fn triplet_targets_span_and_keeps_middle_central() {
    let target = orbfit::gauss::gauss_target_days();
    let t: Vec<f64> = (0..200).map(|i| i as f64).collect();
    let idx: Vec<usize> = (0..200).collect();
    let (a, b, c) = select_triplet(&t, &idx, target).unwrap();
    assert!(((t[c] - t[a]) - target).abs() <= 1.0);
    assert!((t[b] - 0.5 * (t[a] + t[c])).abs() <= 1.0);
    // a short arc is left to the caller's first/middle/last fallback
    assert!(select_triplet(&t[..20], &idx[..20], target).is_none());
}

#[test]
fn residuals_of_the_true_orbit_vanish() {
    let Some(k) = kernel() else { return };
    let a = synthesize(&k, EPOCH, MBA, [0.0; 3], &main_belt_times(EPOCH - 20.0));
    let r = orbfit::residuals(&a, EPOCH, &MBA, &[0.0; 3], &k, &FitOptions::default(), true).unwrap();
    let max = r.resid.iter().fold(0.0f64, |m, v| m.max(v.abs()));
    assert!(max < 1e-10, "max residual {:.3e} rad", max);

    // partials against central differences, relative to each column's largest entry
    let h = [1e-6, 1e-6, 1e-6, 1e-8, 1e-8, 1e-8];
    let opts = FitOptions::default();
    for j in 0..6 {
        let (mut sp, mut sm) = (MBA, MBA);
        sp[j] += h[j];
        sm[j] -= h[j];
        let rp = orbfit::residuals(&a, EPOCH, &sp, &[0.0; 3], &k, &opts, false).unwrap();
        let rm = orbfit::residuals(&a, EPOCH, &sm, &[0.0; 3], &k, &opts, false).unwrap();
        let scale = (0..r.resid.len()).fold(0.0f64, |m, row| m.max(r.jacobian[row * 6 + j].abs()));
        for row in 0..r.resid.len() {
            let fd = (rp.resid[row] - rm.resid[row]) / (2.0 * h[j]);
            let an = r.jacobian[row * 6 + j];
            assert!((fd - an).abs() <= 1e-5 * scale, "d r[{}]/d x[{}]: {:.6e} vs {:.6e}", row, j, an, fd);
        }
    }
}

#[test]
fn determines_a_main_belt_orbit_from_scratch() {
    let Some(k) = kernel() else { return };
    let a = synthesize(&k, EPOCH, MBA, [0.0; 3], &main_belt_times(EPOCH - 20.0));
    let fit = determine_orbit(&a, None, &k, &FitOptions::default());
    assert_eq!(fit.flag, FitFlag::Converged, "{:?}", fit.flag);
    assert!(fit.chi2 < 1e-6, "chi2 {}", fit.chi2);
    assert_eq!(fit.ndof, 2 * a.len() as i64 - 6);
    // compare with the truth at the fit's epoch
    let t0 = Time::new(EPOCH, "tdb", "jd").unwrap();
    let mut sim = SpiceSimulation::horizons(&t0, &k).unwrap();
    sim.integrator.set_epsilon(1e-12);
    sim.add(SpaceRock::from_xyz("rock", MBA[0], MBA[1], MBA[2], MBA[3], MBA[4], MBA[5], t0, "J2000", "SSB").unwrap()).unwrap();
    sim.integrate_jd(fit.epoch, &k).unwrap();
    let p = &sim.state.particles[0];
    let dpos = (p.position - Vector3::new(fit.state[0], fit.state[1], fit.state[2])).norm();
    assert!(dpos < 1e-8, "position error {:.3e} AU", dpos);
    // the formal uncertainty is sensible for 1" astrometry
    let sigma_x = fit.state_covariance()[0][0].sqrt();
    assert!(sigma_x > 1e-9 && sigma_x < 1e-4, "sigma x {:.3e}", sigma_x);
    assert_eq!(fit.residuals.len(), a.len());
    assert!(fit.residuals.iter().all(|r| r[0].is_finite() && r[1].is_finite() && r[2].is_nan() && r[4].is_nan()));
}

#[test]
fn recovers_a2_jointly_with_the_state() {
    let Some(k) = kernel() else { return };
    let a2 = -5e-13;
    let times: Vec<f64> = (0..40).map(|i| EPOCH + i as f64 * 36.5).collect();
    let neo = [-0.6886, 0.7853, 0.2744, -0.012364, -0.007927, -0.003264];
    let a = synthesize(&k, EPOCH, neo, [0.0, a2, 0.0], &times).with_sigma(&[0.1 * orbfit::ARCSEC; 40], &[0.1 * orbfit::ARCSEC; 40]).unwrap();

    let mut start = neo;
    start[0] += 1e-7;
    start[3] += 1e-9;
    let gravity = fit_orbit(&a, EPOCH, &start, &[0.0; 3], &k, &FitOptions::default(), 100);
    let opts = FitOptions { fit_nongrav: [false, true, false], ..FitOptions::default() };
    let fit = determine_orbit(&a, Some((EPOCH, start, [0.0; 3])), &k, &opts);
    assert_eq!(fit.flag, FitFlag::Converged);
    assert_eq!(fit.npar, 7);
    assert!((fit.nongrav[1] - a2).abs() < 1e-2 * a2.abs(), "A2 {:.4e}", fit.nongrav[1]);
    assert!(fit.nongrav_sigma()[1] > 0.0 && fit.nongrav_sigma()[0].is_nan());
    assert!(fit.chi2 < 1e-3 * gravity.chi2, "chi2 {} vs gravity-only {}", fit.chi2, gravity.chi2);
}

#[test]
fn too_few_detections_are_not_attempted() {
    let Some(k) = kernel() else { return };
    let a = synthesize(&k, EPOCH, MBA, [0.0; 3], &[EPOCH, EPOCH + 1.0]);
    assert_eq!(determine_orbit(&a, None, &k, &FitOptions::default()).flag, FitFlag::NotAttempted);
}

const TNO: [f64; 6] = [35.0, 20.0, 5.0, -0.00135, 0.00236, 0.0002];

fn nights(start: f64, nights: &[f64], per_night: usize) -> Vec<f64> {
    nights.iter().flat_map(|n| (0..per_night).map(move |i| start + n + 0.03 * i as f64)).collect()
}

#[test]
fn bk_iod_recovers_distance_and_motion_of_a_distant_object() {
    let Some(k) = kernel() else { return };
    let times = nights(EPOCH - 1.0, &[0.0, 1.0, 3.0], 3);
    let a = synthesize(&k, EPOCH, TNO, [0.0; 3], &times);
    let mid = a.epoch[a.len() / 2];
    let s = orbfit::bk_iod(&a, mid).unwrap();
    let r_true = (TNO[0] * TNO[0] + TNO[1] * TNO[1] + TNO[2] * TNO[2]).sqrt();
    let r = (s[0] * s[0] + s[1] * s[1] + s[2] * s[2]).sqrt();
    assert!((r - r_true).abs() < 0.01 * r_true, "r {} vs {}", r, r_true);
    // the sky-plane velocity is recovered too (the radial part is pinned to zero)
    let v_true = Vector3::new(TNO[3], TNO[4], TNO[5]);
    let n = Vector3::new(s[0], s[1], s[2]).normalize();
    let v = Vector3::new(s[3], s[4], s[5]);
    let vt = v_true - n * n.dot(&v_true);
    assert!((v - vt).norm() < 0.05 * vt.norm(), "tangential velocity {:?} vs {:?}", v, vt);
    // too few detections
    assert!(orbfit::bk_iod(&a.subset(&[0, 1]), mid).is_none());
}

/// A TNO for which Gauss's method finds no real root on two nights from X05.
const TNO_NO_GAUSS: [f64; 6] = [
    -9.225475943135795, -39.13592432366724, 3.8411281325542905,
    0.0024400553922535283, -0.0007506509951382269, 0.0005231263829549057,
];

#[test]
fn bk_fallback_rescues_a_two_night_tno_arc() {
    let Some(k) = kernel() else { return };
    let a = synthesize(&k, EPOCH, TNO_NO_GAUSS, [0.0; 3], &nights(EPOCH, &[0.0, 1.0], 3));
    let without = determine_orbit(&a, None, &k, &FitOptions { bk_fallback: false, ..FitOptions::default() });
    assert_eq!(without.flag, FitFlag::NoSolution);
    let with = determine_orbit(&a, None, &k, &FitOptions::default());
    assert_eq!(with.flag, FitFlag::Converged, "{:?}", with.flag);
    let r_true = TNO_NO_GAUSS[..3].iter().map(|v| v * v).sum::<f64>().sqrt();
    let r = with.state[..3].iter().map(|v| v * v).sum::<f64>().sqrt();
    assert!((r - r_true).abs() < 0.01 * r_true, "r {} vs {}", r, r_true);
}

// ---- rates and radar, against layup's independently generated fixtures ----------------------

fn fixture(name: &str) -> serde_json::Value {
    let p = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data").join(name);
    serde_json::from_str(&std::fs::read_to_string(p).unwrap()).unwrap()
}

fn v3(v: &serde_json::Value) -> [f64; 3] {
    [v[0].as_f64().unwrap(), v[1].as_f64().unwrap(), v[2].as_f64().unwrap()]
}

fn state6(v: &serde_json::Value) -> [f64; 6] {
    let mut s = [0.0; 6];
    for (k, x) in s.iter_mut().enumerate() {
        *x = v[k].as_f64().unwrap();
    }
    s
}

fn streak_fixture() -> (Astrometry, f64, [f64; 6]) {
    let d = fixture("layup_streak_synthetic.json");
    let obs = d["observations"].as_array().unwrap();
    let f = |k: &str| obs.iter().map(|o| o[k].as_f64().unwrap()).collect::<Vec<f64>>();
    let n = obs.len();
    let (su, sr) = (d["ra_unc_rad"].as_f64().unwrap(), d["rate_unc_radday"].as_f64().unwrap());
    let a = Astrometry::new(f("epoch"), f("ra"), f("dec"), obs.iter().map(|o| v3(&o["observer_position"])).collect())
        .unwrap()
        .with_sigma(&vec![su; n], &vec![su; n])
        .unwrap()
        .with_observer_velocity(obs.iter().map(|o| v3(&o["observer_velocity"])).collect())
        .unwrap()
        .with_rates(f("ra_rate"), f("dec_rate"), &vec![sr; n], &vec![sr; n])
        .unwrap();
    (a, d["epoch"].as_f64().unwrap(), state6(&d["true_state"]))
}

fn radar_fixture(delay: bool, doppler: bool) -> (Astrometry, f64, [f64; 6]) {
    let d = fixture("layup_radar_synthetic.json");
    let obs = d["observations"].as_array().unwrap();
    let f = |k: &str, on: bool| obs.iter().map(|o| if on { o[k].as_f64().unwrap() } else { f64::NAN }).collect::<Vec<f64>>();
    let n = obs.len();
    let a = Astrometry::new(f("epoch", true), vec![f64::NAN; n], vec![f64::NAN; n], obs.iter().map(|o| v3(&o["observer_position"])).collect())
        .unwrap()
        .with_observer_velocity(obs.iter().map(|o| v3(&o["observer_velocity"])).collect())
        .unwrap()
        .with_radar(f("delay", delay), f("doppler", doppler), &vec![d["delay_unc_days"].as_f64().unwrap(); n], &vec![d["doppler_unc_audy"].as_f64().unwrap(); n])
        .unwrap();
    (a, d["epoch"].as_f64().unwrap(), state6(&d["true_state"]))
}

/// Residuals at the true orbit, in units of each row's sigma.
fn normalized_residuals(k: &SpiceKernel, a: &Astrometry, epoch: f64, s: &[f64; 6]) -> Vec<f64> {
    let r = orbfit::residuals(a, epoch, s, &[0.0; 3], k, &FitOptions::default(), false).unwrap();
    r.resid.iter().zip(r.row_detection.iter().zip(&r.row_kind)).map(|(v, (&i, &kind))| v / a.sigma(i, kind)).collect()
}

/// Every partial against central differences, relative to its column's largest entry.
///
/// Rows of optical astrometry agree to the differencing noise. The rate and radar partials are
/// layup's approximations: the rate Jacobian leaves out the light-time-rate factor
/// (1 - d(light time)/dt) and its derivative, an error of order |v|/c times the rate's velocity
/// partials; the radar Jacobian treats the up and down legs as parallel (they differ by
/// ~|observer displacement| / distance). So those rows are checked against the largest partial of
/// the same kind over all parameters. The residuals themselves are exact, so fits still converge
/// to the right orbit; only the step direction is slightly off.
fn check_partials(k: &SpiceKernel, a: &Astrometry, epoch: f64, s: &[f64; 6]) {
    let opts = FitOptions::default();
    let r = orbfit::residuals(a, epoch, s, &[0.0; 3], k, &opts, true).unwrap();
    let h = [1e-6, 1e-6, 1e-6, 1e-8, 1e-8, 1e-8];
    for j in 0..6 {
        let (mut sp, mut sm) = (*s, *s);
        sp[j] += h[j];
        sm[j] -= h[j];
        let rp = orbfit::residuals(a, epoch, &sp, &[0.0; 3], k, &opts, false).unwrap();
        let rm = orbfit::residuals(a, epoch, &sm, &[0.0; 3], k, &opts, false).unwrap();
        for kind in [orbfit::RowKind::Ra, orbfit::RowKind::RaRate, orbfit::RowKind::DecRate, orbfit::RowKind::Delay, orbfit::RowKind::Doppler] {
            let rows: Vec<usize> = (0..r.resid.len()).filter(|&q| r.row_kind[q] == kind).collect();
            let cols: Vec<usize> = if kind == orbfit::RowKind::Ra { vec![j] } else { (0..6).collect() };
            let scale = rows.iter().fold(0.0f64, |m, &q| cols.iter().fold(m, |m, &c| m.max(r.jacobian[q * 6 + c].abs())));
            for &q in &rows {
                let fd = (rp.resid[q] - rm.resid[q]) / (2.0 * h[j]);
                let an = r.jacobian[q * 6 + j];
                let tol = if kind == orbfit::RowKind::Ra { 1e-4 } else { 1e-3 };
                assert!((fd - an).abs() <= tol * scale, "{:?} row {} param {}: {:.6e} vs {:.6e}", kind, q, j, an, fd);
            }
        }
    }
}

#[test]
fn rate_residuals_match_layups_synthetic_streaks() {
    let Some(k) = kernel() else { return };
    let (a, epoch, truth) = streak_fixture();
    assert_eq!(a.n_rows(), 4 * a.len());
    let r = normalized_residuals(&k, &a, epoch, &truth);
    let max = r.iter().fold(0.0f64, |m, v| m.max(v.abs()));
    assert!(max < 1e-4, "largest residual at the truth {:.2e} sigma", max);
    check_partials(&k, &a, epoch, &truth);

    let mut start = truth;
    for v in &mut start[3..] {
        *v *= 1.001;
    }
    let fit = determine_orbit(&a, Some((epoch, start, [0.0; 3])), &k, &FitOptions::default());
    assert_eq!(fit.flag, FitFlag::Converged);
    assert_eq!(fit.ndof, 4 * a.len() as i64 - 6);
    let dpos = (0..3).map(|i| (fit.state[i] - truth[i]).powi(2)).sum::<f64>().sqrt();
    assert!(dpos < 1e-8, "position error {:.2e} AU", dpos);
    assert!(fit.residuals.iter().all(|r| r[2].is_finite() && r[3].is_finite() && r[4].is_nan()));

    // a corrupted rate shows up in chi-square
    let mut bad = a.clone();
    bad.ra_rate[0] += 50.0 * bad.sigma_ra_rate[0];
    let fit = determine_orbit(&bad, Some((epoch, truth, [0.0; 3])), &k, &FitOptions::default());
    assert!(fit.chi2 > 2000.0, "chi2 {}", fit.chi2);
}

#[test]
fn radar_residuals_match_layups_synthetic_radar() {
    let Some(k) = kernel() else { return };
    for (delay, doppler) in [(true, true), (true, false), (false, true)] {
        let (a, epoch, truth) = radar_fixture(delay, doppler);
        assert_eq!(a.n_rows(), a.len() * (delay as usize + doppler as usize));
        let r = normalized_residuals(&k, &a, epoch, &truth);
        let max = r.iter().fold(0.0f64, |m, v| m.max(v.abs()));
        assert!(max < 2e-3, "largest residual at the truth {:.2e} sigma (delay {}, doppler {})", max, delay, doppler);
        check_partials(&k, &a, epoch, &truth);

        // layup's test: velocity 0.1% off, recovered to 1e-7 AU and 1e-9 AU/day.
        let mut start = truth;
        for v in &mut start[3..] {
            *v *= 1.001;
        }
        let fit = determine_orbit(&a, Some((epoch, start, [0.0; 3])), &k, &FitOptions::default());
        assert_eq!(fit.flag, FitFlag::Converged, "delay {}, doppler {}", delay, doppler);
        if delay && doppler {
            let dpos = (0..3).map(|i| (fit.state[i] - truth[i]).powi(2)).sum::<f64>().sqrt();
            let dvel = (3..6).map(|i| (fit.state[i] - truth[i]).powi(2)).sum::<f64>().sqrt();
            assert!(dpos < 1e-7 && dvel < 1e-9, "error {:.2e} AU, {:.2e} AU/day", dpos, dvel);
        }
    }
}

#[test]
fn explicit_transmitter_agrees_with_monostatic_extrapolation() {
    let Some(k) = kernel() else { return };
    let (a, epoch, truth) = radar_fixture(true, true);
    let mono = orbfit::residuals(&a, epoch, &truth, &[0.0; 3], &k, &FitOptions::default(), false).unwrap();
    // The receiving station, moved back by the round-trip delay along its velocity (what the
    // monostatic model does with zero acceleration), given as a separate transmitter.
    let tx: Vec<[f64; 6]> = (0..a.len())
        .map(|i| {
            let (p, v, tau) = (a.observer[i], a.observer_velocity[i], a.delay[i]);
            [p[0] - v[0] * tau, p[1] - v[1] * tau, p[2] - v[2] * tau, v[0], v[1], v[2]]
        })
        .collect();
    let bi = a.clone().with_transmitter(tx).unwrap();
    let bist = orbfit::residuals(&bi, epoch, &truth, &[0.0; 3], &k, &FitOptions::default(), false).unwrap();
    for q in 0..mono.resid.len() {
        let s = a.sigma(mono.row_detection[q], mono.row_kind[q]);
        assert!((mono.resid[q] - bist.resid[q]).abs() < 0.05 * s, "row {}: {:.3e} vs {:.3e}", q, mono.resid[q], bist.resid[q]);
    }
}

#[test]
fn earth_fixed_state_matches_observatory_at() {
    let Some(k) = kernel() else { return };
    for code in ["251", "253", "X05"] {
        let site = Observatory::from_obscode(code).unwrap();
        let p = site.earth_fixed_position().unwrap();
        for jd in [2456323.5, 2460400.123] {
            let o = site.at(&Time::new(jd, "tdb", "jd").unwrap(), "J2000", "SSB", &k).unwrap();
            let s = spacerocks::observing::earth_fixed_state(&p, jd, &k).unwrap();
            let v = o.velocity.unwrap();
            for c in 0..3 {
                assert!((s[c] - o.position[c]).abs() < 1e-14, "{} position", code);
                assert!((s[c + 3] - v[c]).abs() < 1e-14, "{} velocity", code);
            }
        }
    }
}

// ---- Robust mode ----------------------------------------------------------------------------

/// The fit's state propagated to the truth's epoch, as a position error (AU).
fn position_error(k: &SpiceKernel, fit: &orbfit::OrbitFit, epoch: f64, truth: [f64; 6]) -> f64 {
    let t0 = Time::new(epoch, "tdb", "jd").unwrap();
    let mut sim = SpiceSimulation::horizons(&t0, k).unwrap();
    sim.integrator.set_epsilon(1e-12);
    sim.add(SpaceRock::from_xyz("rock", truth[0], truth[1], truth[2], truth[3], truth[4], truth[5], t0, "J2000", "SSB").unwrap()).unwrap();
    sim.integrate_jd(fit.epoch, k).unwrap();
    let p = &sim.state.particles[0];
    (p.position - Vector3::new(fit.state[0], fit.state[1], fit.state[2])).norm()
}

#[test]
fn robust_mode_rejects_outliers_and_recovers_the_orbit() {
    let Some(k) = kernel() else { return };
    let mut a = synthesize(&k, EPOCH, MBA, [0.0; 3], &main_belt_times(EPOCH - 20.0));
    let bad = [3, 10, 17];
    for &i in &bad {
        a.ra[i] += 50.0 * orbfit::ARCSEC / a.dec[i].cos();
    }
    let layup = determine_orbit(&a, None, &k, &FitOptions::default());
    assert_ne!(layup.flag, FitFlag::Converged, "the outliers should spoil the plain fit");

    let fit = determine_orbit(&a, None, &k, &FitOptions { robust: true, ..FitOptions::default() });
    assert_eq!(fit.flag, FitFlag::Converged, "{:?}", fit.flag);
    let rejected: Vec<usize> = (0..a.len()).filter(|&i| !fit.used[i]).collect();
    assert_eq!(rejected, bad);
    assert_eq!(fit.ndof, 2 * (a.len() - bad.len()) as i64 - 6);
    assert!(fit.chi2 < 1e-6, "chi2 {}", fit.chi2);
    // residuals cover every detection, the rejected ones included
    assert_eq!(fit.residuals.len(), a.len());
    for &i in &bad {
        assert!((fit.residuals[i][0] / orbfit::ARCSEC - 50.0).abs() < 1e-3, "outlier residual {}", fit.residuals[i][0] / orbfit::ARCSEC);
    }
    let err = position_error(&k, &fit, EPOCH, MBA);
    assert!(err < 1e-8, "position error {:.3e} AU", err);
}

#[test]
fn robust_mode_leaves_clean_fits_unchanged() {
    let Some(k) = kernel() else { return };
    let a = synthesize(&k, EPOCH, MBA, [0.0; 3], &main_belt_times(EPOCH - 20.0));
    let layup = determine_orbit(&a, None, &k, &FitOptions::default());
    let fit = determine_orbit(&a, None, &k, &FitOptions { robust: true, ..FitOptions::default() });
    assert_eq!(fit.flag, FitFlag::Converged);
    assert!(fit.used.iter().all(|&u| u));
    assert_eq!(layup.used, fit.used);
    assert_eq!(fit.epoch, layup.epoch);
    for j in 0..6 {
        assert!((fit.state[j] - layup.state[j]).abs() < 1e-12, "component {}", j);
    }
}

/// Bennu's MPC astrometry from 2010 on (tests/data). Its 313-day 2011-12 arc defeats layup's
/// pipeline, which fits each Gauss orbit to that whole arc at once.
fn bennu(k: &SpiceKernel) -> Astrometry {
    let text = std::fs::read_to_string(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/bennu_mpc_2010_2025.csv")).unwrap();
    let mut a = Astrometry::default();
    let mut stn = Vec::new();
    for line in text.lines().skip(1) {
        let c: Vec<&str> = line.split(',').collect();
        let v = |i: usize| c[i].parse::<f64>().unwrap();
        a.epoch.push(v(0));
        a.ra.push(v(1));
        a.dec.push(v(2));
        a.sigma_ra.push(v(3));
        a.sigma_dec.push(v(4));
        stn.push(c[5].to_string());
    }
    a.observer = orbfit::observer_positions(&stn, &a.epoch, k).unwrap();
    a
}

#[test]
fn robust_mode_fits_bennu_where_layups_pipeline_fails() {
    let Some(k) = kernel() else { return };
    let a = bennu(&k);
    let layup = determine_orbit(&a, None, &k, &FitOptions::default());
    assert_eq!(layup.flag, FitFlag::NoRootConverged);
    assert!(layup.used.is_empty() && layup.residuals.len() == a.len());

    let fit = determine_orbit(&a, None, &k, &FitOptions { robust: true, ..FitOptions::default() });
    assert_eq!(fit.flag, FitFlag::Converged, "{:?}", fit.flag);
    let n_used = fit.used.iter().filter(|&&u| u).count();
    assert!(n_used >= a.len() - 5, "used {} of {}", n_used, a.len());
    assert!(fit.chi2 / (fit.ndof as f64) < 1.0);
    // heliocentric ecliptic elements, as published for Bennu: a 1.126 AU, e 0.204, i 6.03 deg
    let t = Time::new(fit.epoch, "tdb", "jd").unwrap();
    let s = fit.state;
    let mut rock = SpaceRock::from_xyz("bennu", s[0], s[1], s[2], s[3], s[4], s[5], t, "J2000", "SSB").unwrap();
    rock.to_helio(&k).unwrap();
    rock.change_reference_plane("ECLIPJ2000").unwrap();
    assert!((rock.a() - 1.126).abs() < 1e-3, "a {}", rock.a());
    assert!((rock.e() - 0.2037).abs() < 1e-3, "e {}", rock.e());
    assert!((rock.inc().to_degrees() - 6.035).abs() < 1e-2, "i {}", rock.inc().to_degrees());
}

#[test]
fn robust_mode_reports_failure_honestly() {
    let Some(k) = kernel() else { return };
    // Bennu from 2015 on: a few scattered nights, which no method can fit.
    let full = bennu(&k);
    let idx: Vec<usize> = (0..full.len()).filter(|&i| full.epoch[i] > 2457023.5).collect();
    let a = full.subset(&idx);
    let fit = determine_orbit(&a, None, &k, &FitOptions { robust: true, ..FitOptions::default() });
    assert!(!matches!(fit.flag, FitFlag::Converged | FitFlag::Chi2TooLarge), "{:?}", fit.flag);
    assert!(fit.used.is_empty());
}

// ---- Automatic non-gravitational model selection ----------------------------------------------

#[test]
fn nongrav_auto_keeps_gravity_when_it_fits() {
    let Some(k) = kernel() else { return };
    let a = synthesize(&k, EPOCH, MBA, [0.0; 3], &main_belt_times(EPOCH - 20.0));
    let opts = FitOptions { nongrav_auto: Some(orbfit::NongravAuto::default()), ..FitOptions::default() };
    let fit = determine_orbit(&a, None, &k, &opts);
    assert_eq!(fit.flag, FitFlag::Converged);
    assert_eq!(fit.npar, 6);
    assert_eq!(fit.fit_nongrav, [false; 3]);
}

#[test]
fn nongrav_auto_picks_a2_for_a_yarkovsky_drifter() {
    let Some(k) = kernel() else { return };
    // Strong enough that the gravity-only fit fails layup's reduced chi-square gate (1.5).
    let a2 = -1e-11;
    let times: Vec<f64> = (0..40).map(|i| EPOCH + i as f64 * 36.5).collect();
    let neo = [-0.6886, 0.7853, 0.2744, -0.012364, -0.007927, -0.003264];
    let a = synthesize(&k, EPOCH, neo, [0.0, a2, 0.0], &times).with_sigma(&[0.1 * orbfit::ARCSEC; 40], &[0.1 * orbfit::ARCSEC; 40]).unwrap();
    let mut start = neo;
    start[0] += 1e-7;
    start[3] += 1e-9;
    let opts = FitOptions { nongrav_auto: Some(orbfit::NongravAuto::default()), ..FitOptions::default() };
    let fit = determine_orbit(&a, Some((EPOCH, start, [0.0; 3])), &k, &opts);
    assert_eq!(fit.flag, FitFlag::Converged);
    assert_eq!(fit.fit_nongrav, [false, true, false]);
    assert!((fit.nongrav[1] - a2).abs() < 1e-2 * a2.abs(), "A2 {:.4e}", fit.nongrav[1]);

    // The same choice, made directly from the gravity-only fit.
    let gravity = fit_orbit(&a, EPOCH, &start, &[0.0; 3], &k, &FitOptions::default(), 100);
    assert!(gravity.chi2 / gravity.ndof as f64 > 1.5, "gravity-only reduced chi2 {}", gravity.chi2 / gravity.ndof as f64);
    let chosen = orbfit::select_nongrav(&a, &gravity, &k, &FitOptions::default(), &orbfit::NongravAuto::default());
    assert_eq!(chosen.fit_nongrav, [false, true, false]);
    // Raising the bar past the A2 detection keeps the gravity-only orbit.
    let strict = orbfit::NongravAuto { nsigma: 1e9, ..orbfit::NongravAuto::default() };
    let kept = orbfit::select_nongrav(&a, &gravity, &k, &FitOptions::default(), &strict);
    assert_eq!(kept.npar, 6);
}

// ---- Scaled convergence (conv_frac) -------------------------------------------------------------

#[test]
fn conv_frac_stops_within_a_fraction_of_sigma() {
    let Some(k) = kernel() else { return };
    let mut a = synthesize(&k, EPOCH, MBA, [0.0; 3], &main_belt_times(EPOCH - 20.0));
    // Deterministic 0.5" "noise", so the formal uncertainties are meaningful.
    for i in 0..a.len() {
        a.ra[i] += 0.5 * orbfit::ARCSEC * ((i as f64) * 1.7).sin() / a.dec[i].cos();
        a.dec[i] += 0.5 * orbfit::ARCSEC * ((i as f64) * 2.3).cos();
    }
    let mut start = MBA;
    start[0] += 1e-4;
    start[4] += 1e-6;
    let absolute = fit_orbit(&a, EPOCH, &start, &[0.0; 3], &k, &FitOptions::default(), 100);
    let scaled = fit_orbit(&a, EPOCH, &start, &[0.0; 3], &k, &FitOptions { conv_frac: 0.1, ..FitOptions::default() }, 100);
    assert_eq!(absolute.flag, FitFlag::Converged);
    assert_eq!(scaled.flag, FitFlag::Converged);
    assert!(scaled.niter <= absolute.niter, "{} vs {}", scaled.niter, absolute.niter);
    let cov = absolute.state_covariance();
    for j in 0..6 {
        let d = (scaled.state[j] - absolute.state[j]).abs() / cov[j][j].sqrt();
        assert!(d < 0.1, "component {}: {} sigma", j, d);
    }
    // conv_frac = 0 is the absolute test exactly.
    let zero = fit_orbit(&a, EPOCH, &start, &[0.0; 3], &k, &FitOptions { conv_frac: 0.0, ..FitOptions::default() }, 100);
    assert_eq!(zero.state, absolute.state);
    assert_eq!(zero.niter, absolute.niter);
}

// ---- Vereš et al. (2017) uncertainties ---------------------------------------------------------

#[test]
fn veres_sigma_follows_layups_table() {
    use orbfit::veres_sigma;
    assert_eq!(veres_sigma("703", 2456658.5, None, None), 1.0); // on the date split: early
    assert_eq!(veres_sigma("703", 2456658.6, None, None), 0.8);
    assert_eq!(veres_sigma("F51", 2460000.5, None, None), 0.2);
    assert_eq!(veres_sigma("W84", 2460000.5, None, None), 0.5); // layup's first W84 branch wins
    assert_eq!(veres_sigma("568", 2460000.5, Some("Gaia3E"), None), 0.1);
    assert_eq!(veres_sigma("568", 2460000.5, Some("X"), None), 0.1); // MPC code for Gaia EDR3
    assert_eq!(veres_sigma("G83", 2460000.5, Some("UCAC4"), Some("2")), 0.3);
    assert_eq!(veres_sigma("G83", 2460000.5, Some("UCAC4"), None), 1.0); // generic, with catalog
    assert_eq!(veres_sigma("309", 2460000.5, Some("UCAC2"), Some("&")), 1.5);
    assert_eq!(veres_sigma("X05", 2460000.5, None, None), 1.5);
    assert_eq!(veres_sigma("X05", 2460000.5, Some(""), None), 1.5); // empty = absent
    assert_eq!(veres_sigma("X05", 2460000.5, Some("Gaia2"), None), 1.0);
}

// ---- Star-catalog debiasing ---------------------------------------------------------------------

#[test]
fn healpix_nested_pixels_match_healpy() {
    use orbfit::debias::ang2pix_nest_lonlat;
    // healpy.ang2pix(256, lon, lat, nest=True, lonlat=True)
    for (lon, lat, pix) in [
        (0.0, 0.0, 311296),
        (10.0, 20.0, 317814),
        (359.9, -45.0, 742485),
        (123.4, 89.99, 131071),
        (200.0, -89.995, 655360),
        (45.0, 41.8103148957786, 49152),
        (300.0, -30.0, 767385),
    ] {
        assert_eq!(ang2pix_nest_lonlat(256, lon, lat), pix, "({}, {})", lon, lat);
    }
}

#[test]
fn debias_applies_offsets_and_proper_motion() {
    // A one-catalog-per-column synthetic table at nside 1 (12 pixels): every pixel has UCAC4
    // offsets of +0.2" (RA cos Dec) and -0.1" (Dec) and proper motions of 10 and -20 mas/yr.
    let dir = std::env::temp_dir().join(format!("spacerocks_debias_test_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let dat = dir.join("bias.dat");
    let q = orbfit::debias::BIAS_CATALOGS.iter().position(|(n, _)| *n == "UCAC4").unwrap();
    let mut text = String::new();
    for i in 0..23 {
        text.push_str(&format!("! header {}\n", i));
    }
    for _ in 0..12 {
        let mut row = vec![0.0f64; 104];
        row[4 * q..4 * q + 4].copy_from_slice(&[0.2, -0.1, 10.0, -20.0]);
        text.push_str(&row.iter().map(|v| format!("{:.3}", v)).collect::<Vec<_>>().join(" "));
        text.push('\n');
    }
    std::fs::write(&dat, text).unwrap();
    let table = orbfit::BiasTable::load(&dat).unwrap();
    assert_eq!(table.nside(), 1);
    let (ra, dec) = (1.0f64, 0.5f64);
    let jd = 2451545.0 + 10.0 * 365.25; // 10 Julian years after J2000
    for cat in ["UCAC4", "q"] {
        let (r, d) = table.debias(ra, dec, jd, cat);
        let dra = (r - ra) * dec.cos() / orbfit::ARCSEC;
        let ddec = (d - dec) / orbfit::ARCSEC;
        assert!((dra - -(0.2 + 0.1)).abs() < 1e-6, "dRA cos Dec {}", dra); // 0.2" + 10 yr * 10 mas/yr
        assert!((ddec - -(-0.1 - 0.2)).abs() < 1e-6, "dDec {}", ddec); // -0.1" - 10 yr * 20 mas/yr
    }
    // catalogs without a bias model are left alone
    assert_eq!(table.debias(ra, dec, jd, "Gaia2"), (ra, dec));
    assert_eq!(table.debias(ra, dec, jd, ""), (ra, dec));
    std::fs::remove_dir_all(&dir).ok();
}

// ---- Observer positions: ADES, roving observers, before 1962 -----------------------------------

#[test]
fn geodetic_positions_follow_wgs84() {
    let au = 149597870700.0;
    let p = orbfit::geodetic_to_earth_fixed(0.0, 0.0, 0.0);
    assert!((p[0] * au - 6378137.0).abs() < 1e-6 && p[1].abs() < 1e-15 && p[2].abs() < 1e-15);
    let p = orbfit::geodetic_to_earth_fixed(123.0, 90.0, 0.0);
    assert!((p[2] * au - 6356752.314245).abs() < 1e-3, "polar radius {}", p[2] * au); // WGS84 b
    let p = orbfit::geodetic_to_earth_fixed(90.0, 0.0, 1000.0);
    assert!(p[0].abs() * au < 1e-6 && (p[1] * au - 6379137.0).abs() < 1e-6);
}

#[test]
fn ades_positions_are_geocentric_offsets() {
    let Some(k) = kernel() else { return };
    let jd = 2460400.5;
    let earth = k.state_au(399, 0, jd).unwrap();
    let km = 1.0 / 149597870.7;
    let s = orbfit::ades_observer_state("ICRF_KM", 399, [7000.0, -1000.0, 300.0], Some([1.0, 7.0, 0.5]), jd, &k).unwrap();
    for (j, (off, rate)) in [(7000.0, 1.0), (-1000.0, 7.0), (300.0, 0.5)].iter().enumerate() {
        assert!((s[j] - earth[j] - off * km).abs() < 1e-15);
        assert!((s[j + 3] - earth[j + 3] - rate * km * 86400.0).abs() < 1e-15);
    }
    // without a velocity the observer moves with the Earth's center; AU units
    let s = orbfit::ades_observer_state("ICRF_AU", 399, [1e-4, 0.0, 0.0], None, jd, &k).unwrap();
    assert_eq!(&s[3..], &earth[3..]);
    assert!((s[0] - earth[0] - 1e-4).abs() < 1e-15);
    // WGS84 is a point on the rotating Earth, like a ground station
    let s = orbfit::ades_observer_state("WGS84", 399, [10.0, 45.0, 100.0], None, jd, &k).unwrap();
    let t = spacerocks::observing::earth_fixed_state(&orbfit::geodetic_to_earth_fixed(10.0, 45.0, 100.0), jd, &k).unwrap();
    assert_eq!(s, t);
    let rot = ((s[3] - earth[3]).powi(2) + (s[4] - earth[4]).powi(2)).sqrt() * 149597870.7 / 86400.0;
    assert!(rot > 0.3 && rot < 0.35, "rotation speed {} km/s at 45 deg", rot);
    // other frames and centers are refused
    assert!(orbfit::ades_observer_state("ICRF_KM", 10, [7000.0, 0.0, 0.0], None, jd, &k).is_err());
    assert!(orbfit::ades_observer_state("ECLIPJ2000", 399, [7000.0, 0.0, 0.0], None, jd, &k).is_err());
}

#[test]
fn stations_before_1962_use_the_iau_earth_model() {
    let Some(mut k) = kernel() else { return };
    let dir = PathBuf::from(std::env::var("SPACEROCKS_KERNELS").unwrap());
    let jd_1938 = 2429231.5;
    let site = Observatory::from_obscode("024").unwrap();
    let t = Time::new(jd_1938, "tdb", "jd").unwrap();
    if !dir.join("pck00010.tpc").exists() {
        return;
    }
    // Without the IAU model: an error that says what to load.
    let err = site.at(&t, "J2000", "SSB", &k).unwrap_err().to_string();
    assert!(err.contains("pck00010"), "{}", err);
    k.load(dir.join("pck00010.tpc")).unwrap();
    let o = site.at(&t, "J2000", "SSB", &k).unwrap();
    let earth = k.state_au(399, 0, jd_1938).unwrap();
    let d = ((o.position.x - earth[0]).powi(2) + (o.position.y - earth[1]).powi(2) + (o.position.z - earth[2]).powi(2)).sqrt() * 149597870.7;
    assert!(d > 6300.0 && d < 6400.0, "geocentric distance {} km", d);
    // Continuous across the 1962 switch to ITRF93: the IAU model's prime meridian runs on TDB,
    // and corrected for UT1 (Delta T) it matches ITRF93 to ~0.3 km (the nutation it omits).
    // Uncorrected, as in layup, the jump is ~9 km.
    let (a, b) = (spacerocks::spice::ITRF93_START_ET, 1.0);
    let before = Time::new(2451545.0 + (a - b) / 86400.0, "tdb", "jd").unwrap();
    let after = Time::new(2451545.0 + (a + b) / 86400.0, "tdb", "jd").unwrap();
    let (p0, p1) = (site.at(&before, "J2000", "SSB", &k).unwrap(), site.at(&after, "J2000", "SSB", &k).unwrap());
    let jump = (p1.position - p0.position - p0.velocity.unwrap() * (2.0 * b / 86400.0)).norm() * 149597870.7;
    assert!(jump < 0.5, "jump at the switch {} km", jump);
}

#[test]
fn occultation_offsets_apply_on_the_tangent_plane() {
    let as_ = orbfit::ARCSEC;
    // no offset: the star
    assert_eq!(orbfit::occultation_radec(1.0, 0.3, 0.0, 0.0), (1.0, 0.3));
    // small offsets: layup's (intended) linear formula to second order
    for &(ra, dec) in &[(0.1, 0.0), (2.0, 0.7), (5.5, -1.2), (3.0, 1.5)] {
        let (dra, ddec) = (0.8 * as_, -0.5 * as_);
        let (r, d) = orbfit::occultation_radec(ra, dec, dra, ddec);
        let lin_r = ra + dra / f64::cos(dec);
        let lin_d = dec + ddec;
        // they differ at second order, ~ offset^2 tan(dec) / 2
        let bound = 2.0 * (as_ * as_) * dec.tan().abs().max(1.0);
        assert!(((r - lin_r) * dec.cos()).abs() < bound && (d - lin_d).abs() < bound, "({}, {})", ra, dec);
    }
    // RA wraps into [0, 2 pi)
    let (r, _) = orbfit::occultation_radec(std::f64::consts::TAU - 1e-7, 0.0, 2e-7, 0.0);
    assert!(r > 0.0 && r < 2e-7);
}

// ---- Per-arc non-gravitational parameters --------------------------------------------------------

#[test]
fn per_arc_recovers_distinct_amplitudes_before_and_after_the_epoch() {
    let Some(k) = kernel() else { return };
    let neo = [-0.6886, 0.7853, 0.2744, -0.012364, -0.007927, -0.003264];
    let (a2_before, a2_after) = (1e-12, 5e-12);
    let before: Vec<f64> = (1..=20).map(|i| EPOCH - i as f64 * 50.0).collect();
    let after: Vec<f64> = (1..=20).map(|i| EPOCH + i as f64 * 50.0).collect();
    let a = synthesize(&k, EPOCH, neo, [0.0, a2_before, 0.0], &before);
    let b = synthesize(&k, EPOCH, neo, [0.0, a2_after, 0.0], &after);
    let mut all = a.clone();
    all.epoch.extend(&b.epoch);
    all.ra.extend(&b.ra);
    all.dec.extend(&b.dec);
    all.observer.extend(&b.observer);
    all.sigma_ra = vec![0.1 * orbfit::ARCSEC; 40];
    all.sigma_dec = vec![0.1 * orbfit::ARCSEC; 40];

    let opts = FitOptions { fit_nongrav: [false, true, false], nongrav_per_arc: true, ..FitOptions::default() };
    let fit = determine_orbit(&all, Some((EPOCH, neo, [0.0; 3])), &k, &opts);
    assert_eq!(fit.flag, FitFlag::Converged);
    assert!(fit.per_arc);
    assert_eq!(fit.npar, 8);
    assert_eq!(fit.ndof, 80 - 8);
    assert!((fit.nongrav[1] - a2_before).abs() < 0.02 * a2_before, "arc A {:.4e}", fit.nongrav[1]);
    assert!((fit.nongrav_arc2[1] - a2_after).abs() < 0.02 * a2_after, "arc B {:.4e}", fit.nongrav_arc2[1]);
    assert!(fit.nongrav_sigma()[1] > 0.0 && fit.nongrav_arc2_sigma()[1] > 0.0);
    assert!(fit.nongrav_arc2_sigma()[0].is_nan());

    // A shared amplitude can't do both: its chi-square stays far higher.
    let shared = determine_orbit(&all, Some((EPOCH, neo, [0.0; 3])), &k, &FitOptions { nongrav_per_arc: false, ..opts });
    assert!(!shared.per_arc && shared.npar == 7);
    assert!(shared.chi2 > 100.0 * fit.chi2.max(1e-3), "shared {} vs per-arc {}", shared.chi2, fit.chi2);
}

#[test]
fn herget_iod_finds_the_ranges_of_a_main_belt_arc() {
    let Some(k) = kernel() else { return };
    let times = main_belt_times(EPOCH - 20.0);
    let a = synthesize(&k, EPOCH, MBA, [0.0; 3], &times[..16]);
    let idx: Vec<usize> = (0..a.len()).collect();
    let (epoch, s) = orbfit::herget_iod(&a, &idx, &k, &FitOptions::default()).expect("Herget found no orbit");
    assert_eq!(epoch, a.epoch[0]);
    // the truth at the first detection
    let t0 = Time::new(EPOCH, "tdb", "jd").unwrap();
    let mut sim = SpiceSimulation::horizons(&t0, &k).unwrap();
    sim.integrator.set_epsilon(1e-12);
    sim.add(SpaceRock::from_xyz("rock", MBA[0], MBA[1], MBA[2], MBA[3], MBA[4], MBA[5], t0, "J2000", "SSB").unwrap()).unwrap();
    sim.integrate_jd(epoch, &k).unwrap();
    let p = &sim.state.particles[0];
    let dpos = (p.position - Vector3::new(s[0], s[1], s[2])).norm();
    let dvel = (p.velocity - Vector3::new(s[3], s[4], s[5])).norm();
    // An initial orbit: the ranges stop once they move by less than 0.003 AU per iteration, the
    // returned state is the one before that last move, and the velocity is two-body.
    assert!(dpos < 0.03, "position error {:.3e} AU", dpos);
    assert!(dvel < 0.03 * p.velocity.norm(), "velocity error {:.3e} of {:.3e} AU/day", dvel, p.velocity.norm());
}

#[test]
fn herget_seeds_the_full_pipeline() {
    let Some(k) = kernel() else { return };
    let a = synthesize(&k, EPOCH, MBA, [0.0; 3], &main_belt_times(EPOCH - 20.0));
    let fit = determine_orbit(&a, None, &k, &FitOptions { herget: true, ..FitOptions::default() });
    assert_eq!(fit.flag, FitFlag::Converged, "{:?}", fit.flag);
    assert!(fit.chi2 < 1e-6, "chi2 {}", fit.chi2);
    assert_eq!(fit.epoch, a.epoch[0]);
}

#[test]
fn ut1_corrected_iau_earth_follows_itrf93() {
    // The correction that places stations before 1962, checked where both models exist.
    let Some(mut k) = kernel() else { return };
    let dir = PathBuf::from(std::env::var("SPACEROCKS_KERNELS").unwrap());
    if k.load(dir.join("pck00010.tpc")).is_err() {
        return;
    }
    use spacerocks::spice::{frames, iau_earth_rotation_delay};
    let stations = [[6378.137, 0.0, 0.0], [0.0, 4510.0, 4510.0], [-3189.0, -3189.0, 4510.0]];
    let (mut worst, mut worst_raw) = (0.0f64, 0.0f64);
    let mut year = 1962.1;
    while year < 2019.0 {
        let et = (year - 2000.0) * 365.25 * 86400.0;
        let truth = k.frame_to_j2000(frames::ITRF93, et).unwrap().rotation;
        let model = k.frame_to_j2000(frames::IAU_EARTH, et - iau_earth_rotation_delay(et)).unwrap().rotation;
        let raw = k.frame_to_j2000(frames::IAU_EARTH, et).unwrap().rotation;
        for r in stations {
            let apply = |m: &[[f64; 3]; 3]| [0, 1, 2].map(|i| (0..3).map(|j| m[i][j] * r[j]).sum::<f64>());
            let (a, b, c) = (apply(&truth), apply(&model), apply(&raw));
            worst = worst.max((0..3).map(|i| (a[i] - b[i]).powi(2)).sum::<f64>().sqrt());
            worst_raw = worst_raw.max((0..3).map(|i| (a[i] - c[i]).powi(2)).sum::<f64>().sqrt());
        }
        year += 0.37;
    }
    assert!(worst < 0.4, "corrected IAU_EARTH off ITRF93 by {} km", worst);
    assert!(worst_raw > 5.0, "uncorrected {} km", worst_raw);
}

// ---- Sequential updates and incremental fitting (layup issue #419) ------------------------------

#[test]
fn sequential_update_matches_a_full_refit() {
    let Some(k) = kernel() else { return };
    let times = main_belt_times(EPOCH - 20.0);
    let all = synthesize(&k, EPOCH, MBA, [0.0; 3], &times);
    let (n_old, n) = (16, all.len());
    let old_idx: Vec<usize> = (0..n_old).collect();
    let new_idx: Vec<usize> = (n_old..n).collect();
    let opts = FitOptions::default();
    let prior = determine_orbit(&all.subset(&old_idx), None, &k, &opts);
    assert!(prior.converged());

    let (seq, accepted) = orbfit::sequential_update(&all.subset(&new_idx), &prior, None, &k, &opts);
    assert!(accepted && seq.converged(), "{:?}", seq.flag);
    assert_eq!(seq.epoch, prior.epoch);
    assert_eq!(seq.ndof, 2 * (n - n_old) as i64);
    let full = fit_orbit(&all, prior.epoch, &prior.state, &[0.0; 3], &k, &opts, 100);
    let d = orbfit::update_mahalanobis(&full, &seq);
    assert!(d < 1e-3, "sequential vs full refit: {} sigma", d);
    // the posterior is tighter than the prior
    assert!(seq.state_covariance()[0][0] < prior.state_covariance()[0][0]);
}

#[test]
fn sequential_update_flags_what_it_cannot_do() {
    let Some(k) = kernel() else { return };
    let times = main_belt_times(EPOCH - 20.0);
    let all = synthesize(&k, EPOCH, MBA, [0.0; 3], &times);
    let old: Vec<usize> = (0..16).collect();
    let new: Vec<usize> = (16..all.len()).collect();
    let opts = FitOptions::default();
    let prior = determine_orbit(&all.subset(&old), None, &k, &opts);
    // a move beyond the gate (a prior nudged by one sigma, and the gate at 0.01 sigma)
    let gate = FitOptions { max_update_sigma: 0.01, ..opts.clone() };
    let mut prior = prior;
    prior.state[0] += prior.state_covariance()[0][0].sqrt();
    let (f, accepted) = orbfit::sequential_update(&all.subset(&new), &prior, None, &k, &gate);
    assert!(!accepted);
    assert_eq!(f.flag, FitFlag::NonlinearUpdate);
    // ... and with all detections it refits instead
    let (f, accepted) = orbfit::sequential_update(&all.subset(&new), &prior, Some(&all), &k, &gate);
    assert!(!accepted && f.converged() && f.ndof == 2 * all.len() as i64 - 6);
    // a covariance that is not positive definite
    let mut bad = prior.clone();
    bad.covariance[0] = -1.0;
    let (f, _) = orbfit::sequential_update(&all.subset(&new), &bad, None, &k, &opts);
    assert_eq!(f.flag, FitFlag::PriorNotPositiveDefinite);
}

#[test]
fn update_orbit_routes_like_layups_incremental_fit() {
    use orbfit::{update_orbit, PriorFit, UpdateRoute};
    let Some(k) = kernel() else { return };
    let times = main_belt_times(EPOCH - 20.0);
    let all = synthesize(&k, EPOCH, MBA, [0.0; 3], &times);
    let opts = FitOptions::default();
    let old = all.subset(&(0..16).collect::<Vec<_>>());
    let fit = determine_orbit(&old, None, &k, &opts);
    let prior = PriorFit { fit: fit.clone(), keys: old.detection_keys() };

    // unchanged (in another order): the prior itself
    let shuffled = old.subset(&(0..16).rev().collect::<Vec<_>>());
    let (f, r) = update_orbit(&shuffled, Some(&prior), &k, &opts);
    assert_eq!(r, UpdateRoute::Skip);
    assert_eq!(f.state, fit.state);
    assert_eq!(f.residuals[0][..2], fit.residuals[15][..2]);
    // added detections
    let (f, r) = update_orbit(&all, Some(&prior), &k, &opts);
    assert_eq!(r, UpdateRoute::Sequential);
    assert!(f.converged() && f.residuals.len() == all.len() && f.residuals[0][0].is_nan() && f.residuals[20][0].is_finite());
    // one removed
    let (f, r) = update_orbit(&all.subset(&(1..all.len()).collect::<Vec<_>>()), Some(&prior), &k, &opts);
    assert_eq!(r, UpdateRoute::Full);
    assert!(f.converged());
    // no prior
    assert_eq!(update_orbit(&all, None, &k, &opts).1, UpdateRoute::Cold);
    // the fingerprint ignores order but not content
    assert_eq!(old.fingerprint(), shuffled.fingerprint());
    let mut changed = old.clone();
    changed.ra[3] += 1e-9;
    assert_ne!(old.fingerprint(), changed.fingerprint());
}

// ---- The Bernstein-Khushalani fitting engine (layup's engine="bk_native") ------------------------

#[test]
fn bk_engine_agrees_with_cartesian_on_a_well_observed_orbit() {
    let Some(k) = kernel() else { return };
    let a = synthesize(&k, EPOCH, MBA, [0.0; 3], &main_belt_times(EPOCH - 20.0));
    let cart = determine_orbit(&a, None, &k, &FitOptions::default());
    let bk = determine_orbit(&a, None, &k, &FitOptions { engine: orbfit::Engine::BkNative, ..FitOptions::default() });
    assert!(cart.converged() && bk.converged(), "{:?} {:?}", cart.flag, bk.flag);
    // the bound-orbit prior is weak here: the same orbit to well within its uncertainty
    let d = orbfit::update_mahalanobis(&cart, &bk);
    assert!(d < 1e-3, "{} sigma apart", d);
    assert_eq!(bk.ndof, cart.ndof);
}

#[test]
fn bk_engine_fits_a_short_tno_arc_directly() {
    let Some(k) = kernel() else { return };
    let times = nights(EPOCH - 1.0, &[0.0, 1.0, 3.0], 3);
    let a = synthesize(&k, EPOCH, TNO, [0.0; 3], &times);
    let mid = a.epoch[a.len() / 2];
    let seed = orbfit::bk_iod(&a, mid).unwrap();
    let fit = orbfit::fit_orbit_bk(&a, mid, &seed, &k, &FitOptions::default(), 100);
    assert!(fit.converged(), "{:?}", fit.flag);
    let r_true = (TNO[0] * TNO[0] + TNO[1] * TNO[1] + TNO[2] * TNO[2]).sqrt();
    let r = fit.state[..3].iter().map(|v| v * v).sum::<f64>().sqrt();
    assert!((r - r_true).abs() < 0.01 * r_true, "r {} vs {}", r, r_true);
    // bound, and with a finite covariance
    let v2: f64 = fit.state[3..].iter().map(|v| v * v).sum();
    assert!(v2 < 2.0 * orbfit::bk_fit::MU_SUN / r);
    assert!(fit.state_covariance().iter().flatten().all(|c| c.is_finite()));
}

// ---- Prediction with uncertainty, and comet original/future orbits -----------------------------

#[test]
fn predictions_reproduce_the_detections_and_map_the_covariance() {
    let Some(k) = kernel() else { return };
    let a = synthesize(&k, EPOCH, MBA, [0.0; 3], &main_belt_times(EPOCH - 20.0));
    let fit = determine_orbit(&a, None, &k, &FitOptions::default());
    let opts = FitOptions::default();
    let p = orbfit::predict_astrometry(&fit, &a, &k, &opts).unwrap();
    for (i, q) in p.iter().enumerate() {
        let d = ((q.ra - a.ra[i]) * a.dec[i].cos()).hypot(q.dec - a.dec[i]);
        assert!(d < 1e-9, "detection {}: {} rad", i, d);
    }
    // One parameter's variance maps to the square of the finite-difference sky partial.
    let s = 1e-6;
    let mut one = fit.clone();
    one.covariance = vec![0.0; fit.npar * fit.npar];
    one.covariance[0] = s * s;
    let mut shifted = fit.clone();
    shifted.state[0] += s;
    let (p0, p1) = (orbfit::predict_astrometry(&one, &a, &k, &opts).unwrap(), orbfit::predict_astrometry(&shifted, &a, &k, &opts).unwrap());
    for i in [0, 10, a.len() - 1] {
        let dx = (p1[i].ra - p0[i].ra) * p0[i].dec.cos();
        let dy = p1[i].dec - p0[i].dec;
        let c = p0[i].cov;
        assert!((c[0][0].sqrt() - dx.abs()).abs() < 1e-3 * dx.abs(), "{} vs {}", c[0][0].sqrt(), dx.abs());
        assert!((c[1][1].sqrt() - dy.abs()).abs() < 1e-3 * dy.abs());
        assert!((c[0][1] - dx * dy).abs() < 2e-3 * (dx * dy).abs());
        let (maj, min, _) = p0[i].ellipse();
        assert!((maj - dx.hypot(dy)).abs() < 1e-3 * maj && min < 1e-3 * maj);
    }
}

#[test]
fn comet_original_and_future_orbits() {
    use orbfit::comet::{barycentric_elements, comet_orbit, REFERENCE_DISTANCE};
    let Some(k) = kernel() else { return };
    // A near-parabolic comet at 3 AU inbound (barycentric J2000), q ~ 2 AU, inclined.
    let r: Vector3<f64> = Vector3::new(1.2, -2.5, 1.1);
    let vesc: f64 = (2.0 * orbfit::MU_TOTAL / r.norm()).sqrt();
    let dir = (Vector3::new(-0.2, 0.9, 0.35) - r.normalize() * 0.55).normalize();
    let v = dir * vesc * 0.9995;
    let s = [r.x, r.y, r.z, v.x, v.y, v.z];
    let opts = FitOptions::default();
    let orig = comet_orbit(EPOCH, &s, &[0.0; 3], false, REFERENCE_DISTANCE, &k, &opts).expect("original");
    let fut = comet_orbit(EPOCH, &s, &[0.0; 3], true, REFERENCE_DISTANCE, &k, &opts).expect("future");
    assert!(orig.reached && fut.reached);
    assert!((orig.distance - 250.0).abs() < 1e-3 && (fut.distance - 250.0).abs() < 1e-3);
    assert!(orig.epoch < EPOCH && fut.epoch > EPOCH);
    // inside the planetary region the osculating barycentric 1/a differs by ~1e-3 / AU
    let (inv_a, ..) = barycentric_elements(&s, orbfit::MU_TOTAL);
    assert!((orig.inv_a - inv_a).abs() < 3e-3 && (fut.inv_a - inv_a).abs() < 3e-3 && orig.inv_a != fut.inv_a, "{} {} {}", orig.inv_a, inv_a, fut.inv_a);
    // a main-belt orbit never gets to 250 AU
    assert!(comet_orbit(EPOCH, &MBA, &[0.0; 3], false, REFERENCE_DISTANCE, &k, &opts).is_none());
}
