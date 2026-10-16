//! The orbit-fit dashboard. The page and its controls are tested without kernels; the fits
//! need a planetary ephemeris (`de440.bsp`, `de440s.bsp` or `de421.bsp`) in
//! `SPACEROCKS_KERNELS` and are skipped otherwise (without `sb441-n16.bsp` they leave the
//! asteroids out of the force model).

use std::io::{Read, Write};
use std::net::TcpStream;
use std::path::PathBuf;
use std::time::Duration;

use nalgebra::Vector3;

use spacerocks::assist::SpiceSimulation;
use spacerocks::orbfit::{self, Astrometry, FitDashboard, FitFlag, FitOptions, ForceModel, ASTEROID_PERTURBERS};
use spacerocks::{SpaceRock, SpiceKernel, Time};

const C: f64 = 173.1446326742403;

fn get(url: &str, path: &str) -> (String, String) {
    request(url, "GET", path, "")
}

fn control(url: &str, query: &str) -> (String, String) {
    request(url, "POST", &format!("/control?{query}"), "X-Spacerocks: 1\r\n")
}

fn json(url: &str) -> serde_json::Value {
    serde_json::from_str(&get(url, "/state.json").1).unwrap()
}

fn request(url: &str, method: &str, path: &str, headers: &str) -> (String, String) {
    let addr = url.trim_start_matches("http://").trim_end_matches('/');
    let mut stream = TcpStream::connect(addr).unwrap();
    write!(stream, "{method} {path} HTTP/1.1\r\nHost: {addr}\r\n{headers}\r\n").unwrap();
    let mut bytes = Vec::new();
    stream.read_to_end(&mut bytes).unwrap();
    let response = String::from_utf8_lossy(&bytes);
    let (head, body) = response.split_once("\r\n\r\n").unwrap();
    (head.lines().next().unwrap().to_string(), body.to_string())
}

/// A kernel with a planetary ephemeris, and the force model it supports.
fn kernel() -> Option<(SpiceKernel, ForceModel)> {
    let dir = PathBuf::from(std::env::var("SPACEROCKS_KERNELS").ok()?);
    let mut k = SpiceKernel::new();
    let planets = ["de440.bsp", "de440s.bsp", "de421.bsp"].into_iter().find(|f| dir.join(f).exists())?;
    k.load(dir.join(planets)).ok()?;
    let forces = if k.load(dir.join("sb441-n16.bsp")).is_ok() {
        ForceModel::default()
    } else {
        ForceModel { exclude: ASTEROID_PERTURBERS.to_vec(), ..Default::default() }
    };
    Some((k, forces))
}

/// 40 geocentric detections (four a night, ten nights over 81 days) of a main-belt orbit,
/// with uniform noise of up to one sigma and one detection 30 sigma off in Dec.
fn detections(k: &SpiceKernel, forces: &ForceModel) -> Astrometry {
    let epoch = 2460000.5;
    let state = [1.2, -2.1, -0.3, 0.0085, 0.0052, 0.0011];
    let t0 = Time::new(epoch, "tdb", "jd").unwrap();
    let mut a = Astrometry::default();
    let mut seed = 12345u64;
    let mut noise = || {
        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        ((seed >> 11) as f64 / (1u64 << 53) as f64 - 0.5) * 2.0 * orbfit::DEFAULT_SIGMA
    };
    for i in 0..40 {
        let t = epoch - 30.0 + (i / 4) as f64 * 9.0 + (i % 4) as f64 * 0.02;
        let earth = k.state_au(399, 0, t).unwrap();
        let obs = Vector3::new(earth[0], earth[1], earth[2]);
        let mut sim = SpiceSimulation::horizons_except(&t0, k, &forces.exclude).unwrap();
        sim.integrator.set_epsilon(1e-12);
        sim.add(SpaceRock::from_xyz("rock", state[0], state[1], state[2], state[3], state[4], state[5], t0.clone(), "J2000", "SSB").unwrap()).unwrap();
        let (mut lt, mut rho) = (0.0, Vector3::zeros());
        for _ in 0..6 {
            sim.integrate_rel(t - lt - epoch, k).unwrap();
            rho = sim.state.particles[0].position - obs;
            lt = rho.norm() / C;
        }
        let rho = rho / rho.norm();
        a.epoch.push(t);
        a.ra.push(rho.y.atan2(rho.x).rem_euclid(std::f64::consts::TAU) + noise() / rho.z.asin().cos());
        a.dec.push(rho.z.asin() + noise() + if i == 17 { 30.0 * orbfit::DEFAULT_SIGMA } else { 0.0 });
        a.observer.push([obs.x, obs.y, obs.z]);
        a.sigma_ra.push(orbfit::DEFAULT_SIGMA);
        a.sigma_dec.push(orbfit::DEFAULT_SIGMA);
    }
    a
}

/// Three detections and no kernel: every fit fails, but the page works.
fn kernel_free() -> FitDashboard {
    let a = Astrometry::new(vec![2460000.5, 2460001.5, 2460002.5], vec![1.0, 1.01, 1.02], vec![0.1, 0.1, 0.1], vec![[1.0, 0.0, 0.0]; 3]).unwrap();
    FitDashboard::start(0, "test rock", a, SpiceKernel::new(), FitOptions::default(), None).unwrap()
}

#[test]
fn page_serves_its_state_and_the_shared_assets() {
    let d = kernel_free();
    assert!(d.wait_idle(Duration::from_secs(10)));
    let (status, page) = get(d.url(), "/");
    assert!(status.contains("200") && page.contains("orbit fit"));
    for asset in ["/retro.css", "/retro.js", "/font.ttf"] {
        assert!(get(d.url(), asset).0.contains("200"), "{asset}");
    }
    let s = json(d.url());
    assert_eq!(s["name"], "test rock");
    assert_eq!(s["detections"].as_array().unwrap().len(), 3);
    assert_eq!(s["runs"], 1);
    assert_eq!(s["fit"]["converged"], false);
    assert_eq!(s["settings"]["asteroids"], true);
    assert_eq!(s["settings"]["relativity"], "eih");
    assert!(d.fit().is_some_and(|f| !f.converged()));
}

#[test]
fn controls_need_the_page_header() {
    let d = kernel_free();
    let (status, _) = request(d.url(), "POST", "/control?action=set&key=robust&value=1", "");
    assert!(status.contains("403"));
    assert!(!d.options().robust);
}

#[test]
fn settings_change_the_options_and_refit() {
    let d = kernel_free();
    assert!(d.wait_idle(Duration::from_secs(10)));
    for (key, value) in [
        ("asteroids", "0"),
        ("Jupiter", "0"),
        ("earth_harmonics", "0"),
        ("solar_j2", "0"),
        ("relativity", "simple"),
        ("a2", "1"),
        ("comet", "1"),
        ("engine", "bk"),
        ("iod", "herget"),
        ("robust", "1"),
        ("outlier_sigma", "3.5"),
        ("epsilon", "1e-10"),
        ("max_iter", "50"),
        ("chi2_threshold", "20"),
        ("arc_gap", "30"),
    ] {
        let (status, body) = control(d.url(), &format!("action=set&key={key}&value={value}"));
        assert!(status.contains("200"), "{key}: {body}");
    }
    let o = d.options();
    assert!(ASTEROID_PERTURBERS.iter().all(|c| o.forces.exclude.contains(c)) && o.forces.exclude.contains(&5));
    assert!(!o.forces.earth_harmonics && !o.forces.solar_j2);
    assert_eq!(o.forces.relativity, orbfit::Relativity::Simple);
    assert_eq!(o.fit_nongrav, [false, true, false]);
    assert!(o.gofr.is_some() && o.engine == orbfit::Engine::BkNative && o.herget && o.robust);
    assert_eq!((o.outlier_sigma, o.epsilon, o.max_iter, o.chi2_threshold, o.arc_gap), (3.5, 1e-10, 50, 20.0, 30.0));

    // Auto picks the non-gravitational model, so it clears A2; A2 clears it.
    control(d.url(), "action=set&key=nongrav_auto&value=1");
    assert!(d.options().nongrav_auto.is_some() && d.options().fit_nongrav == [false; 3]);
    control(d.url(), "action=set&key=a1&value=1");
    assert!(d.options().nongrav_auto.is_none() && d.options().fit_nongrav[0]);

    assert!(control(d.url(), "action=set&key=warp&value=1").0.contains("400"));
    assert!(control(d.url(), "action=set&key=epsilon&value=-1").0.contains("400"));
    assert!(d.wait_idle(Duration::from_secs(10)));
    let s = json(d.url());
    assert!(s["runs"].as_u64().unwrap() >= 2, "auto refit: {s}");
    assert!(s["log"][0]["settings"].as_str().unwrap().contains("no asteroids"));
}

#[test]
fn leaving_detections_out_keeps_three() {
    let d = kernel_free();
    assert!(control(d.url(), "action=exclude&index=1&on=1").0.contains("400"));
    let a = Astrometry::new(vec![2460000.5, 2460001.5, 2460002.5, 2460003.5], vec![1.0; 4], vec![0.1; 4], vec![[1.0, 0.0, 0.0]; 4]).unwrap();
    let d = FitDashboard::start(0, "rock", a, SpiceKernel::new(), FitOptions::default(), None).unwrap();
    assert!(control(d.url(), "action=exclude&index=2&on=1").0.contains("200"));
    assert_eq!(d.excluded(), vec![false, false, true, false]);
    assert!(control(d.url(), "action=exclude&index=9&on=1").0.contains("400"));
    control(d.url(), "action=include_all");
    assert_eq!(d.excluded(), vec![false; 4]);
}

#[test]
fn done_ends_wait_and_an_interrupt_stops_waiting() {
    let d = kernel_free();
    assert!(d.wait(Some(&|| true)).is_none());
    let url = d.url().to_string();
    let presser = std::thread::spawn(move || {
        std::thread::sleep(Duration::from_millis(200));
        control(&url, "action=done");
    });
    let fit = d.wait(None);
    presser.join().unwrap();
    assert!(fit.is_some());
    assert_eq!(json(d.url())["done"], true);
}

#[test]
fn dashboard_fits_and_robust_mode_rejects_the_outlier() {
    let Some((k, forces)) = kernel() else { return };
    let a = detections(&k, &forces);
    let d = FitDashboard::start(0, "main belt", a, k, FitOptions { forces, ..Default::default() }, None).unwrap();
    assert!(d.wait_idle(Duration::from_secs(60)));
    // The outlier pushes chi-square per degree of freedom over 10.
    let fit = d.fit().unwrap();
    assert!(!fit.converged(), "{:?}", fit.flag);
    let s = json(d.url());
    assert!(!s["live"].as_array().unwrap().is_empty(), "iterations are reported");

    control(d.url(), "action=set&key=robust&value=1");
    assert!(d.wait_idle(Duration::from_secs(60)));
    let fit = d.fit().unwrap();
    assert_eq!(fit.flag, FitFlag::Converged);
    assert!(!fit.used[17] && fit.used.iter().filter(|&&u| u).count() == 39);
    let s = json(d.url());
    let el = &s["fit"]["elements"];
    assert!((el["a"].as_f64().unwrap() - 2.08).abs() < 0.05, "{el}");
    assert!(el["sigma"][0].as_f64().unwrap() > 0.0);
    assert!(s["fit"]["rms"][2].as_f64().unwrap() < 1.0);

    // Leaving the outlier out by hand instead also converges, from the last orbit.
    control(d.url(), "action=set&key=robust&value=0");
    control(d.url(), "action=seed&on=1");
    control(d.url(), "action=exclude&index=17&on=1");
    assert!(d.wait_idle(Duration::from_secs(60)));
    let fit = d.fit().unwrap();
    assert_eq!(fit.flag, FitFlag::Converged);
    assert!(!fit.used[17] && fit.residuals[17][1].abs() > 20.0 * orbfit::DEFAULT_SIGMA);
}

#[test]
fn stop_ends_a_fit_early() {
    let Some((k, forces)) = kernel() else { return };
    let a = detections(&k, &forces);
    // A tight tolerance and many iterations, so the fit is still running when Stop arrives.
    let opts = FitOptions { forces, robust: true, epsilon: 1e-13, ..Default::default() };
    let d = FitDashboard::start(0, "slow", a, k, opts, None).unwrap();
    while json(d.url())["live"].as_array().unwrap().is_empty() {
        std::thread::sleep(Duration::from_millis(5));
    }
    control(d.url(), "action=stop");
    assert!(d.wait_idle(Duration::from_secs(20)));
    assert_eq!(json(d.url())["state"], "stopped");
}

/// Load and Random swap in another object's detections (here from a stand-in for the MPC),
/// fit them from scratch, and keep the current object when the fetch fails.
#[test]
fn load_and_random_swap_the_object() {
    let d = kernel_free();
    assert!(d.wait_idle(Duration::from_secs(10)));
    d.set_loader(|designation, _kernel| match designation {
        Some("nothing") => Err("the MPC has no observations of that object".into()),
        Some(name) => Astrometry::new(vec![2460100.5; 5], vec![2.0; 5], vec![-0.2; 5], vec![[1.0, 0.0, 0.0]; 5]).map(|a| (name.to_string(), a)),
        None => Astrometry::new(vec![2460200.5; 4], vec![3.0; 4], vec![0.3; 4], vec![[1.0, 0.0, 0.0]; 4]).map(|a| ("(12345)".to_string(), a)),
    });
    // From the page: a designation, then Random.
    control(d.url(), "action=exclude&index=0&on=1");
    assert!(control(d.url(), "action=load&designation=2020+AA").0.contains("200"));
    assert!(d.wait_idle(Duration::from_secs(10)));
    let s = json(d.url());
    assert_eq!(s["name"], "2020 AA");
    assert_eq!(s["detections"].as_array().unwrap().len(), 5);
    assert_eq!(s["excluded"], serde_json::json!(vec![false; 5]));
    assert_eq!((s["runs"].as_u64(), s["log"].as_array().unwrap().len()), (Some(1), 1));
    assert_eq!(s["loading"], false);
    assert!(s["message"].as_str().unwrap().contains("loaded 2020 AA"));
    assert_eq!(d.fit().unwrap().residuals.len(), 5);

    assert!(control(d.url(), "action=load").0.contains("200"));
    assert!(d.wait_idle(Duration::from_secs(10)));
    assert_eq!((d.name().as_str(), d.astrometry().len()), ("(12345)", 4));

    // A failed fetch keeps the object and says why.
    d.load(Some("nothing"));
    assert!(d.wait_idle(Duration::from_secs(10)));
    let s = json(d.url());
    assert_eq!((s["name"].as_str(), s["detections"].as_array().unwrap().len()), (Some("(12345)"), 4));
    assert!(s["message"].as_str().unwrap().contains("no observations"), "{}", s["message"]);
    assert_eq!(s["state"], "done");
    assert!(!d.loading());
}
