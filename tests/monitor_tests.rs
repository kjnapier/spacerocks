//! The live progress page of an nbody Simulation (no kernels needed).

use std::io::{Read, Write};
use std::net::TcpStream;

use spacerocks::constants::GRAVITATIONAL_CONSTANT;
use spacerocks::nbody::{Integrator, Simulation, WisdomHolman, IAS15};
use spacerocks::{SpaceRock, Time};

const T0: f64 = 2460000.5;

fn rock(name: &str, mass: f64, x: [f64; 3], v: [f64; 3]) -> SpaceRock {
    let mut r = SpaceRock::from_xyz(name, x[0], x[1], x[2], v[0], v[1], v[2], Time::new(T0, "tdb", "jd").unwrap(), "J2000", "SSB").unwrap();
    r.set_mass(mass);
    r
}

fn outer_system(integrator: Box<dyn Integrator + Send + Sync>) -> Simulation {
    let mut sim = Simulation::new(&Time::new(T0, "tdb", "jd").unwrap(), "J2000", "SSB").unwrap();
    sim.integrator = integrator;
    let vc = |a: f64| (GRAVITATIONAL_CONSTANT / a).sqrt();
    sim.add(rock("sun", 1.0, [0.0; 3], [0.0; 3])).unwrap();
    sim.add(rock("jupiter", 9.547919e-4, [5.2, 0.0, 0.0], [0.0, vc(5.2) * 1.02, 2e-4])).unwrap();
    sim.add(rock("saturn", 2.858860e-4, [0.0, 9.55, 0.3], [-vc(9.55) * 0.97, 0.0, 0.0])).unwrap();
    for k in 0..50 {
        let a = 2.0 + 0.05 * k as f64;
        sim.add(rock(&format!("tp{k}"), 0.0, [a, 0.0, 0.0], [0.0, vc(a), 1e-4])).unwrap();
    }
    sim
}

/// GET `path` from the page at `url`; returns (status line, body).
fn get(url: &str, path: &str) -> (String, String) {
    request(url, "GET", path, "")
}

/// POST a control action as the page does.
fn control(url: &str, query: &str) -> String {
    request(url, "POST", &format!("/control?{query}"), "X-Spacerocks: 1\r\n").0
}

fn json(url: &str) -> serde_json::Value {
    serde_json::from_str(&get(url, "/progress.json").1).unwrap()
}

/// Wait (up to a few seconds) until the page's JSON satisfies `ok`.
fn wait_for(url: &str, ok: impl Fn(&serde_json::Value) -> bool) -> serde_json::Value {
    for _ in 0..300 {
        let s = json(url);
        if ok(&s) {
            return s;
        }
        std::thread::sleep(std::time::Duration::from_millis(10));
    }
    panic!("page never got there: {}", json(url));
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

#[test]
fn page_reports_a_finished_integration() {
    let mut sim = outer_system(Box::new(IAS15::new(1.0)));
    sim.move_to_center_of_mass().unwrap();
    let url = sim.monitor(0).unwrap();
    assert!(url.starts_with("http://127.0.0.1:"));
    assert_eq!(sim.monitor_url(), Some(url.as_str()));

    let (status, body) = get(&url, "/progress.json");
    assert!(status.contains("200"));
    let s: serde_json::Value = serde_json::from_str(&body).unwrap();
    assert_eq!(s["state"], "idle");

    sim.integrate(&Time::new(T0 + 3650.0, "tdb", "jd").unwrap());
    let (_, body) = get(&url, "/progress.json");
    let s: serde_json::Value = serde_json::from_str(&body).unwrap();
    assert_eq!(s["state"], "done");
    assert_eq!(s["t"], 3650.0);
    assert_eq!(s["t_target"], 3650.0);
    assert_eq!(s["jd_ref"], T0);
    assert!(s["steps"].as_u64().unwrap() > 100);
    assert_eq!(s["n_particles"], 53);
    assert_eq!(s["n_massive"], 3);
    assert_eq!(s["massive"].as_array().unwrap().len(), 3);
    assert_eq!(s["massive"][1][0], "jupiter");
    assert_eq!(s["test_particles"][0][0], "tp0");
    assert_eq!(s["test_particles"].as_array().unwrap().len(), 50);
    assert_eq!(s["integrator"], "IAS15");
    assert!(s["energy_error"].as_f64().unwrap() < 1e-12);

    let (status, body) = get(&url, "/");
    assert!(status.contains("200"));
    assert!(body.contains("<title>Simulation Monitor</title>"));
    let (status, body) = get(&url, "/font.ttf");
    assert!(status.contains("200") && body.len() > 1000);
    let (status, _) = get(&url, "/nope");
    assert!(status.contains("404"));
}

#[test]
fn monitor_target_spans_several_calls() {
    let mut sim = outer_system(Box::new(IAS15::new(1.0)));
    let url = sim.monitor(0).unwrap();
    sim.set_monitor_target(&Time::new(T0 + 1000.0, "tdb", "jd").unwrap());
    for k in 1..=4 {
        sim.integrate_or_interpolate(&Time::new(T0 + 250.0 * k as f64, "tdb", "jd").unwrap());
        let (_, body) = get(&url, "/progress.json");
        let s: serde_json::Value = serde_json::from_str(&body).unwrap();
        assert_eq!(s["t_start"], 0.0);
        assert_eq!(s["t_target"], 1000.0);
        assert_eq!(s["state"], if k < 4 { "running" } else { "done" });
    }
}

#[test]
fn monitoring_does_not_change_the_integration() {
    // Wisdom–Holman hands `steps` to the integrator in chunks while monitored.
    let target = Time::new(T0 + 20.0 * 2500.0, "tdb", "jd").unwrap();
    let mut plain = outer_system(Box::new(WisdomHolman::new(20.0)));
    let mut watched = plain.clone();
    watched.monitor(0).unwrap();
    plain.integrate(&target);
    watched.integrate(&target);
    assert_eq!(plain.states(), watched.states());
    assert_eq!(plain.t(), watched.t());
}

#[test]
fn server_stops_with_the_simulation() {
    let mut sim = outer_system(Box::new(IAS15::new(1.0)));
    let url = sim.monitor(0).unwrap();
    let addr = url.trim_start_matches("http://").trim_end_matches('/').to_string();
    sim.stop_monitor();
    assert_eq!(sim.monitor_url(), None);
    let mut stopped = false;
    for _ in 0..100 {
        if TcpStream::connect(&addr).is_err() {
            stopped = true;
            break;
        }
        std::thread::sleep(std::time::Duration::from_millis(20));
    }
    assert!(stopped);
}

#[test]
fn page_pauses_steps_and_switches_integrator() {
    let mut sim = outer_system(Box::new(IAS15::new(1.0)));
    let url = sim.monitor(0).unwrap();
    sim.progress_monitor().unwrap().set_paused(true);
    let target = Time::new(T0 + 2000.0, "tdb", "jd").unwrap();
    let run = std::thread::spawn(move || {
        sim.integrate(&target);
        sim
    });

    let s = wait_for(&url, |s| s["state"] == "paused");
    assert_eq!(s["paused"], true);
    assert_eq!(s["steps"], 0);

    for k in 1..=3 {
        assert!(control(&url, "action=step").contains("200"));
        wait_for(&url, |s| s["steps"] == k && s["state"] == "paused");
    }

    assert!(control(&url, "action=integrator&name=WisdomHolman&timestep=10").contains("200"));
    let s = wait_for(&url, |s| s["integrator"] == "WisdomHolman");
    assert_eq!(s["timestep"], 10.0);
    assert_eq!(s["steps"], 3);
    assert!(control(&url, "action=timestep&timestep=5").contains("200"));
    wait_for(&url, |s| s["timestep"] == 5.0);
    let t = json(&url)["t"].as_f64().unwrap();
    control(&url, "action=step");
    let s = wait_for(&url, |s| s["steps"] == 4);
    assert!((s["t"].as_f64().unwrap() - t - 5.0).abs() < 1e-9);

    control(&url, "action=resume");
    let sim = run.join().unwrap();
    assert_eq!(sim.t(), 2000.0);
    let s = json(&url);
    assert_eq!(s["state"], "done");
    assert_eq!(s["paused"], false);
    assert_eq!(s["integrator"], "WisdomHolman");
}

#[test]
fn control_needs_the_page() {
    let mut sim = outer_system(Box::new(IAS15::new(1.0)));
    let url = sim.monitor(0).unwrap();
    // No marker header, as a form on another site would send it.
    assert!(request(&url, "POST", "/control?action=pause", "").0.contains("403"));
    assert_eq!(json(&url)["paused"], false);
    assert!(control(&url, "action=launch").contains("400"));
    assert!(control(&url, "action=integrator&name=RK4&timestep=1").contains("400"));
    assert!(control(&url, "action=timestep&timestep=0").contains("400"));
    assert!(control(&url, "action=pause").contains("200"));
    assert_eq!(json(&url)["paused"], true);
}

#[test]
fn page_tracks_a_particles_elements() {
    let mut sim = outer_system(Box::new(IAS15::new(1.0)));
    let url = sim.monitor(0).unwrap();
    // Before any run: a particle drawn on the page is tracked from the last snapshot at once.
    assert!(control(&url, "action=track&name=tp10").contains("200"));
    let s = json(&url);
    assert_eq!(s["tracked"], "tp10");
    let el = &s["elements"];
    assert_eq!(el["name"], "tp10");
    assert_eq!(el["central"], "sun");
    assert!((el["a"].as_f64().unwrap() - 2.5).abs() < 0.01, "{el}");
    assert!(el["e"].as_f64().unwrap() < 0.01);
    assert!((el["inc"].as_f64().unwrap() - 0.53).abs() < 0.05, "{el}");
    assert!(el["period"].as_f64().unwrap() > 1400.0);
    assert_eq!(s["element_history"].as_array().unwrap().len(), 1);

    sim.integrate(&Time::new(T0 + 1000.0, "tdb", "jd").unwrap());
    let s = json(&url);
    let history = s["element_history"].as_array().unwrap();
    assert!(history.len() >= 2);
    assert_eq!(history.last().unwrap()[0], 1000.0);

    // Names are percent-decoded; unknown names and the central body have no elements.
    control(&url, "action=track&name=no%20such+rock");
    let s = json(&url);
    assert_eq!(s["tracked"], "no such rock");
    assert!(s["elements"].is_null());
    control(&url, "action=track&name=sun");
    assert!(json(&url)["elements"].is_null());
    sim.progress_monitor().unwrap().track(Some("jupiter"));
    assert!((json(&url)["elements"]["a"].as_f64().unwrap() - 5.4).abs() < 0.5);
    control(&url, "action=track&name=");
    let s = json(&url);
    assert!(s["tracked"].is_null() && s["element_history"].as_array().unwrap().is_empty());
}

#[test]
fn interrupt_stops_running_and_paused_runs() {
    use std::sync::atomic::{AtomicBool, Ordering};
    use std::sync::Arc;
    use std::time::{Duration, Instant};

    for paused in [false, true] {
        let mut sim = outer_system(Box::new(WisdomHolman::new(1.0)));
        let url = sim.monitor(0).unwrap();
        sim.progress_monitor().unwrap().set_paused(paused);
        let stop = Arc::new(AtomicBool::new(false));
        let flag = stop.clone();
        sim.set_interrupt(Some(Arc::new(move || flag.load(Ordering::Relaxed))));
        let setter = stop.clone();
        std::thread::spawn(move || {
            std::thread::sleep(Duration::from_millis(200));
            setter.store(true, Ordering::Relaxed);
        });

        let start = Instant::now();
        sim.integrate(&Time::new(T0 + 1e9, "tdb", "jd").unwrap());
        assert!(start.elapsed() < Duration::from_secs(3));
        assert!(sim.interrupted());
        assert!(sim.t() < 1e9);
        if paused {
            assert_eq!(sim.t(), 0.0);
        }
        assert_eq!(json(&url)["state"], "stopped");

        // Removing the check lets the simulation run again.
        sim.set_interrupt(None);
        sim.progress_monitor().unwrap().set_paused(false);
        let t = sim.t();
        sim.integrate(&Time::new(T0 + t + 10.0, "tdb", "jd").unwrap());
        assert!((sim.t() - t - 10.0).abs() < 1e-9);
        assert!(!sim.interrupted());
    }
}

#[test]
fn page_steps_back_and_moves_the_end() {
    let mut sim = outer_system(Box::new(WisdomHolman::new(10.0)));
    let url = sim.monitor(0).unwrap();
    sim.progress_monitor().unwrap().set_paused(true);
    let target = Time::new(T0 + 2000.0, "tdb", "jd").unwrap();
    let run = std::thread::spawn(move || {
        sim.integrate(&target);
        sim
    });

    let s = wait_for(&url, |s| s["state"] == "paused");
    assert_eq!(s["end_editable"], true);
    control(&url, "action=step");
    wait_for(&url, |s| s["steps"] == 1 && s["state"] == "paused");
    control(&url, "action=step");
    let s = wait_for(&url, |s| s["steps"] == 2 && s["state"] == "paused");
    assert!((s["t"].as_f64().unwrap() - 20.0).abs() < 1e-9);
    assert!(control(&url, "action=back").contains("200"));
    let s = wait_for(&url, |s| s["steps"] == 3);
    assert!((s["t"].as_f64().unwrap() - 10.0).abs() < 1e-9);
    assert_eq!(s["timestep"], 10.0, "Back leaves the run going forwards");

    assert!(control(&url, &format!("action=end&jd={}", T0 + 105.0)).contains("200"));
    let s = wait_for(&url, |s| s["t_target"] == 105.0);
    assert_eq!(s["state"], "paused");
    control(&url, "action=resume");
    let sim = run.join().unwrap();
    assert!((sim.t() - 105.0).abs() < 1e-12);
    assert_eq!(json(&url)["state"], "done");
}

#[test]
fn open_ended_runs_go_until_stop() {
    let mut sim = outer_system(Box::new(WisdomHolman::new(1.0)));
    let url = sim.monitor(0).unwrap();
    let run = std::thread::spawn(move || {
        sim.run();
        sim
    });
    let s = wait_for(&url, |s| s["state"] == "running" && s["steps"].as_u64().unwrap() > 100);
    assert!(s["t_target"].is_null());
    assert!(control(&url, "action=stop").contains("200"));
    let sim = run.join().unwrap();
    assert!(sim.t() > 100.0);
    let s = json(&url);
    assert_eq!(s["state"], "stopped");
    assert_eq!(s["t"].as_f64().unwrap(), sim.t());

    // An integrate whose end the page takes away runs on until Stop too.
    let mut sim = sim;
    let url = sim.monitor(0).unwrap();
    sim.progress_monitor().unwrap().set_paused(true);
    let start = sim.t();
    let target = Time::new(T0 + start + 10.0, "tdb", "jd").unwrap();
    let run = std::thread::spawn(move || {
        sim.integrate(&target);
        sim
    });
    wait_for(&url, |s| s["state"] == "paused");
    control(&url, "action=end");
    wait_for(&url, |s| s["t_target"].is_null());
    control(&url, "action=resume");
    wait_for(&url, |s| s["t"].as_f64().unwrap() > start + 50.0);
    control(&url, "action=stop");
    let sim = run.join().unwrap();
    assert!(sim.t() > start + 50.0);
}

#[test]
fn page_toggles_gr_and_perturbers() {
    let mut sim = outer_system(Box::new(IAS15::new(1.0)));
    assert!(!sim.solar_gr());
    assert_eq!(sim.perturbers(), vec![("jupiter".to_string(), true), ("saturn".to_string(), true)]);
    let url = sim.monitor(0).unwrap();
    sim.progress_monitor().unwrap().set_paused(true);
    let target = Time::new(T0 + 100.0, "tdb", "jd").unwrap();
    let run = std::thread::spawn(move || {
        sim.integrate(&target);
        sim
    });
    let s = wait_for(&url, |s| s["state"] == "paused");
    assert_eq!(s["solar_gr"], false);
    assert_eq!(s["n_massive"], 3);

    assert!(control(&url, "action=gr&on=1").contains("200"));
    wait_for(&url, |s| s["solar_gr"] == true);
    assert!(control(&url, "action=perturber&name=saturn&on=0").contains("200"));
    let s = wait_for(&url, |s| s["n_massive"] == 2);
    assert_eq!(s["perturbers"], serde_json::json!([["jupiter", true], ["saturn", false]]));
    // Saturn is still drawn with the massive bodies, just without gravity.
    assert!(s["massive"].as_array().unwrap().iter().any(|p| p[0] == "saturn"));
    // The central body can't be turned off; a missing body is ignored.
    control(&url, "action=perturber&name=sun&on=0");
    control(&url, "action=perturber&name=pluto&on=0");
    assert!(control(&url, "action=perturber&name=saturn&on=maybe").contains("400"));
    control(&url, "action=step");
    let s = wait_for(&url, |s| s["steps"] == 1);
    assert_eq!(s["n_massive"], 2);

    control(&url, "action=resume");
    let mut sim = run.join().unwrap();
    assert!(sim.solar_gr());
    assert_eq!(sim.masses()[2], 0.0);
    assert_eq!(sim.get_particle("saturn").unwrap().mass(), 2.858860e-4);
    sim.set_perturber("saturn", true).unwrap();
    assert_eq!(sim.masses()[2], 2.858860e-4);
    sim.set_solar_gr(false);
    assert!(!sim.solar_gr());
    assert_eq!(sim.forces.len(), 1);
}

#[test]
fn muted_perturber_matches_a_simulation_without_it() {
    let mut with = outer_system(Box::new(IAS15::new(1.0)));
    with.set_perturber("saturn", false).unwrap();
    let mut without = outer_system(Box::new(IAS15::new(1.0)));
    without.remove("saturn").unwrap();
    let target = Time::new(T0 + 500.0, "tdb", "jd").unwrap();
    with.integrate(&target);
    without.integrate(&target);
    let a = with.get_particle("jupiter").unwrap();
    let b = without.get_particle("jupiter").unwrap();
    assert!((a.position - b.position).norm() < 1e-12, "{}", (a.position - b.position).norm());
}
