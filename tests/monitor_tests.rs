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
    let addr = url.trim_start_matches("http://").trim_end_matches('/');
    let mut stream = TcpStream::connect(addr).unwrap();
    write!(stream, "GET {path} HTTP/1.1\r\nHost: {addr}\r\n\r\n").unwrap();
    let mut response = String::new();
    stream.read_to_string(&mut response).unwrap();
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
    assert_eq!(s["test_particles"].as_array().unwrap().len(), 50);
    assert_eq!(s["integrator"], "IAS15");
    assert!(s["energy_error"].as_f64().unwrap() < 1e-12);

    let (status, body) = get(&url, "/");
    assert!(status.contains("200"));
    assert!(body.contains("<title>Simulation Monitor</title>"));
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
