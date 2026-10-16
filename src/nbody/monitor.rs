//! A live progress page for a running [`Simulation`].
//!
//! [`Simulation::monitor`] starts a small web server on `127.0.0.1` and returns its URL. Open it
//! in a browser (or a notebook iframe) to watch an integration: how far it has got towards its
//! target epoch, the current epoch and timestep, steps per second, an estimate of the time left,
//! the relative energy error of the massive bodies, and a top-down view of the particles.
//!
//! The simulation publishes a snapshot at most every [`PUBLISH_INTERVAL`] while it steps, so a
//! monitored run costs one clock read per step plus the occasional snapshot. The server runs on
//! its own thread with the standard library only, and stops once the simulation (and every
//! clone of it) is dropped or [`Simulation::stop_monitor`] is called.
use std::io::{Read, Write};
use std::net::{TcpListener, TcpStream};
use std::sync::{Arc, Mutex, Weak};
use std::time::{Duration, Instant};

use serde::Serialize;

use crate::constants::GRAVITATIONAL_CONSTANT;
use crate::nbody::Simulation;
use crate::state;

/// The least time between two snapshots.
pub const PUBLISH_INTERVAL: Duration = Duration::from_millis(100);

/// Most test particles drawn in the top-down view (every massive body is always drawn).
const MAX_DRAWN: usize = 2000;

/// Most points kept in the energy-error history.
const MAX_HISTORY: usize = 2000;

/// A shared handle to a progress page. Clones (and clones of the simulation) share the page.
#[derive(Clone)]
pub struct ProgressMonitor {
    shared: Arc<Shared>,
    url: String,
}

struct Shared {
    status: Mutex<Status>,
}

/// What the page shows; serialized as `/progress.json`.
#[derive(Serialize, Clone)]
struct Status {
    /// "idle" before the first run, "running" during one, "done" after it.
    state: &'static str,
    /// TDB Julian date the simulation's clock counts from.
    jd_ref: f64,
    /// Days since `jd_ref` when the run started, now, and at its target.
    t_start: f64,
    t: f64,
    t_target: Option<f64>,
    /// Steps since the run started, and every step since the monitor was attached.
    steps: u64,
    total_steps: u64,
    timestep: f64,
    /// Wall-clock seconds since the run started.
    elapsed: f64,
    n_particles: usize,
    n_massive: usize,
    integrator: String,
    /// Energy of the massive bodies when the monitor was attached.
    energy0: f64,
    /// `|E - E0| / |E0|` of the massive bodies.
    energy_error: f64,
    /// `(t, energy_error)` pairs, oldest first.
    history: Vec<(f64, f64)>,
    /// Positions `[x, y]` (AU) of the massive bodies, with their names, then of up to
    /// [`MAX_DRAWN`] test particles.
    massive: Vec<(String, f64, f64)>,
    test_particles: Vec<(f64, f64)>,
    reference_plane: String,
    origin: String,
    /// The target set with [`ProgressMonitor::set_target`] overrides each run's own.
    #[serde(skip)]
    fixed_target: Option<(f64, f64)>,
    #[serde(skip)]
    run_started: Option<Instant>,
    #[serde(skip)]
    last_publish: Option<Instant>,
}

impl ProgressMonitor {
    /// Starts the page's server on `127.0.0.1:port` (`port` 0 picks a free port).
    fn start(port: u16) -> std::io::Result<ProgressMonitor> {
        let listener = TcpListener::bind(("127.0.0.1", port))?;
        listener.set_nonblocking(true)?;
        let url = format!("http://{}/", listener.local_addr()?);
        let shared = Arc::new(Shared {
            status: Mutex::new(Status {
                state: "idle",
                jd_ref: 0.0,
                t_start: 0.0,
                t: 0.0,
                t_target: None,
                steps: 0,
                total_steps: 0,
                timestep: 0.0,
                elapsed: 0.0,
                n_particles: 0,
                n_massive: 0,
                integrator: String::new(),
                energy0: 0.0,
                energy_error: 0.0,
                history: Vec::new(),
                massive: Vec::new(),
                test_particles: Vec::new(),
                reference_plane: String::new(),
                origin: String::new(),
                fixed_target: None,
                run_started: None,
                last_publish: None,
            }),
        });
        let weak = Arc::downgrade(&shared);
        std::thread::Builder::new()
            .name("spacerocks-monitor".into())
            .spawn(move || serve(listener, weak))?;
        Ok(ProgressMonitor { shared, url })
    }

    /// The page's address, like `http://127.0.0.1:8765/`.
    pub fn url(&self) -> &str {
        &self.url
    }

    fn lock(&self) -> std::sync::MutexGuard<'_, Status> {
        self.shared.status.lock().unwrap_or_else(|e| e.into_inner())
    }

    /// Measure progress from `t_start` to `t_target` (days since the simulation's reference
    /// epoch) instead of over each call's own span, e.g. across a loop of
    /// [`Simulation::integrate_or_interpolate`] calls.
    fn set_target(&self, t_start: f64, t_target: f64) {
        self.lock().fixed_target = Some((t_start, t_target));
    }

    /// A run (one call of `integrate`, `integrate_or_interpolate` or `steps`) starts.
    pub(crate) fn begin_run(&self, sim: &Simulation, t_target: Option<f64>) {
        let mut s = self.lock();
        let (t_start, t_target) = match s.fixed_target {
            Some((start, target)) => (start, Some(target)),
            None => (sim.t(), t_target),
        };
        if s.state != "running" || s.fixed_target.is_none() {
            s.run_started = Some(Instant::now());
            s.steps = 0;
        }
        s.t_start = t_start;
        s.t_target = t_target;
        s.state = "running";
        snapshot(&mut s, sim);
    }

    /// `n` steps were taken; publish a snapshot if the last one is old enough.
    pub(crate) fn tick(&self, sim: &Simulation, n: usize) {
        let mut s = self.lock();
        s.steps += n as u64;
        s.total_steps += n as u64;
        if s.last_publish.map_or(true, |last| last.elapsed() >= PUBLISH_INTERVAL) {
            snapshot(&mut s, sim);
        }
    }

    /// The run ended. With a fixed target it is only done once that target is reached.
    pub(crate) fn end_run(&self, sim: &Simulation) {
        let mut s = self.lock();
        let finished = match s.fixed_target {
            Some((start, target)) => (sim.t() - target) * (target - start).signum() >= -1e-9,
            None => true,
        };
        if finished {
            s.state = "done";
        }
        snapshot(&mut s, sim);
    }
}

/// Copy what the page shows out of the simulation.
fn snapshot(s: &mut Status, sim: &Simulation) {
    let energy = massive_energy(sim);
    if s.energy0 == 0.0 {
        s.energy0 = energy;
    }
    s.energy_error = if s.energy0 != 0.0 { ((energy - s.energy0) / s.energy0).abs() } else { 0.0 };
    if s.history.last().map_or(true, |&(t, _)| t != sim.t()) {
        if s.history.len() == MAX_HISTORY {
            // Keep every other point, so the history still spans the whole run.
            let thinned = s.history.iter().step_by(2).copied().collect();
            s.history = thinned;
        }
        s.history.push((sim.t(), s.energy_error));
    }

    s.jd_ref = sim.jd_ref();
    s.t = sim.t();
    s.timestep = sim.integrator.timestep();
    s.elapsed = s.run_started.map_or(0.0, |start| start.elapsed().as_secs_f64());
    s.n_particles = sim.len();
    s.integrator = integrator_name(sim);
    s.reference_plane = sim.reference_plane.to_string();
    s.origin = sim.origin.to_string();

    let masses = sim.masses();
    let states = sim.states();
    let names = sim.names();
    s.n_massive = masses.iter().filter(|&&m| m != 0.0).count();
    s.massive = (0..states.len())
        .filter(|&i| masses[i] != 0.0)
        .map(|i| (names[i].clone(), states[i][0], states[i][1]))
        .collect();
    let n_test = states.len() - s.n_massive;
    let stride = n_test.div_ceil(MAX_DRAWN).max(1);
    s.test_particles = (0..states.len())
        .filter(|&i| masses[i] == 0.0)
        .step_by(stride)
        .map(|i| (states[i][0], states[i][1]))
        .collect();
    s.last_publish = Some(Instant::now());
}

/// Energy of the massive bodies (test particles have none), so it costs nothing per test
/// particle.
fn massive_energy(sim: &Simulation) -> f64 {
    let masses = sim.masses();
    let states = sim.states();
    let massive: Vec<usize> = (0..masses.len()).filter(|&i| masses[i] != 0.0).collect();
    let mut energy = 0.0;
    for (k, &i) in massive.iter().enumerate() {
        let (x_i, v_i) = state::pv(&states[i]);
        energy += 0.5 * masses[i] * v_i.norm_squared();
        for &j in &massive[k + 1..] {
            let r = (x_i - state::position(&states[j])).norm();
            energy -= GRAVITATIONAL_CONSTANT * masses[i] * masses[j] / r;
        }
    }
    energy
}

/// The integrator's type name, like `IAS15`.
fn integrator_name(sim: &Simulation) -> String {
    let full = sim.integrator.name();
    full.rsplit("::").next().unwrap_or(full).to_string()
}

/// Answer requests until the monitor is dropped.
fn serve(listener: TcpListener, shared: Weak<Shared>) {
    loop {
        if shared.strong_count() == 0 {
            return;
        }
        match listener.accept() {
            Ok((stream, _)) => {
                let Some(shared) = shared.upgrade() else { return };
                let _ = respond(stream, &shared);
            }
            Err(e) if e.kind() == std::io::ErrorKind::WouldBlock => {
                std::thread::sleep(Duration::from_millis(25));
            }
            Err(_) => std::thread::sleep(Duration::from_millis(25)),
        }
    }
}

fn respond(mut stream: TcpStream, shared: &Shared) -> std::io::Result<()> {
    stream.set_nonblocking(false)?;
    stream.set_read_timeout(Some(Duration::from_secs(2)))?;
    let mut buf = [0u8; 2048];
    let n = stream.read(&mut buf)?;
    let request = String::from_utf8_lossy(&buf[..n]);
    let path = request.split_whitespace().nth(1).unwrap_or("/");
    let path = path.split('?').next().unwrap_or("/");

    let (status, content_type, body) = match path {
        "/" | "/index.html" => ("200 OK", "text/html; charset=utf-8", PAGE.to_string()),
        "/progress.json" => {
            let status = shared.status.lock().unwrap_or_else(|e| e.into_inner()).clone();
            ("200 OK", "application/json", serde_json::to_string(&status).unwrap_or_default())
        }
        _ => ("404 Not Found", "text/plain", "not found".to_string()),
    };
    write!(
        stream,
        "HTTP/1.1 {status}\r\nContent-Type: {content_type}\r\nContent-Length: {}\r\nCache-Control: no-store\r\nConnection: close\r\n\r\n",
        body.len()
    )?;
    stream.write_all(body.as_bytes())?;
    stream.flush()
}

impl Simulation {
    /// Serve a live progress page for this simulation on `127.0.0.1:port` and return its URL
    /// (`port` 0 picks a free port). Replaces any page already attached.
    ///
    /// ```no_run
    /// # use spacerocks::{Simulation, Time};
    /// # let mut sim = Simulation::default();
    /// let url = sim.monitor(8765).unwrap();
    /// println!("watch at {url}");
    /// sim.integrate(&Time::new(2470000.5, "tdb", "jd").unwrap());
    /// ```
    pub fn monitor(&mut self, port: u16) -> std::io::Result<String> {
        let monitor = ProgressMonitor::start(port)?;
        let url = monitor.url().to_string();
        snapshot(&mut monitor.lock(), self);
        self.progress_monitor = Some(monitor);
        Ok(url)
    }

    /// Measure the progress page's progress up to `epoch` from the current epoch, across
    /// however many calls it takes to get there (e.g. a loop of
    /// [`Simulation::integrate_or_interpolate`] over output epochs), instead of over each
    /// call. Does nothing without a page attached.
    pub fn set_monitor_target(&mut self, epoch: &crate::Time) {
        let target = epoch.tdb().jd() - self.jd_ref();
        if let Some(monitor) = &self.progress_monitor {
            monitor.set_target(self.t(), target);
        }
    }

    /// The progress page's URL, if one is attached.
    pub fn monitor_url(&self) -> Option<&str> {
        self.progress_monitor.as_ref().map(|m| m.url())
    }

    /// Detach the progress page; its server stops once no clone of this simulation holds it.
    pub fn stop_monitor(&mut self) {
        self.progress_monitor = None;
    }
}

const PAGE: &str = include_str!("monitor.html");
