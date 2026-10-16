//! A live progress page for a running [`Simulation`].
//!
//! [`Simulation::monitor`] starts a small web server on `127.0.0.1` and returns its URL. Open it
//! in a browser (or a notebook iframe) to watch an integration: how far it has got towards its
//! target epoch, the current epoch and timestep, steps per second, an estimate of the time left,
//! the relative energy error of the massive bodies, a top-down view of the particles (hover one
//! for its name and position) and the osculating elements of one tracked particle over time
//! (click one, or [`ProgressMonitor::track`]).
//!
//! The simulation publishes a snapshot at most every [`PUBLISH_INTERVAL`] while it steps, so a
//! monitored run costs one clock read per step plus the occasional snapshot. The server runs on
//! its own thread with the standard library only, and stops once the simulation (and every
//! clone of it) is dropped or [`Simulation::stop_monitor`] is called.
//!
//! The page also steers the run. Pause holds the simulation before its next step (the calling
//! thread waits, so `integrate` returns later, not early); Step lets one step through and pauses
//! again, and Back takes one step backwards; the integrator and timestep can be swapped, solar
//! general relativity turned on or off, and each perturber's gravity turned off or on, taking
//! effect at the next step. The run's end can be moved, or removed so it runs until Stop (which
//! ends it where it is, so `integrate` returns early). The page asks with `POST /control?action=...`, which needs an `X-Spacerocks` header and
//! a local `Host`, so other web pages open in the browser can't send it.
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Condvar, Mutex};
use std::time::{Duration, Instant};

use serde::Serialize;

use crate::constants::GRAVITATIONAL_CONSTANT;
use crate::nbody::Simulation;
use crate::state::{self, State};
use crate::web::{self, Request, Response};

/// The least time between two snapshots.
pub const PUBLISH_INTERVAL: Duration = Duration::from_millis(100);

/// Most test particles drawn in the top-down view (every massive body is always drawn).
const MAX_DRAWN: usize = 2000;

/// Most points kept in the energy-error history.
const MAX_HISTORY: usize = 2000;

/// Osculating elements of the tracked particle about the most massive body. Angles in degrees.
#[derive(Serialize, Clone)]
struct TrackedElements {
    name: String,
    central: String,
    a: f64,
    e: f64,
    inc: f64,
    node: f64,
    arg: f64,
    mean_anomaly: f64,
    true_anomaly: f64,
    q: f64,
    /// Aphelion distance and period (days), for bound orbits.
    big_q: Option<f64>,
    period: Option<f64>,
    /// Distance from the central body.
    r: f64,
}

/// A shared handle to a progress page. Clones (and clones of the simulation) share the page.
#[derive(Clone)]
pub struct ProgressMonitor {
    shared: Arc<Shared>,
    url: String,
}

struct Shared {
    status: Mutex<Status>,
    control: Mutex<Control>,
    /// Signalled when `control` changes.
    wake: Condvar,
    /// Whether `control` is paused or holds commands, so a step can skip its lock otherwise.
    attention: AtomicBool,
}

/// What the page has asked for.
#[derive(Default)]
struct Control {
    paused: bool,
    /// Steps let through while paused, and steps backwards asked for.
    step_credits: u64,
    back_credits: u64,
    commands: Vec<Command>,
}

/// What a paused simulation may do next ([`ProgressMonitor::wait`]).
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) enum Release {
    /// Run on.
    Run,
    /// Take one step, then wait again.
    One,
    /// Take one step backwards, then wait again.
    Back,
}

/// A change to the simulation asked for on the page.
pub(crate) enum Command {
    /// Switch to the integrator of this name (`IAS15`, `WisdomHolman`, `Leapfrog`, `Trace`) with
    /// this timestep (days).
    Integrator(String, f64),
    /// Set the integrator's timestep (days).
    Timestep(f64),
    /// End the run at this time (days since the reference epoch), or never (`None`): it then
    /// runs until stopped.
    End(Option<f64>),
    /// End the run now.
    Stop,
    /// Turn solar general relativity ([`SolarGR`](crate::nbody::forces::SolarGR)) on or off.
    SolarGR(bool),
    /// Turn the gravity of the massive particle of this name on or off.
    Perturber(String, bool),
}

/// The integrators the page offers.
const INTEGRATORS: [&str; 4] = ["IAS15", "WisdomHolman", "Leapfrog", "Trace"];

/// What the page shows; serialized as `/progress.json`.
#[derive(Serialize, Clone)]
struct Status {
    /// "idle" before the first run, "running" during one, "paused" while held in one, "done"
    /// after it, "stopped" after one cut short by an interrupt.
    state: &'static str,
    /// Whether the page's Pause is on (the next step will wait).
    paused: bool,
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
    /// Names and positions `(name, x, y, z)` (AU) of the massive bodies, then of up to
    /// [`MAX_DRAWN`] test particles.
    massive: Vec<(String, f64, f64, f64)>,
    test_particles: Vec<(String, f64, f64, f64)>,
    /// Whether the run's end can be changed on the page (not while
    /// [`Simulation::set_monitor_target`] spans several calls).
    end_editable: bool,
    /// Whether solar general relativity is on.
    solar_gr: bool,
    /// Every particle with a mass but the central body, and whether its gravity is on.
    perturbers: Vec<(String, bool)>,
    /// The particle whose elements are followed, if any, and its elements now.
    tracked: Option<String>,
    elements: Option<TrackedElements>,
    /// `[t, a, e, inc, q]` of the tracked particle, oldest first.
    element_history: Vec<[f64; 5]>,
    /// States and masses of the drawn particles (massive first), so a particle clicked on the
    /// page can be tracked at once, and the central body's name, state and mass.
    #[serde(skip)]
    drawn: Vec<(State, f64)>,
    #[serde(skip)]
    central: Option<(String, State, f64)>,
    reference_plane: String,
    origin: String,
    /// The target set with [`ProgressMonitor::set_target`] overrides each run's own.
    #[serde(skip)]
    fixed_target: Option<(f64, f64)>,
    #[serde(skip)]
    run_started: Option<Instant>,
    #[serde(skip)]
    last_publish: Option<Instant>,
    /// Time spent paused in this run, left out of `elapsed`.
    #[serde(skip)]
    paused_total: Duration,
    #[serde(skip)]
    paused_since: Option<Instant>,
    /// Whether the page's Stop ended the run.
    #[serde(skip)]
    stop_asked: bool,
}

impl ProgressMonitor {
    /// Starts the page's server on `127.0.0.1:port` (`port` 0 picks a free port).
    fn start(port: u16) -> std::io::Result<ProgressMonitor> {
        let shared = Arc::new(Shared {
            status: Mutex::new(Status {
                state: "idle",
                paused: false,
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
                end_editable: false,
                solar_gr: false,
                perturbers: Vec::new(),
                tracked: None,
                elements: None,
                element_history: Vec::new(),
                drawn: Vec::new(),
                central: None,
                reference_plane: String::new(),
                origin: String::new(),
                fixed_target: None,
                run_started: None,
                last_publish: None,
                paused_total: Duration::ZERO,
                paused_since: None,
                stop_asked: false,
            }),
            control: Mutex::new(Control::default()),
            wake: Condvar::new(),
            attention: AtomicBool::new(false),
        });
        let url = web::start(port, "spacerocks-monitor", Arc::downgrade(&shared), respond)?;
        Ok(ProgressMonitor { shared, url })
    }

    /// The page's address, like `http://127.0.0.1:8765/`.
    pub fn url(&self) -> &str {
        &self.url
    }

    fn lock(&self) -> std::sync::MutexGuard<'_, Status> {
        self.shared.status.lock().unwrap_or_else(|e| e.into_inner())
    }

    fn control(&self) -> std::sync::MutexGuard<'_, Control> {
        self.shared.control.lock().unwrap_or_else(|e| e.into_inner())
    }

    /// Hold the simulation while the page is paused. Returns the commands the page sent, and
    /// what the simulation may do next, or `None` to come back after applying the commands, or
    /// after a quarter of a second of waiting (so the caller can check for an interrupt such as
    /// Ctrl-C).
    pub(crate) fn wait(&self, sim: &Simulation) -> (Vec<Command>, Option<Release>) {
        if !self.shared.attention.load(Ordering::Acquire) {
            return (Vec::new(), Some(Release::Run));
        }
        let mut control = self.control();
        let mut announced = false;
        loop {
            let commands = std::mem::take(&mut control.commands);
            self.shared.attention.store(control.paused || control.back_credits > 0, Ordering::Release);
            let released = if control.back_credits > 0 {
                control.back_credits -= 1;
                Some(Release::Back)
            } else if !control.paused {
                Some(Release::Run)
            } else if control.step_credits > 0 {
                control.step_credits -= 1;
                Some(Release::One)
            } else {
                None
            };
            if !commands.is_empty() || released.is_some() {
                drop(control);
                if announced || released.is_some() {
                    let mut s = self.lock();
                    if s.state == "paused" {
                        s.state = "running";
                    }
                    if let Some(since) = s.paused_since.take() {
                        s.paused_total += since.elapsed();
                    }
                }
                return (commands, released);
            }
            if !announced {
                drop(control);
                {
                    let mut s = self.lock();
                    if s.state == "running" {
                        s.state = "paused";
                    }
                    if s.paused_since.is_none() {
                        s.paused_since = Some(Instant::now());
                    }
                    snapshot(&mut s, sim);
                }
                announced = true;
                control = self.control();
                continue;
            }
            control = self.shared.wake.wait_timeout(control, Duration::from_millis(250)).unwrap_or_else(|e| e.into_inner()).0;
            if control.paused && control.step_credits == 0 && control.back_credits == 0 && control.commands.is_empty() {
                // Still paused: let the caller look for an interrupt, then wait again.
                return (Vec::new(), None);
            }
        }
    }

    /// Publish a snapshot now.
    pub(crate) fn publish(&self, sim: &Simulation) {
        snapshot(&mut self.lock(), sim);
    }

    /// Follow the osculating elements of the particle called `name` on the page (None stops), as
    /// clicking it there does.
    pub fn track(&self, name: Option<&str>) {
        set_tracked(&mut self.lock(), name.map(str::to_string));
    }

    /// Pause or resume the simulation, as the page's Pause button does.
    pub fn set_paused(&self, paused: bool) {
        let mut c = self.control();
        c.paused = paused;
        self.shared.attention.store(c.paused || !c.commands.is_empty(), Ordering::Release);
        drop(c);
        self.shared.wake.notify_all();
    }

    /// Measure progress from `t_start` to `t_target` (days since the simulation's reference
    /// epoch) instead of over each call's own span, e.g. across a loop of
    /// [`Simulation::integrate_or_interpolate`] calls.
    fn set_target(&self, t_start: f64, t_target: f64) {
        self.lock().fixed_target = Some((t_start, t_target));
    }

    /// The page moved the run's end to `t_target` (days since the reference epoch; `None` runs
    /// on until stopped).
    pub(crate) fn set_end(&self, t_target: Option<f64>) {
        self.lock().t_target = t_target;
    }

    /// The page's Stop ended the run.
    pub(crate) fn note_stop(&self) {
        let mut s = self.lock();
        s.stop_asked = true;
        s.t_target = None;
    }

    /// A run (one call of `integrate`, `integrate_or_interpolate`, `run` or `steps`) starts.
    /// `editable_end` says whether the run's loop follows an end set on the page (`steps` doesn't).
    pub(crate) fn begin_run(&self, sim: &Simulation, t_target: Option<f64>, editable_end: bool) {
        let mut s = self.lock();
        let (t_start, t_target) = match s.fixed_target {
            Some((start, target)) => (start, Some(target)),
            None => (sim.t(), t_target),
        };
        if s.state != "running" || s.fixed_target.is_none() {
            s.run_started = Some(Instant::now());
            s.paused_total = Duration::ZERO;
            s.steps = 0;
        }
        s.t_start = t_start;
        s.t_target = t_target;
        s.end_editable = s.fixed_target.is_none() && editable_end;
        s.stop_asked = false;
        s.state = "running";
        snapshot(&mut s, sim);
    }

    /// `n` steps were taken; publish a snapshot if the last one is old enough.
    pub(crate) fn tick(&self, sim: &Simulation, n: usize) {
        let mut s = self.lock();
        s.steps += n as u64;
        s.total_steps += n as u64;
        if s.last_publish.map_or(true, |last| last.elapsed() >= PUBLISH_INTERVAL) || self.shared.attention.load(Ordering::Acquire) {
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
        if sim.interrupted() || s.stop_asked {
            s.state = "stopped";
        } else if finished {
            s.state = "done";
        }
        s.end_editable = false;
        if let Some(since) = s.paused_since.take() {
            s.paused_total += since.elapsed();
        }
        snapshot(&mut s, sim);
    }

    /// The energy the page measures its error from changed (a perturber was turned on or off):
    /// start the error over.
    pub(crate) fn reset_energy(&self) {
        let mut s = self.lock();
        s.energy0 = 0.0;
        s.history.clear();
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
    let paused = s.paused_total + s.paused_since.map_or(Duration::ZERO, |since| since.elapsed());
    s.elapsed = s.run_started.map_or(0.0, |start| start.elapsed().saturating_sub(paused).as_secs_f64());
    s.n_particles = sim.len();
    s.integrator = integrator_name(sim);
    s.reference_plane = sim.reference_plane.to_string();
    s.origin = sim.origin.to_string();

    let masses = sim.masses();
    let states = sim.states();
    let names = sim.names();
    s.n_massive = masses.iter().filter(|&&m| m != 0.0).count();
    let n_test = states.len() - s.n_massive;
    let stride = n_test.div_ceil(MAX_DRAWN).max(1);
    let perturbers = sim.perturbers();
    let muted = |i: usize| masses[i] == 0.0 && perturbers.iter().any(|(name, on)| !on && *name == names[i]);
    let massive: Vec<usize> = (0..states.len()).filter(|&i| masses[i] != 0.0 || muted(i)).collect();
    let test: Vec<usize> = (0..states.len()).filter(|&i| masses[i] == 0.0 && !muted(i)).step_by(stride).collect();
    let row = |&i: &usize| (names[i].clone(), states[i][0], states[i][1], states[i][2]);
    s.massive = massive.iter().map(row).collect();
    s.test_particles = test.iter().map(row).collect();
    s.drawn = massive.iter().chain(&test).map(|&i| (states[i], masses[i])).collect();
    s.central = massive.first().map(|&i| (names[i].clone(), states[i], masses[i]));
    s.perturbers = perturbers;
    s.solar_gr = sim.solar_gr();

    if let Some(name) = s.tracked.clone() {
        let found = sim.particle_index_map.get(&name).map(|&i| (states[i], masses[i]));
        record_elements(s, &name, found, sim.t());
    }
    s.last_publish = Some(Instant::now());
}

/// Work out the elements of `name` (with its state and mass, if it exists) at `t` and add them to
/// the history.
fn record_elements(s: &mut Status, name: &str, found: Option<(State, f64)>, t: f64) {
    s.elements = match (found, &s.central) {
        (Some((st, m)), Some((central, cs, cm))) if central != name => {
            let rel: State = std::array::from_fn(|k| st[k] - cs[k]);
            let mu = GRAVITATIONAL_CONSTANT * (cm + m);
            let el = state::elements(&rel, mu);
            let bound = el.e < 1.0 && el.a > 0.0;
            let deg = f64::to_degrees;
            Some(TrackedElements {
                name: name.to_string(),
                central: central.clone(),
                a: el.a,
                e: el.e,
                inc: deg(el.inc),
                node: deg(el.node),
                arg: deg(el.arg),
                mean_anomaly: deg(el.mean_anomaly),
                true_anomaly: deg(el.true_anomaly),
                q: el.q,
                big_q: bound.then(|| el.a * (1.0 + el.e)),
                period: bound.then(|| std::f64::consts::TAU * (el.a.powi(3) / mu).sqrt()),
                r: state::position(&rel).norm(),
            })
        }
        _ => None,
    };
    if let Some(el) = &s.elements {
        let point = [t, el.a, el.e, el.inc, el.q];
        if s.element_history.last().map_or(true, |p| p[0] != t) {
            if s.element_history.len() == MAX_HISTORY {
                let thinned = s.element_history.iter().step_by(2).copied().collect();
                s.element_history = thinned;
            }
            s.element_history.push(point);
        }
    }
}

/// Follow the elements of the particle called `name` (None stops). Uses the last snapshot when
/// the particle is drawn on the page, so the page shows them at once even while nothing steps.
fn set_tracked(s: &mut Status, name: Option<String>) {
    s.element_history.clear();
    s.elements = None;
    s.tracked = name.clone();
    if let Some(name) = name {
        let found = s
            .massive
            .iter()
            .chain(&s.test_particles)
            .position(|p| p.0 == name)
            .map(|i| s.drawn[i]);
        if found.is_some() {
            let t = s.t;
            record_elements(s, &name, found, t);
        }
    }
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

fn respond(request: &Request, shared: &Shared) -> Response {
    match (request.method.as_str(), request.path.as_str()) {
        ("GET", "/" | "/index.html") => Response::html(PAGE),
        ("GET", "/progress.json") => {
            let mut status = shared.status.lock().unwrap_or_else(|e| e.into_inner()).clone();
            status.paused = shared.control.lock().unwrap_or_else(|e| e.into_inner()).paused;
            Response::json(&status)
        }
        ("POST", "/control") => Response::outcome(control(shared, request)),
        _ => Response::not_found(),
    }
}

/// Apply `action=pause|resume|step|back|stop|integrator|timestep|end|gr|perturber|track` (with
/// `name=`, `timestep=`, `jd=` (none for no end) and `on=1|0`).
fn control(shared: &Shared, request: &Request) -> Result<(), String> {
    let param = |key: &str| request.param(key);
    if param("action").as_deref() == Some("track") {
        let name = param("name").filter(|n| !n.is_empty());
        set_tracked(&mut shared.status.lock().unwrap_or_else(|e| e.into_inner()), name);
        return Ok(());
    }
    let timestep = || -> Result<f64, String> {
        let h: f64 = param("timestep").ok_or("missing timestep")?.parse().map_err(|_| "bad timestep")?;
        if h.is_finite() && h != 0.0 { Ok(h) } else { Err("timestep must be finite and non-zero".into()) }
    };
    let on = || match param("on").as_deref() {
        Some("1" | "true") => Ok(true),
        Some("0" | "false") => Ok(false),
        _ => Err("on must be 1 or 0".to_string()),
    };
    let jd_ref = shared.status.lock().unwrap_or_else(|e| e.into_inner()).jd_ref;
    let mut c = shared.control.lock().unwrap_or_else(|e| e.into_inner());
    match param("action").as_deref().unwrap_or("") {
        "pause" => c.paused = true,
        "resume" => {
            c.paused = false;
            c.step_credits = 0;
            c.back_credits = 0;
        }
        "step" => {
            c.paused = true;
            c.step_credits += 1;
        }
        "back" => {
            c.paused = true;
            c.back_credits += 1;
        }
        "end" => {
            let end = match param("jd").as_deref() {
                None | Some("") => None,
                Some(jd) => {
                    let jd: f64 = jd.parse().map_err(|_| "bad end date")?;
                    if !jd.is_finite() {
                        return Err("end date must be finite".into());
                    }
                    Some(jd - jd_ref)
                }
            };
            c.commands.push(Command::End(end));
        }
        "stop" => {
            c.paused = false;
            c.step_credits = 0;
            c.back_credits = 0;
            c.commands.push(Command::Stop);
        }
        "gr" => c.commands.push(Command::SolarGR(on()?)),
        "perturber" => {
            let name = param("name").ok_or("missing name")?;
            c.commands.push(Command::Perturber(name, on()?));
        }
        "integrator" => {
            let name = param("name").ok_or("missing name")?;
            if !INTEGRATORS.contains(&name.as_str()) {
                return Err(format!("unknown integrator {name}"));
            }
            c.commands.push(Command::Integrator(name, timestep()?));
        }
        "timestep" => c.commands.push(Command::Timestep(timestep()?)),
        other => return Err(format!("unknown action {other}")),
    }
    shared.attention.store(c.paused || c.back_credits > 0 || !c.commands.is_empty(), Ordering::Release);
    drop(c);
    shared.wake.notify_all();
    Ok(())
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

    /// The attached progress page, e.g. to pause it from code ([`ProgressMonitor::set_paused`]).
    pub fn progress_monitor(&self) -> Option<&ProgressMonitor> {
        self.progress_monitor.as_ref()
    }

    /// Detach the progress page; its server stops once no clone of this simulation holds it.
    pub fn stop_monitor(&mut self) {
        self.progress_monitor = None;
    }
}

const PAGE: &str = include_str!("monitor.html");

