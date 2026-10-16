//! A live orbit-fitting dashboard: a localhost web page to fit an orbit interactively.
//!
//! [`FitDashboard::start`] takes the detections, a kernel and starting [`FitOptions`], serves
//! the page on `127.0.0.1` (in the same 8-bit style as the simulation monitor) and fits on a
//! background thread. From the page:
//!
//! - **Forces**: the 16 massive asteroids, each planet, the Moon and Pluto, Earth J2–J4, solar
//!   J2, the relativistic model (EIH, simple, potential or off), the non-gravitational
//!   parameters A1, A2, A3 (or the automatic choice) and the comet g(r).
//! - **Solver**: the engine (Cartesian or Bernstein–Khushalani), the initial orbit (Gauss with
//!   the BK fallback, Gauss only, or Herget), robust outlier rejection and its threshold, the
//!   IAS15 tolerance, the iteration limit, the convergence tolerance, the chi-square threshold
//!   and the arc gap.
//! - **Fit** runs the pipeline (from scratch, or from the last orbit with "start from last
//!   fit"); **Stop** ends a fit early; with **Auto refit** on, every change refits at once.
//!   Clicking a detection on the residual plot leaves it out of the fit (or puts it back).
//! - **Done** ends [`FitDashboard::wait`], which returns the last fit.
//! - **Load** fetches the observations of the object typed beside it from the Minor Planet
//!   Center, and **Random** those of a random numbered asteroid (see [`crate::orbfit::mpc`]),
//!   then fits them. [`FitDashboard::set_loader`] swaps in another source.
//!
//! The page shows each Levenberg–Marquardt iteration's chi-square while a fit runs, then the
//! residuals of every detection (in arcseconds, against time and against each other), the
//! heliocentric ecliptic elements with their 1-sigma uncertainties, a top-down map of the orbit
//! and the planets, and a log of the fits so far.
//!
//! The page changes things with `POST /control?action=...`, which needs an `X-Spacerocks` header
//! and a local `Host`, as the monitor's does.
use std::sync::{Arc, Condvar, Mutex, MutexGuard};
use std::time::{Duration, Instant};

use serde_json::{json, Value};

use crate::constants::GRAVITATIONAL_CONSTANT;
use crate::orbfit::{
    determine_orbit, fetch_astrometry, residuals, Astrometry, Engine, FitIteration, FitOptions, FitProgress, ForceModel, NongravAuto, OrbitFit, Relativity,
    ASTEROID_PERTURBERS,
};
use crate::spice::SpiceKernel;
use crate::state;
use crate::web::{self, Request, Response};

const PAGE: &str = include_str!("dashboard.html");

/// Radians to arcseconds.
const ARCSEC: f64 = 206264.80624709636;

/// Most iterations kept for the live chi-square plot.
const MAX_LIVE: usize = 2000;

/// Most fits kept in the log.
const MAX_LOG: usize = 12;

/// The major bodies the page can switch off: (NAIF id, name).
const MAJOR: [(i32, &str); 10] =
    [(1, "Mercury"), (2, "Venus"), (399, "Earth"), (301, "Moon"), (4, "Mars"), (5, "Jupiter"), (6, "Saturn"), (7, "Uranus"), (8, "Neptune"), (9, "Pluto")];

/// The planets drawn on the map (barycenters, as the kernel has them).
const MAP_PLANETS: [(i32, &str); 8] = [(1, "Mercury"), (2, "Venus"), (399, "Earth"), (4, "Mars"), (5, "Jupiter"), (6, "Saturn"), (7, "Uranus"), (8, "Neptune")];

/// The comet g(r) of Marsden et al. (1973), as `FitOptions::gofr` takes it: `[alpha, m, n, k, r0]`.
const COMET_GOFR: [f64; 5] = [0.1112620426, 2.15, 5.093, 4.6142, 2.808];

/// Most detections the default loader keeps of an object (spread over its arc).
pub const LOAD_MAX_DETECTIONS: usize = 1000;

/// Where Load and Random get detections: a designation (`None` for a random object) and the
/// kernel to place the observers with, to the object's label and its astrometry.
pub type Loader = Arc<dyn Fn(Option<&str>, &SpiceKernel) -> Result<(String, Astrometry), String> + Send + Sync>;

/// A running dashboard. Dropping it stops the fitting thread and the page's server.
pub struct FitDashboard {
    shared: Arc<Shared>,
    url: String,
}

struct Shared {
    status: Mutex<Status>,
    /// Signalled when a fit is asked for, a fit ends, Done is pressed or the dashboard closes.
    wake: Condvar,
}

struct Status {
    name: String,
    /// Where Load and Random get detections.
    loader: Loader,
    /// A load asked for: the designation, empty for a random object.
    load: Option<String>,
    astrometry: Astrometry,
    options: FitOptions,
    /// An initial orbit (epoch, state, nongrav) to start the first fit from, if given.
    initial: Option<(f64, [f64; 6], [f64; 3])>,
    /// "idle", "fitting", "done" (a fit finished), "stopped" (one was stopped) or "loading".
    state: &'static str,
    /// A fit has been asked for.
    pending: bool,
    /// The running fit should stop.
    cancel: bool,
    closed: bool,
    /// Done was pressed.
    done: bool,
    /// Refit whenever something changes.
    auto: bool,
    /// Start from the last fit's orbit rather than an initial orbit determination.
    seed_last: bool,
    /// The settings changed since the last fit started.
    changed: bool,
    /// Detections left out by hand.
    excluded: Vec<bool>,
    /// `[iteration, chi2, detections, npar, accepted]` of the running (or last) fit.
    live: Vec<[f64; 5]>,
    started: Option<Instant>,
    elapsed: f64,
    runs: usize,
    /// The last fit, with per-detection outputs over every detection.
    fit: Option<OrbitFit>,
    /// What the page shows of it.
    fit_view: Value,
    log: Vec<Value>,
    /// Static page data: detections and the heliocentric ecliptic observer positions.
    detections: Value,
    observers: Value,
    message: String,
}

impl FitDashboard {
    /// Serve the dashboard for `astrometry` on `127.0.0.1:port` (`port` 0 picks a free port),
    /// fitting with `kernel` from the settings in `options`, and start a first fit at once
    /// (from `initial`, `(epoch, state, nongrav)`, if given). `name` labels the page.
    pub fn start(port: u16, name: &str, astrometry: Astrometry, kernel: SpiceKernel, options: FitOptions, initial: Option<(f64, [f64; 6], [f64; 3])>) -> std::io::Result<FitDashboard> {
        let n = astrometry.len();
        let (detections, observers) = views(&astrometry, &kernel);
        let shared = Arc::new(Shared {
            status: Mutex::new(Status {
                name: name.to_string(),
                loader: Arc::new(|designation, kernel| fetch_astrometry(designation, kernel, LOAD_MAX_DETECTIONS)),
                load: None,
                astrometry,
                options: FitOptions { progress: None, ..options },
                initial,
                state: "idle",
                pending: true,
                cancel: false,
                closed: false,
                done: false,
                auto: true,
                seed_last: false,
                changed: false,
                excluded: vec![false; n],
                live: Vec::new(),
                started: None,
                elapsed: 0.0,
                runs: 0,
                fit: None,
                fit_view: Value::Null,
                log: Vec::new(),
                detections,
                observers,
                message: String::new(),
            }),
            wake: Condvar::new(),
        });
        let url = web::start(port, "spacerocks-fit-dashboard", Arc::downgrade(&shared), respond)?;
        let worker = shared.clone();
        std::thread::Builder::new().name("spacerocks-fit".into()).spawn(move || work(&worker, &kernel))?;
        Ok(FitDashboard { shared, url })
    }

    /// The page's address, like `http://127.0.0.1:8765/`.
    pub fn url(&self) -> &str {
        &self.url
    }

    fn lock(&self) -> MutexGuard<'_, Status> {
        self.shared.lock()
    }

    /// The last fit, if one has finished. Its residuals and `used` cover every detection;
    /// detections left out on the page are marked unused.
    pub fn fit(&self) -> Option<OrbitFit> {
        self.lock().fit.clone()
    }

    /// The settings as they are on the page now.
    pub fn options(&self) -> FitOptions {
        self.lock().options.clone()
    }

    /// What the page shows of the last fit, as JSON (`Null` before one finishes): its outcome,
    /// the residuals in arcseconds and their RMS, and the heliocentric ecliptic elements with
    /// their 1-sigma uncertainties (`elements`, angles in degrees).
    pub fn fit_summary(&self) -> Value {
        self.lock().fit_view.clone()
    }

    /// Which detections are left out by hand.
    pub fn excluded(&self) -> Vec<bool> {
        self.lock().excluded.clone()
    }

    /// Whether a fit is running or about to.
    pub fn fitting(&self) -> bool {
        let s = self.lock();
        s.pending || s.state == "fitting"
    }

    /// Whether another object's detections are being fetched or about to be.
    pub fn loading(&self) -> bool {
        let s = self.lock();
        s.load.is_some() || s.state == "loading"
    }

    /// The object's label, which changes when another object is loaded.
    pub fn name(&self) -> String {
        self.lock().name.clone()
    }

    /// The detections being fitted, which change when another object is loaded.
    pub fn astrometry(&self) -> Astrometry {
        self.lock().astrometry.clone()
    }

    /// Fetch another object's detections and fit them, as the page's Load (with a designation)
    /// and Random (`None`) buttons do. Stops a running fit. The fetch runs on the fitting
    /// thread; if it fails, the page shows why and keeps the current object.
    pub fn load(&self, designation: Option<&str>) {
        let mut s = self.lock();
        request_load(&mut s, designation.unwrap_or("").trim());
        drop(s);
        self.shared.wake.notify_all();
    }

    /// Get Load and Random's detections from `loader` instead of the Minor Planet Center.
    pub fn set_loader(&self, loader: impl Fn(Option<&str>, &SpiceKernel) -> Result<(String, Astrometry), String> + Send + Sync + 'static) {
        self.lock().loader = Arc::new(loader);
    }

    /// Ask for a fit, as the page's Fit button does.
    pub fn refit(&self) {
        let mut s = self.lock();
        s.pending = true;
        drop(s);
        self.shared.wake.notify_all();
    }

    /// Wait until Done is pressed on the page, then until no fit is running, and return the
    /// last fit. `interrupted` is checked every 50 ms; when it returns true, waiting stops (the
    /// dashboard keeps running) and this returns `None`.
    pub fn wait(&self, interrupted: Option<&dyn Fn() -> bool>) -> Option<OrbitFit> {
        let mut s = self.lock();
        while !(s.done && !s.pending && s.load.is_none() && !matches!(s.state, "fitting" | "loading")) && !s.closed {
            if interrupted.is_some_and(|f| f()) {
                return None;
            }
            s = self.shared.wake.wait_timeout(s, Duration::from_millis(50)).unwrap_or_else(|e| e.into_inner()).0;
        }
        s.fit.clone()
    }

    /// Wait (up to `timeout`) until no fit or load is running or asked for; true if none is.
    pub fn wait_idle(&self, timeout: Duration) -> bool {
        let end = Instant::now() + timeout;
        let mut s = self.lock();
        while s.pending || s.load.is_some() || matches!(s.state, "fitting" | "loading") {
            let now = Instant::now();
            if now >= end {
                return false;
            }
            s = self.shared.wake.wait_timeout(s, (end - now).min(Duration::from_millis(50))).unwrap_or_else(|e| e.into_inner()).0;
        }
        true
    }
}

impl Drop for FitDashboard {
    fn drop(&mut self) {
        let mut s = self.lock();
        s.closed = true;
        s.cancel = true;
        drop(s);
        self.shared.wake.notify_all();
    }
}

impl Shared {
    fn lock(&self) -> MutexGuard<'_, Status> {
        self.status.lock().unwrap_or_else(|e| e.into_inner())
    }
}

/// What the page shows of the detections: `[jd, ra°, dec°, sigma_ra″, sigma_dec″]` each, and
/// the heliocentric ecliptic positions of their observers.
fn views(a: &Astrometry, kernel: &SpiceKernel) -> (Value, Value) {
    let n = a.len();
    let detections = json!((0..n).map(|i| json!([a.epoch[i], a.ra[i].to_degrees(), a.dec[i].to_degrees(), a.sigma_ra[i] * ARCSEC, a.sigma_dec[i] * ARCSEC])).collect::<Vec<_>>());
    let observers = json!((0..n)
        .map(|i| {
            let sun = kernel.state_au(10, 0, a.epoch[i]).unwrap_or([f64::NAN; 6]);
            let o = a.observer[i];
            ecliptic(&[o[0] - sun[0], o[1] - sun[1], o[2] - sun[2]]).to_vec()
        })
        .collect::<Vec<_>>());
    (detections, observers)
}

/// Ask the fitting thread to load `designation` (empty for a random object), stopping a
/// running fit.
fn request_load(s: &mut Status, designation: &str) {
    if s.state == "fitting" {
        s.cancel = true;
    }
    s.load = Some(designation.to_string());
    s.done = false;
}

/// Fetch the detections asked for and swap them in (outside the lock), or say why not.
fn load(shared: &Arc<Shared>, kernel: &SpiceKernel) {
    let mut s = shared.lock();
    let Some(designation) = s.load.take() else { return };
    let before = s.state;
    s.state = "loading";
    s.message = if designation.is_empty() { "fetching a random object from the MPC...".into() } else { format!("fetching {designation} from the MPC...") };
    let loader = s.loader.clone();
    drop(s);
    shared.wake.notify_all();

    let loaded = loader((!designation.is_empty()).then_some(designation.as_str()), kernel).and_then(|(name, a)| {
        if a.len() < 3 {
            return Err(format!("{name} has only {} usable detections", a.len()));
        }
        let (detections, observers) = views(&a, kernel);
        Ok((name, a, detections, observers))
    });

    let mut s = shared.lock();
    match loaded {
        Ok((name, a, detections, observers)) => {
            s.message = format!("loaded {name}: {} detections", a.len());
            s.name = name;
            s.excluded = vec![false; a.len()];
            s.astrometry = a;
            s.detections = detections;
            s.observers = observers;
            s.initial = None;
            s.fit = None;
            s.fit_view = Value::Null;
            s.log.clear();
            s.live.clear();
            s.runs = 0;
            s.elapsed = 0.0;
            s.changed = false;
            s.state = "idle";
            s.pending = true;
        }
        Err(e) => {
            s.message = format!("couldn't load: {e}");
            s.state = if before == "fitting" { "stopped" } else { before };
        }
    }
    drop(s);
    shared.wake.notify_all();
}

/// J2000 equatorial to ecliptic (obliquity 84381.448").
fn ecliptic(v: &[f64; 3]) -> [f64; 3] {
    let eps = (84381.448f64 / 3600.0).to_radians();
    let (s, c) = eps.sin_cos();
    [v[0], c * v[1] + s * v[2], -s * v[1] + c * v[2]]
}

/// The fitting thread: run a fit whenever one is asked for, until the dashboard is dropped.
fn work(shared: &Arc<Shared>, kernel: &SpiceKernel) {
    loop {
        let mut s = shared.lock();
        while !s.pending && s.load.is_none() && !s.closed {
            s = shared.wake.wait(s).unwrap_or_else(|e| e.into_inner());
        }
        if s.closed {
            return;
        }
        if s.load.is_some() {
            drop(s);
            load(shared, kernel);
            continue;
        }
        s.pending = false;
        s.cancel = false;
        s.changed = false;
        s.state = "fitting";
        s.live.clear();
        s.started = Some(Instant::now());
        let all = s.astrometry.clone();
        let excluded = s.excluded.clone();
        let mut options = s.options.clone();
        let initial = match (&s.fit, s.seed_last) {
            (Some(f), true) if f.state.iter().all(|v| v.is_finite()) => Some((f.epoch, f.state, f.nongrav)),
            _ if s.runs == 0 => s.initial,
            _ => None,
        };
        drop(s);

        let progress = Arc::downgrade(shared);
        options.progress = Some(FitProgress::new(move |it: &FitIteration| {
            let Some(shared) = progress.upgrade() else { return false };
            let mut s = shared.lock();
            if s.live.len() == MAX_LIVE {
                let thinned = s.live.iter().step_by(2).copied().collect();
                s.live = thinned;
            }
            let k = s.live.len() as f64;
            s.live.push([k, it.chi2, it.detections as f64, it.npar as f64, if it.accepted { 1.0 } else { 0.0 }]);
            !s.cancel
        }));
        let idx: Vec<usize> = (0..all.len()).filter(|&i| !excluded[i]).collect();
        let sub = all.subset(&idx);
        let fit = determine_orbit(&sub, initial, kernel, &options);
        let full = expand(&all, &idx, fit, kernel, &options);
        let view = fit_view(&full, &all, &excluded, kernel);

        let mut s = shared.lock();
        s.elapsed = s.started.map_or(0.0, |t| t.elapsed().as_secs_f64());
        s.started = None;
        s.state = if s.cancel { "stopped" } else { "done" };
        s.runs += 1;
        let entry = json!({
            "run": s.runs,
            "status": if s.cancel { "stopped" } else { full.flag.description() },
            "converged": full.converged(),
            "chi2_nu": if full.ndof > 0 { full.chi2 / full.ndof as f64 } else { f64::NAN },
            "rms": view["rms"],
            "used": full.used.iter().filter(|&&u| u).count(),
            "settings": summary(&options, idx.len(), all.len()),
            "seconds": s.elapsed,
        });
        s.log.insert(0, entry);
        s.log.truncate(MAX_LOG);
        s.fit = Some(full);
        s.fit_view = view;
        drop(s);
        shared.wake.notify_all();
    }
}

/// Put a fit to the detections `idx` of `all` back in terms of every detection: residuals at
/// the final orbit for all of them (the left-out ones included) and `used` false for those left
/// out.
fn expand(all: &Astrometry, idx: &[usize], mut fit: OrbitFit, kernel: &SpiceKernel, options: &FitOptions) -> OrbitFit {
    let n = all.len();
    let mut used = vec![false; n];
    for (k, &i) in idx.iter().enumerate() {
        used[i] = fit.used.get(k).copied().unwrap_or(false);
    }
    let mut resid = vec![[f64::NAN; 6]; n];
    let quiet = FitOptions { progress: None, ..options.clone() };
    let at_orbit = fit.state.iter().all(|v| v.is_finite()) && fit.epoch.is_finite();
    match at_orbit.then(|| residuals(all, fit.epoch, &fit.state, &fit.nongrav, kernel, &quiet, false)) {
        Some(Ok(r)) => resid = r.per_detection(n),
        _ => {
            for (k, &i) in idx.iter().enumerate() {
                if let Some(r) = fit.residuals.get(k) {
                    resid[i] = *r;
                }
            }
        }
    }
    fit.residuals = resid;
    fit.used = if fit.used.is_empty() { Vec::new() } else { used };
    fit
}

/// A short description of the settings of a fit, for the log.
fn summary(o: &FitOptions, used: usize, n: usize) -> String {
    let mut parts = Vec::new();
    let f = &o.forces;
    let asteroids = ASTEROID_PERTURBERS.iter().all(|c| !f.exclude.contains(c));
    parts.push(if asteroids { "asteroids".to_string() } else { "no asteroids".to_string() });
    let off: Vec<&str> = MAJOR.iter().filter(|(c, _)| f.exclude.contains(c)).map(|(_, n)| *n).collect();
    if !off.is_empty() {
        parts.push(format!("no {}", off.join("/")));
    }
    if !f.earth_harmonics {
        parts.push("no earth J2".into());
    }
    if !f.solar_j2 {
        parts.push("no solar J2".into());
    }
    parts.push(format!("gr {}", relativity_name(f.relativity)));
    if o.nongrav_auto.is_some() {
        parts.push("nongrav auto".into());
    } else if o.fit_nongrav.iter().any(|&b| b) {
        let names: Vec<String> = (0..3).filter(|&k| o.fit_nongrav[k]).map(|k| format!("A{}", k + 1)).collect();
        parts.push(names.join(""));
    }
    if o.gofr.is_some() {
        parts.push("comet g(r)".into());
    }
    parts.push(if o.engine == Engine::BkNative { "bk".into() } else { "cartesian".into() });
    parts.push(iod_name(o).into());
    if o.robust {
        parts.push(format!("robust {}σ", o.outlier_sigma));
    }
    if used < n {
        parts.push(format!("{} left out", n - used));
    }
    parts.join(" · ")
}

fn relativity_name(r: Relativity) -> &'static str {
    match r {
        Relativity::Eih => "eih",
        Relativity::Simple => "simple",
        Relativity::Potential => "potential",
        Relativity::Off => "off",
    }
}

fn iod_name(o: &FitOptions) -> &'static str {
    if o.herget {
        "herget"
    } else if o.bk_fallback {
        "auto"
    } else {
        "gauss"
    }
}

/// Wrap an angle difference into (-pi, pi].
fn wrap(d: f64) -> f64 {
    let t = std::f64::consts::TAU;
    let w = d.rem_euclid(t);
    if w > std::f64::consts::PI { w - t } else { w }
}

/// Heliocentric ecliptic elements `[a, e, inc, node, arg, mean anomaly, q]` (radians) of the
/// heliocentric ecliptic state `s`.
fn element_vector(s: &[f64; 6]) -> [f64; 7] {
    let el = state::elements(s, GRAVITATIONAL_CONSTANT);
    [el.a, el.e, el.inc, el.node, el.arg, el.mean_anomaly, el.q]
}

/// What the page shows of a fit: its outcome, residuals in arcseconds, and its heliocentric
/// ecliptic elements with uncertainties.
fn fit_view(fit: &OrbitFit, all: &Astrometry, excluded: &[bool], kernel: &SpiceKernel) -> Value {
    let n = all.len();
    let resid: Vec<Value> = (0..n).map(|i| json!([fit.residuals[i][0] * ARCSEC, fit.residuals[i][1] * ARCSEC])).collect();
    // A fit with no orbit marks nothing as used: count every detection it was given.
    let used: Vec<u8> = (0..n).map(|i| if fit.used.is_empty() { !excluded[i] as u8 } else { fit.used[i] as u8 }).collect();
    let (mut sum, mut count) = ([0.0; 2], 0usize);
    for i in 0..n {
        let r = fit.residuals[i];
        if used[i] == 1 && r[0].is_finite() && r[1].is_finite() {
            sum[0] += (r[0] * ARCSEC).powi(2);
            sum[1] += (r[1] * ARCSEC).powi(2);
            count += 1;
        }
    }
    let rms = if count > 0 { json!([(sum[0] / count as f64).sqrt(), (sum[1] / count as f64).sqrt(), ((sum[0] + sum[1]) / count as f64).sqrt()]) } else { json!([f64::NAN, f64::NAN, f64::NAN]) };

    let has_orbit = fit.state.iter().all(|v| v.is_finite()) && fit.epoch.is_finite();
    let mut elements = Value::Null;
    let mut position = Value::Null;
    let mut planets = Value::Null;
    if has_orbit {
        if let Ok(sun) = kernel.state_au(10, 0, fit.epoch) {
            let rel = std::array::from_fn::<f64, 3, _>(|k| fit.state[k] - sun[k]);
            let vel = std::array::from_fn::<f64, 3, _>(|k| fit.state[k + 3] - sun[k + 3]);
            let (r, v) = (ecliptic(&rel), ecliptic(&vel));
            let s: [f64; 6] = [r[0], r[1], r[2], v[0], v[1], v[2]];
            let el = element_vector(&s);
            // Linear propagation of the state covariance: rotate it to the ecliptic, then
            // through the elements' Jacobian (central differences).
            let cov = fit.state_covariance();
            let rot = |m: &[[f64; 6]; 6]| -> [[f64; 6]; 6] {
                let eps = (84381.448f64 / 3600.0).to_radians();
                let (se, ce) = eps.sin_cos();
                let mut rm = [[0.0; 6]; 6];
                for b in [0, 3] {
                    rm[b][b] = 1.0;
                    rm[b + 1][b + 1] = ce;
                    rm[b + 1][b + 2] = se;
                    rm[b + 2][b + 1] = -se;
                    rm[b + 2][b + 2] = ce;
                }
                let mut out = [[0.0; 6]; 6];
                for i in 0..6 {
                    for j in 0..6 {
                        out[i][j] = (0..6).map(|k| (0..6).map(|l| rm[i][k] * m[k][l] * rm[j][l]).sum::<f64>()).sum();
                    }
                }
                out
            };
            let c = rot(&cov);
            let mut jac = [[0.0; 6]; 7];
            for j in 0..6 {
                let h = if j < 3 { 1e-7 } else { 1e-9 } * s[j].abs().max(if j < 3 { 1.0 } else { 0.01 });
                let (mut p, mut m) = (s, s);
                p[j] += h;
                m[j] -= h;
                let (ep, em) = (element_vector(&p), element_vector(&m));
                for i in 0..7 {
                    let d = if (2..6).contains(&i) { wrap(ep[i] - em[i]) } else { ep[i] - em[i] };
                    jac[i][j] = d / (2.0 * h);
                }
            }
            let sigma: Vec<f64> = (0..7)
                .map(|i| (0..6).map(|k| (0..6).map(|l| jac[i][k] * c[k][l] * jac[i][l]).sum::<f64>()).sum::<f64>().max(0.0).sqrt())
                .collect();
            let deg = |x: f64| x.to_degrees();
            let bound = el[1] < 1.0 && el[0] > 0.0;
            let ng_sigma = fit.nongrav_sigma();
            elements = json!({
                "a": el[0], "e": el[1], "inc": deg(el[2]), "node": deg(el[3].rem_euclid(std::f64::consts::TAU)),
                "arg": deg(el[4].rem_euclid(std::f64::consts::TAU)), "mean_anomaly": deg(el[5].rem_euclid(std::f64::consts::TAU)), "q": el[6],
                "sigma": [sigma[0], sigma[1], deg(sigma[2]), deg(sigma[3]), deg(sigma[4]), deg(sigma[5]), sigma[6]],
                "big_q": if bound { el[0] * (1.0 + el[1]) } else { f64::NAN },
                "period": if bound { std::f64::consts::TAU * (el[0].powi(3) / GRAVITATIONAL_CONSTANT).sqrt() } else { f64::NAN },
                "nongrav": fit.nongrav, "nongrav_sigma": ng_sigma, "fit_nongrav": fit.fit_nongrav,
                "n": (GRAVITATIONAL_CONSTANT / el[0].abs().powi(3)).sqrt().to_degrees(),
            });
            position = json!(r);
        }
        planets = json!(MAP_PLANETS
            .iter()
            .filter_map(|&(code, name)| kernel.state_au(code, 10, fit.epoch).ok().map(|p| {
                let e = ecliptic(&[p[0], p[1], p[2]]);
                json!([name, e[0], e[1], e[2]])
            }))
            .collect::<Vec<_>>());
    }
    json!({
        "flag": fit.flag.code(),
        "status": fit.flag.description(),
        "converged": fit.converged(),
        "chi2": fit.chi2,
        "ndof": fit.ndof,
        "chi2_nu": if fit.ndof > 0 { fit.chi2 / fit.ndof as f64 } else { f64::NAN },
        "niter": fit.niter,
        "epoch": fit.epoch,
        "state": fit.state,
        "resid": resid,
        "used": used,
        "rms": rms,
        "elements": elements,
        "position": position,
        "planets": planets,
    })
}

/// The page's settings, from the options.
fn settings(o: &FitOptions) -> Value {
    let f = &o.forces;
    json!({
        "asteroids": ASTEROID_PERTURBERS.iter().all(|c| !f.exclude.contains(c)),
        "bodies": MAJOR.iter().map(|(c, name)| json!([name, !f.exclude.contains(c)])).collect::<Vec<_>>(),
        "earth_harmonics": f.earth_harmonics,
        "solar_j2": f.solar_j2,
        "relativity": relativity_name(f.relativity),
        "nongrav": o.fit_nongrav,
        "nongrav_auto": o.nongrav_auto.is_some(),
        "comet": o.gofr.is_some(),
        "engine": if o.engine == Engine::BkNative { "bk" } else { "cartesian" },
        "iod": iod_name(o),
        "robust": o.robust,
        "outlier_sigma": o.outlier_sigma,
        "epsilon": o.epsilon,
        "max_iter": o.max_iter,
        "tolerance": o.tolerance,
        "chi2_threshold": o.chi2_threshold,
        "arc_gap": o.arc_gap,
    })
}

fn respond(request: &Request, shared: &Shared) -> Response {
    match (request.method.as_str(), request.path.as_str()) {
        ("GET", "/" | "/index.html") => Response::html(PAGE),
        ("GET", "/state.json") => {
            let s = shared.lock();
            let elapsed = s.started.map_or(s.elapsed, |t| t.elapsed().as_secs_f64());
            Response::json(&json!({
                "name": s.name,
                "state": s.state,
                "pending": s.pending,
                "loading": s.load.is_some() || s.state == "loading",
                "done": s.done,
                "auto": s.auto,
                "seed_last": s.seed_last,
                "changed": s.changed,
                "runs": s.runs,
                "elapsed": elapsed,
                "excluded": s.excluded,
                "settings": settings(&s.options),
                "detections": s.detections,
                "observers": s.observers,
                "live": s.live,
                "fit": s.fit_view,
                "log": s.log,
                "message": s.message,
            }))
        }
        ("POST", "/control") => {
            let r = control(shared, request);
            shared.wake.notify_all();
            Response::outcome(r)
        }
        _ => Response::not_found(),
    }
}

/// Apply `action=fit|stop|done|auto|seed|exclude|include_all|set|load` (with `on=1|0`,
/// `index=`, `key=`/`value=` and `designation=`, none for a random object).
fn control(shared: &Shared, request: &Request) -> Result<(), String> {
    let param = |key: &str| request.param(key);
    let on = || match param("on").as_deref() {
        Some("1" | "true") => Ok(true),
        Some("0" | "false") => Ok(false),
        _ => Err("on must be 1 or 0".to_string()),
    };
    let mut s = shared.lock();
    s.message.clear();
    match param("action").as_deref().unwrap_or("") {
        "fit" => {
            if s.state == "fitting" {
                s.cancel = true;
            }
            s.pending = true;
            s.done = false;
        }
        "stop" => {
            s.cancel = true;
            s.pending = false;
        }
        "done" => s.done = true,
        "load" => request_load(&mut s, param("designation").unwrap_or_default().trim()),
        "auto" => s.auto = on()?,
        "seed" => s.seed_last = on()?,
        "exclude" => {
            let i: usize = param("index").ok_or("missing index")?.parse().map_err(|_| "bad index")?;
            if i >= s.excluded.len() {
                return Err(format!("no detection {i}"));
            }
            let out = on()?;
            if s.excluded[i] != out {
                if out && s.excluded.iter().filter(|&&x| !x).count() <= 3 {
                    return Err("a fit needs at least 3 detections".into());
                }
                s.excluded[i] = out;
                changed(&mut s);
            }
        }
        "include_all" => {
            if s.excluded.iter().any(|&x| x) {
                s.excluded.iter_mut().for_each(|x| *x = false);
                changed(&mut s);
            }
        }
        "set" => {
            let key = param("key").ok_or("missing key")?;
            let value = param("value").ok_or("missing value")?;
            let mut o = s.options.clone();
            set(&mut o, &key, &value)?;
            if settings(&o) != settings(&s.options) {
                s.options = o;
                changed(&mut s);
            }
        }
        other => return Err(format!("unknown action {other}")),
    }
    Ok(())
}

/// The settings changed: refit at once with Auto refit on (stopping a fit that is running).
fn changed(s: &mut Status) {
    s.changed = true;
    s.done = false;
    if s.auto {
        if s.state == "fitting" {
            s.cancel = true;
        }
        s.pending = true;
    }
}

/// Change one setting from the page.
fn set(o: &mut FitOptions, key: &str, value: &str) -> Result<(), String> {
    let flag = || match value {
        "1" | "true" => Ok(true),
        "0" | "false" => Ok(false),
        _ => Err(format!("{key} must be 1 or 0")),
    };
    let number = || -> Result<f64, String> {
        let x: f64 = value.parse().map_err(|_| format!("{key}: not a number"))?;
        if x.is_finite() { Ok(x) } else { Err(format!("{key} must be finite")) }
    };
    let positive = || -> Result<f64, String> {
        let x = number()?;
        if x > 0.0 { Ok(x) } else { Err(format!("{key} must be positive")) }
    };
    let f: &mut ForceModel = &mut o.forces;
    let toggle = |exclude: &mut Vec<i32>, codes: &[i32], on: bool| {
        exclude.retain(|c| !codes.contains(c));
        if !on {
            exclude.extend_from_slice(codes);
        }
    };
    match key {
        "asteroids" => toggle(&mut f.exclude, &ASTEROID_PERTURBERS, flag()?),
        "earth_harmonics" => f.earth_harmonics = flag()?,
        "solar_j2" => f.solar_j2 = flag()?,
        "relativity" => {
            f.relativity = match value {
                "eih" => Relativity::Eih,
                "simple" => Relativity::Simple,
                "potential" => Relativity::Potential,
                "off" => Relativity::Off,
                _ => return Err(format!("unknown relativity {value}")),
            }
        }
        "a1" | "a2" | "a3" => {
            let k = key[1..].parse::<usize>().unwrap() - 1;
            o.fit_nongrav[k] = flag()?;
            if o.fit_nongrav[k] {
                o.nongrav_auto = None;
            }
            if !o.fit_nongrav.iter().any(|&b| b) {
                o.nongrav_per_arc = false;
            }
        }
        "nongrav_auto" => {
            o.nongrav_auto = if flag()? { Some(NongravAuto::default()) } else { None };
            if o.nongrav_auto.is_some() {
                o.fit_nongrav = [false; 3];
                o.nongrav_per_arc = false;
            }
        }
        "comet" => o.gofr = if flag()? { Some(COMET_GOFR) } else { None },
        "engine" => {
            o.engine = match value {
                "cartesian" => Engine::Cartesian,
                "bk" => Engine::BkNative,
                _ => return Err(format!("unknown engine {value}")),
            }
        }
        "iod" => {
            (o.bk_fallback, o.herget) = match value {
                "auto" => (true, false),
                "gauss" => (false, false),
                "herget" => (false, true),
                _ => return Err(format!("unknown initial orbit method {value}")),
            }
        }
        "robust" => o.robust = flag()?,
        "outlier_sigma" => o.outlier_sigma = positive()?,
        "epsilon" => o.epsilon = positive()?,
        "max_iter" => {
            let x = positive()?;
            if x.fract() != 0.0 || x > 100000.0 {
                return Err("max_iter must be a whole number".into());
            }
            o.max_iter = x as usize;
        }
        "tolerance" => o.tolerance = positive()?,
        "chi2_threshold" => o.chi2_threshold = positive()?,
        "arc_gap" => o.arc_gap = positive()?,
        body => {
            let Some(&(code, _)) = MAJOR.iter().find(|(_, n)| n.eq_ignore_ascii_case(body)) else {
                return Err(format!("unknown setting {key}"));
            };
            toggle(&mut f.exclude, &[code], flag()?);
        }
    }
    Ok(())
}
