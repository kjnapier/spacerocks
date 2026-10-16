use pyo3::exceptions::PyOSError;
use pyo3::prelude::*;

use numpy::{IntoPyArray, PyArray2, PyArrayMethods};
use pyo3::types::PyDict;
use spacerocks::orbfit::{FitDashboard, McmcOptions, ASTEROID_PERTURBERS};

use crate::py_orbfit::fit::{build_astrometry, err, options, rock_state, Extras, PyOrbitFit};
use crate::py_spice::spicekernel::PySpiceKernel;
use crate::PySpaceRock;

/// A live orbit-fitting dashboard: a page on 127.0.0.1 to fit an orbit interactively (see
/// `spacerocks.orbfit.dashboard`). It keeps running until `close()` or until it is garbage
/// collected.
#[pyclass(name = "FitDashboard", module = "spacerocks.orbfit")]
pub struct PyFitDashboard {
    inner: Option<FitDashboard>,
    url: String,
}

impl PyFitDashboard {
    fn dashboard(&self) -> PyResult<&FitDashboard> {
        self.inner.as_ref().ok_or_else(|| PyOSError::new_err("the dashboard is closed"))
    }

    /// A fit of the dashboard's current object (which Load and Random change). The dashboard
    /// fits optical detections only, so there are no transmit frequencies.
    fn wrap(d: &FitDashboard, fit: spacerocks::orbfit::OrbitFit) -> PyOrbitFit {
        let keys = d.astrometry().detection_keys();
        PyOrbitFit { inner: fit, name: d.name(), frequency: vec![f64::NAN; keys.len()], keys, route: None }
    }
}

#[pymethods]
impl PyFitDashboard {
    /// The page's address.
    #[getter]
    fn url(&self) -> &str {
        &self.url
    }

    /// The last fit (None before the first one finishes). Its residuals and `used` cover every
    /// detection; those left out on the page are marked unused.
    #[getter]
    fn fit(&self) -> PyResult<Option<PyOrbitFit>> {
        let d = self.dashboard()?;
        Ok(d.fit().map(|f| Self::wrap(d, f)))
    }

    /// The object's name (Load and Random on the page change it).
    #[getter]
    fn name(&self) -> PyResult<String> {
        Ok(self.dashboard()?.name())
    }

    /// Whether another object's observations are being fetched.
    #[getter]
    fn loading(&self) -> PyResult<bool> {
        Ok(self.dashboard()?.loading())
    }

    /// Fetch the observations of `designation` (a number, name or designation; None for a random
    /// numbered asteroid) from the Minor Planet Center and fit them, as the page's Load and
    /// Random buttons do. Returns at once; if the fetch fails, the page says why and keeps the
    /// current object.
    #[pyo3(signature = (designation=None))]
    fn load(&self, designation: Option<&str>) -> PyResult<()> {
        self.dashboard()?.load(designation);
        Ok(())
    }

    /// Whether a fit is running or about to.
    #[getter]
    fn fitting(&self) -> PyResult<bool> {
        Ok(self.dashboard()?.fitting())
    }

    /// Which detections are left out on the page.
    #[getter]
    fn excluded(&self) -> PyResult<Vec<bool>> {
        Ok(self.dashboard()?.excluded())
    }

    /// Start a fit with the page's settings, as its Fit button does.
    fn refit(&self) -> PyResult<()> {
        self.dashboard()?.refit();
        Ok(())
    }

    /// Block until Done is pressed on the page (and no fit is running), then return the last
    /// fit. Ctrl-C stops waiting with KeyboardInterrupt; the dashboard keeps running.
    fn wait(&self, py: Python<'_>) -> PyResult<Option<PyOrbitFit>> {
        let d = self.dashboard()?;
        let caught: std::sync::Mutex<Option<PyErr>> = std::sync::Mutex::new(None);
        let fit = py.detach(|| {
            d.wait(Some(&|| {
                Python::attach(|py| match py.check_signals() {
                    Ok(()) => false,
                    Err(e) => {
                        *caught.lock().unwrap() = Some(e);
                        true
                    }
                })
            }))
        });
        if let Some(e) = caught.into_inner().unwrap() {
            return Err(e);
        }
        Ok(fit.map(|f| Self::wrap(d, f)))
    }

    /// Sample the last fit's posterior with a Markov chain, as the page's Sample button does:
    /// `steps` random-walk Metropolis steps starting `spread` sigmas from the fit, the first
    /// `burn_in` of them left out of the posterior summary (`seed` makes it repeatable).
    /// Returns at once; `.chain` has the chain when `.sampling` turns False.
    #[pyo3(signature = (steps=5000, burn_in=1000, spread=3.0, seed=None))]
    fn sample(&self, steps: usize, burn_in: usize, spread: f64, seed: Option<u64>) -> PyResult<()> {
        self.dashboard()?.sample(Some(McmcOptions { steps, burn_in, spread, seed }));
        Ok(())
    }

    /// Whether a chain is running or about to.
    #[getter]
    fn sampling(&self) -> PyResult<bool> {
        Ok(self.dashboard()?.sampling())
    }

    /// The burn-in: the steps the posterior summary drops from the start of the chain.
    #[getter]
    fn burn_in(&self) -> PyResult<usize> {
        Ok(self.dashboard()?.mcmc_options().burn_in)
    }

    #[setter]
    fn set_burn_in(&self, burn_in: usize) -> PyResult<()> {
        self.dashboard()?.set_burn_in(burn_in);
        Ok(())
    }

    /// The last chain (None before one runs, and after a new fit): a dict with `params` (steps x
    /// parameters: the barycentric J2000 state at `epoch`, AU and AU/day, then the fitted
    /// non-gravitational parameters, AU/day²), `names`, `chi2`, `accepted`, `epoch` (TDB JD),
    /// `burn_in` and `posterior` (the params after the burn-in).
    #[getter]
    fn chain<'py>(&self, py: Python<'py>) -> PyResult<Option<Bound<'py, PyDict>>> {
        let d = self.dashboard()?;
        let Some(c) = d.chain() else { return Ok(None) };
        let burn_in = d.mcmc_options().burn_in.min(c.steps.len());
        let npar = c.npar();
        let rows = |from: usize| -> PyResult<Bound<'py, PyArray2<f64>>> {
            let flat: Vec<f64> = c.steps[from..].iter().flat_map(|s| s.params.iter().copied()).collect();
            Ok(flat.into_pyarray(py).reshape([c.steps.len() - from, npar])?)
        };
        let mut names: Vec<String> = ["x", "y", "z", "vx", "vy", "vz"].iter().map(|s| s.to_string()).collect();
        names.extend((0..3).filter(|&k| c.fit_nongrav[k]).map(|k| format!("A{}", k + 1)));
        let out = PyDict::new(py);
        out.set_item("params", rows(0)?)?;
        out.set_item("posterior", rows(burn_in)?)?;
        out.set_item("names", names)?;
        out.set_item("chi2", c.steps.iter().map(|s| s.chi2).collect::<Vec<_>>().into_pyarray(py))?;
        out.set_item("accepted", c.steps.iter().map(|s| s.accepted).collect::<Vec<_>>().into_pyarray(py))?;
        out.set_item("epoch", c.epoch)?;
        out.set_item("burn_in", burn_in)?;
        Ok(Some(out))
    }

    /// Stop the dashboard: its page goes away and any running fit stops.
    fn close(&mut self) {
        self.inner = None;
    }

    fn __repr__(&self) -> String {
        format!("FitDashboard({}{})", self.url, if self.inner.is_none() { ", closed" } else { "" })
    }
}

/// Fit an orbit interactively: serve a live dashboard on 127.0.0.1:`port` (0 picks a free port),
/// open it in the browser, and start a first fit at once. Returns a FitDashboard at once (the
/// fits run in the background); `.wait()` blocks until you press Done on the page and returns
/// the last fit, and `.fit` has it at any time.
///
/// On the page you can switch forces (the 16 massive asteroids, each planet, the Moon and Pluto,
/// Earth J2-J4, solar J2, the relativistic model, the non-gravitational parameters A1, A2, A3 or
/// their automatic choice, and the comet g(r)) and the solver (engine, initial orbit method,
/// robust outlier rejection and its threshold, IAS15 tolerance, iteration limit, convergence
/// tolerance, chi-square threshold and arc gap), refit from scratch or from the last orbit, stop
/// a fit, and leave detections out by clicking them on the residual plot. It shows the
/// chi-square of each iteration live, the residuals, the heliocentric ecliptic elements with
/// their uncertainties, a map of the orbit, and a log of the fits.
///
/// The arguments are those of `fit` (optical astrometry only) for the starting settings, plus
/// asteroids: include the 16 massive asteroids (needs sb441-n16.bsp in the kernel);
/// open_browser: open the page.
#[pyfunction]
#[pyo3(signature = (ra, dec, epoch, observer, kernel, sigma_ra=None, sigma_dec=None, timescale="utc", initial=None, nongrav=None, nongrav_thresholds=None, gofr=None, max_iter=100, epsilon=1e-9, chi2_threshold=10.0, arc_gap=90.0, iod="auto", engine="cartesian", robust=false, outlier_sigma=4.0, conv_frac=0.0, asteroids=true, name="rock", port=0, open_browser=true))]
#[allow(clippy::too_many_arguments)]
pub fn dashboard<'py>(
    py: Python<'py>,
    ra: &Bound<'py, PyAny>,
    dec: &Bound<'py, PyAny>,
    epoch: &Bound<'py, PyAny>,
    observer: &Bound<'py, PyAny>,
    kernel: PyRef<PySpiceKernel>,
    sigma_ra: Option<&Bound<'py, PyAny>>,
    sigma_dec: Option<&Bound<'py, PyAny>>,
    timescale: &str,
    initial: Option<PyRef<PySpaceRock>>,
    nongrav: Option<&Bound<'py, PyAny>>,
    nongrav_thresholds: Option<[f64; 3]>,
    gofr: Option<[f64; 5]>,
    max_iter: usize,
    epsilon: f64,
    chi2_threshold: f64,
    arc_gap: f64,
    iod: &str,
    engine: &str,
    robust: bool,
    outlier_sigma: f64,
    conv_frac: f64,
    asteroids: bool,
    name: &str,
    port: u16,
    open_browser: bool,
) -> PyResult<PyFitDashboard> {
    let k = &kernel.inner;
    let extras = Extras { ra_rate: None, dec_rate: None, sigma_ra_rate: None, sigma_dec_rate: None, delay: None, doppler: None, sigma_delay: None, sigma_doppler: None, frequency: None, transmitter: None, want_velocity: false };
    let (a, _) = build_astrometry(ra, dec, epoch, observer, k, sigma_ra, sigma_dec, timescale, &extras)?;
    let mut opts = options(nongrav, nongrav_thresholds, gofr, max_iter, epsilon, chi2_threshold, arc_gap, iod, engine, robust, outlier_sigma, conv_frac, false)?;
    if !asteroids {
        opts.forces.exclude = ASTEROID_PERTURBERS.to_vec();
    }
    let init = match initial {
        Some(r) => Some(rock_state(&r.inner, k)?),
        None => None,
    };
    let d = FitDashboard::start(port, name, a, k.clone(), opts, init).map_err(err)?;
    let url = d.url().to_string();
    if open_browser {
        py.import("webbrowser")?.call_method1("open", (url.as_str(),))?;
    }
    Ok(PyFitDashboard { inner: Some(d), url })
}

/// The `spacerocks-fit` command: fit an orbit to an object's Minor Planet Center observations
/// in the dashboard. `args` are the command-line arguments (default `sys.argv[1:]`); run
/// `spacerocks-fit --help` for them. With no designation it picks a random numbered asteroid.
/// Prints the fit when Done is pressed on the page.
#[pyfunction]
#[pyo3(signature = (args=None))]
pub fn main(py: Python<'_>, args: Option<Vec<String>>) -> PyResult<()> {
    let args = match args {
        Some(a) => a,
        None => py.import("sys")?.getattr("argv")?.extract::<Vec<String>>()?.into_iter().skip(1).collect(),
    };
    let caught: std::sync::Mutex<Option<PyErr>> = std::sync::Mutex::new(None);
    let open = |url: &str| Python::attach(|py| py.import("webbrowser")?.call_method1("open", (url,)).map(|_| ())).map_err(|e| e.to_string());
    let result = py.detach(|| {
        spacerocks::orbfit::cli::run(
            &args,
            &open,
            Some(&|| {
                Python::attach(|py| match py.check_signals() {
                    Ok(()) => false,
                    Err(e) => {
                        *caught.lock().unwrap() = Some(e);
                        true
                    }
                })
            }),
        )
    });
    if let Some(e) = caught.into_inner().unwrap() {
        // Ctrl-C quits quietly, as a command should.
        if e.is_instance_of::<pyo3::exceptions::PyKeyboardInterrupt>(py) {
            return Ok(());
        }
        return Err(e);
    }
    result.map_err(|e| {
        eprintln!("spacerocks-fit: {e}");
        pyo3::exceptions::PySystemExit::new_err(1)
    })
}
