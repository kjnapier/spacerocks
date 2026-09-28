use std::collections::HashMap;

use numpy::{IntoPyArray, PyArray1, PyArray2, PyReadonlyArray1, PyReadonlyArray2};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use rayon::prelude::*;

use spacerocks::orbfit::{self, Astrometry, Engine, FitFlag, FitOptions, NongravAuto, OrbitFit, PriorFit, UpdateRoute};
use spacerocks::{SpaceRock, Time};

use crate::py_observing::observatory::PyObservatory;
use crate::py_observing::observer::PyObserver;
use crate::py_spice::spicekernel::PySpiceKernel;
use crate::py_time::time::PyTime;
use crate::rockcollection::extract_times;
use crate::PySpaceRock;

pub(crate) fn err<E: std::fmt::Display>(e: E) -> PyErr {
    PyValueError::new_err(e.to_string())
}

pub(crate) fn f64_vec(obj: &Bound<'_, PyAny>, what: &str) -> PyResult<Vec<f64>> {
    if let Ok(a) = obj.extract::<PyReadonlyArray1<f64>>() {
        return Ok(a.as_slice().map(|s| s.to_vec()).unwrap_or_else(|_| a.as_array().to_vec()));
    }
    obj.extract::<Vec<f64>>().map_err(|_| PyValueError::new_err(format!("{} must be a 1-d array of floats", what)))
}

/// A scalar or an array of length n (None -> None).
pub(crate) fn sigma_vec(obj: Option<&Bound<'_, PyAny>>, n: usize, what: &str) -> PyResult<Option<Vec<f64>>> {
    let Some(obj) = obj else { return Ok(None) };
    if obj.is_none() {
        return Ok(None);
    }
    if let Ok(s) = obj.extract::<f64>() {
        return Ok(Some(vec![s; n]));
    }
    let v = f64_vec(obj, what)?;
    if v.len() != n {
        return Err(PyValueError::new_err(format!("{} has {} entries for {} detections", what, v.len(), n)));
    }
    Ok(Some(v))
}

/// Observer positions, velocities and (on request) accelerations, barycentric J2000 (AU, AU/day,
/// AU/day²), from an Observatory, an MPC code, a sequence of MPC codes, a sequence of Observers,
/// or an (n, 3), (n, 6) or (n, 9) array of positions[, velocities[, accelerations]].
type ObserverStates = (Vec<[f64; 3]>, Vec<[f64; 3]>, Vec<[f64; 3]>);

fn observer_states(observer: &Bound<'_, PyAny>, times: &[Time], kernel: &spacerocks::SpiceKernel, accelerations: bool) -> PyResult<ObserverStates> {
    let n = times.len();
    let py = observer.py();
    let jd: Vec<f64> = times.iter().map(|t| t.tdb().jd()).collect();
    let split = |s: Vec<[f64; 6]>| -> (Vec<[f64; 3]>, Vec<[f64; 3]>) { s.iter().map(|s| ([s[0], s[1], s[2]], [s[3], s[4], s[5]])).unzip() };

    // A single site is handled as a sequence of the same site.
    let site = if let Ok(o) = observer.cast::<PyObservatory>() {
        Some(o.borrow().inner.clone())
    } else if let Ok(code) = observer.extract::<String>() {
        Some(spacerocks::Observatory::from_obscode(&code).map_err(err)?)
    } else {
        None
    };
    if let Some(site) = site {
        let states = |jd: &[f64]| -> Result<Vec<[f64; 6]>, String> {
            jd.par_iter()
                .map(|&t| {
                    let t = Time::new(t, "tdb", "jd").map_err(|e| e.to_string())?;
                    let o = site.at(&t, "J2000", "SSB", kernel).map_err(|e| e.to_string())?;
                    let v = o.velocity.unwrap_or_else(nalgebra::Vector3::zeros);
                    Ok([o.position.x, o.position.y, o.position.z, v.x, v.y, v.z])
                })
                .collect()
        };
        return py
            .detach(|| {
                let (pos, vel) = split(states(&jd)?);
                let acc = if accelerations {
                    let dt = 2.0 / 86400.0;
                    let p = states(&jd.iter().map(|t| t + dt).collect::<Vec<_>>())?;
                    let m = states(&jd.iter().map(|t| t - dt).collect::<Vec<_>>())?;
                    p.iter().zip(&m).map(|(p, m)| [(p[3] - m[3]) / (2.0 * dt), (p[4] - m[4]) / (2.0 * dt), (p[5] - m[5]) / (2.0 * dt)]).collect()
                } else {
                    Vec::new()
                };
                Ok((pos, vel, acc))
            })
            .map_err(|e: String| err(e));
    }
    if let Ok(a) = observer.extract::<PyReadonlyArray2<f64>>() {
        let a = a.as_array();
        let w = a.shape()[1];
        if a.shape()[0] != n || !(w == 3 || w == 6 || w == 9) {
            return Err(PyValueError::new_err(format!("observer array has shape {:?}, expected ({}, 3), ({}, 6) or ({}, 9)", a.shape(), n, n, n)));
        }
        let col = |k: usize| -> Vec<[f64; 3]> { if w > 3 * k { (0..n).map(|i| [a[[i, 3 * k]], a[[i, 3 * k + 1]], a[[i, 3 * k + 2]]]).collect() } else { Vec::new() } };
        return Ok((col(0), col(1), col(2)));
    }
    if let Ok(obs) = observer.extract::<Vec<PyObserver>>() {
        if obs.len() != n {
            return Err(PyValueError::new_err(format!("{} observers for {} detections", obs.len(), n)));
        }
        let mut out = (Vec::new(), Vec::new(), Vec::new());
        for o in obs {
            let mut o = o.inner;
            if o.origin.name().to_uppercase() != "SSB" {
                return Err(PyValueError::new_err("observers must be barycentric (origin 'SSB')"));
            }
            o.change_reference_plane("J2000").map_err(err)?;
            out.0.push([o.position.x, o.position.y, o.position.z]);
            if let Some(v) = o.velocity {
                out.1.push([v.x, v.y, v.z]);
            }
        }
        if out.1.len() != n {
            out.1.clear();
        }
        return Ok(out);
    }
    if let Ok(codes) = observer.extract::<Vec<String>>() {
        if codes.len() != n {
            return Err(PyValueError::new_err(format!("{} observatory codes for {} detections", codes.len(), n)));
        }
        return py
            .detach(|| {
                let (pos, vel) = split(orbfit::observer_states(&codes, &jd, kernel)?);
                let acc = if accelerations { orbfit::observer_accelerations(&codes, &jd, kernel)? } else { Vec::new() };
                Ok((pos, vel, acc))
            })
            .map_err(|e: String| err(e));
    }
    Err(PyValueError::new_err(
        "observer must be an Observatory, an MPC code, a sequence of MPC codes or Observers, or an (n, 3), (n, 6) or (n, 9) array of barycentric J2000 positions (AU), velocities (AU/day) and accelerations (AU/day^2)",
    ))
}

/// Earth-fixed positions (ITRF93, AU) of ground stations given as an Observatory, an MPC code or
/// a sequence of MPC codes; NaN for stations that are not on the ground. `None` for other kinds
/// of `observer`.
fn earth_fixed_sites(observer: &Bound<'_, PyAny>, n: usize) -> PyResult<Option<Vec<[f64; 3]>>> {
    let site = |o: &spacerocks::Observatory| o.earth_fixed_position().unwrap_or([f64::NAN; 3]);
    if let Ok(o) = observer.cast::<PyObservatory>() {
        return Ok(Some(vec![site(&o.borrow().inner); n]));
    }
    let codes: Vec<String> = if let Ok(c) = observer.extract::<String>() {
        vec![c; n]
    } else if observer.extract::<PyReadonlyArray2<f64>>().is_ok() {
        return Ok(None);
    } else if let Ok(c) = observer.extract::<Vec<String>>() {
        c
    } else {
        return Ok(None);
    };
    let mut cache: HashMap<String, [f64; 3]> = HashMap::new();
    codes
        .iter()
        .map(|c| {
            if c.trim().is_empty() {
                return Ok([f64::NAN; 3]);
            }
            if let Some(p) = cache.get(c) {
                return Ok(*p);
            }
            let p = site(&spacerocks::Observatory::from_obscode(c).map_err(err)?);
            cache.insert(c.clone(), p);
            Ok(p)
        })
        .collect::<PyResult<Vec<_>>>()
        .map(Some)
}

/// Rate and radar measurements, as passed from Python.
#[derive(Default)]
pub struct Extras<'a, 'py> {
    pub ra_rate: Option<&'a Bound<'py, PyAny>>,
    pub dec_rate: Option<&'a Bound<'py, PyAny>>,
    pub sigma_ra_rate: Option<&'a Bound<'py, PyAny>>,
    pub sigma_dec_rate: Option<&'a Bound<'py, PyAny>>,
    pub delay: Option<&'a Bound<'py, PyAny>>,
    pub doppler: Option<&'a Bound<'py, PyAny>>,
    pub sigma_delay: Option<&'a Bound<'py, PyAny>>,
    pub sigma_doppler: Option<&'a Bound<'py, PyAny>>,
    pub frequency: Option<&'a Bound<'py, PyAny>>,
    pub transmitter: Option<&'a Bound<'py, PyAny>>,
    /// Attach observer velocities whenever the observer gives them (for predicted rates).
    pub want_velocity: bool,
}

const C_AU_DAY: f64 = spacerocks::constants::SPEED_OF_LIGHT;

/// A scalar, an array of length n, or None.
fn column(obj: Option<&Bound<'_, PyAny>>, n: usize, what: &str) -> PyResult<Option<Vec<f64>>> {
    sigma_vec(obj, n, what)
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn build_astrometry(
    ra: &Bound<'_, PyAny>,
    dec: &Bound<'_, PyAny>,
    epoch: &Bound<'_, PyAny>,
    observer: &Bound<'_, PyAny>,
    kernel: &spacerocks::SpiceKernel,
    sigma_ra: Option<&Bound<'_, PyAny>>,
    sigma_dec: Option<&Bound<'_, PyAny>>,
    timescale: &str,
    extras: &Extras,
) -> PyResult<(Astrometry, Vec<f64>)> {
    let ra = f64_vec(ra, "ra")?;
    let dec = f64_vec(dec, "dec")?;
    let times = extract_times(epoch, timescale)?;
    if ra.len() != dec.len() || ra.len() != times.len() {
        return Err(PyValueError::new_err(format!("ra, dec and epoch have lengths {}, {} and {}", ra.len(), dec.len(), times.len())));
    }
    let n = ra.len();
    let ra_rate = column(extras.ra_rate, n, "ra_rate")?;
    let dec_rate = column(extras.dec_rate, n, "dec_rate")?;
    let delay = column(extras.delay, n, "delay")?;
    let doppler = column(extras.doppler, n, "doppler")?;
    let radar = delay.is_some() || doppler.is_some();
    // Where the transmitting antenna is: explicit states, stations given by code, or (failing
    // both) the receiving station extrapolated back, which needs its acceleration.
    let tx = extras.transmitter.filter(|t| !t.is_none());
    let tx_states = tx.and_then(|t| t.extract::<PyReadonlyArray2<f64>>().ok()).is_some();
    let tx_sites = if !radar || tx_states {
        None
    } else if let Some(t) = tx {
        Some(earth_fixed_sites(t, n)?.ok_or_else(|| PyValueError::new_err("transmitter must be an MPC code, a sequence of MPC codes, or an (n, 6) array of states"))?)
    } else {
        earth_fixed_sites(observer, n)?
    };
    let need_acc = radar && !tx_states && tx_sites.as_ref().map(|s| s.iter().any(|p| !p[0].is_finite())).unwrap_or(true);

    let (pos, vel, acc) = observer_states(observer, &times, kernel, need_acc)?;
    let jd: Vec<f64> = times.iter().map(|t| t.tdb().jd()).collect();
    let mut a = Astrometry::new(jd.clone(), ra, dec, pos).map_err(err)?;
    let sr = sigma_vec(sigma_ra, n, "sigma_ra")?;
    let sd = sigma_vec(sigma_dec, n, "sigma_dec")?;
    if sr.is_some() || sd.is_some() {
        let sr = sr.unwrap_or_else(|| a.sigma_ra.clone());
        let sd = sd.unwrap_or_else(|| a.sigma_dec.clone());
        a = a.with_sigma(&sr, &sd).map_err(err)?;
    }
    if extras.want_velocity && vel.len() == n && !(ra_rate.is_some() || dec_rate.is_some() || radar) {
        a = a.with_observer_velocity(vel.clone()).map_err(err)?;
    }
    if ra_rate.is_some() || dec_rate.is_some() || radar {
        if vel.len() != n {
            return Err(PyValueError::new_err("rates and radar need observer velocities: pass an Observatory, MPC codes, Observers with velocities, or an (n, 6) array"));
        }
        a = a.with_observer_velocity(vel).map_err(err)?;
    }
    if ra_rate.is_some() || dec_rate.is_some() {
        let (Some(rr), Some(dr)) = (ra_rate, dec_rate) else {
            return Err(PyValueError::new_err("give both ra_rate and dec_rate"));
        };
        let srr = column(extras.sigma_ra_rate, n, "sigma_ra_rate")?.unwrap_or_else(|| vec![f64::NAN; n]);
        let sdr = column(extras.sigma_dec_rate, n, "sigma_dec_rate")?.unwrap_or_else(|| vec![f64::NAN; n]);
        a = a.with_rates(rr, dr, &srr, &sdr).map_err(err)?;
    }
    let mut frequency = vec![f64::NAN; n];
    if radar {
        // Python units: delay and its sigma in seconds, Doppler and its sigma in Hz at the
        // transmit frequency (Hz); the fit uses days and round-trip range rate (AU/day).
        let delay_s = delay.unwrap_or_else(|| vec![f64::NAN; n]);
        let doppler_hz = doppler.unwrap_or_else(|| vec![f64::NAN; n]);
        if doppler_hz.iter().any(|d| d.is_finite()) {
            frequency = column(extras.frequency, n, "frequency")?.ok_or_else(|| PyValueError::new_err("Doppler measurements need the transmit frequency (Hz)"))?;
        }
        let s_delay = column(extras.sigma_delay, n, "sigma_delay")?.unwrap_or_else(|| vec![f64::NAN; n]);
        let s_doppler = column(extras.sigma_doppler, n, "sigma_doppler")?.unwrap_or_else(|| vec![1.0; n]);
        let delay_d: Vec<f64> = delay_s.iter().map(|s| s / 86400.0).collect();
        let doppler_au: Vec<f64> = doppler_hz.iter().zip(&frequency).map(|(d, f)| -C_AU_DAY * d / f).collect();
        let s_delay_d: Vec<f64> = s_delay.iter().map(|s| s / 86400.0).collect();
        let s_doppler_au: Vec<f64> = s_doppler.iter().zip(&frequency).map(|(s, f)| (C_AU_DAY * s / f).abs()).collect();
        for (i, (d, f)) in doppler_hz.iter().zip(&frequency).enumerate() {
            if d.is_finite() && !(f.is_finite() && *f != 0.0) {
                return Err(PyValueError::new_err(format!("Doppler measurement {} has no transmit frequency", i)));
            }
        }
        a = a.with_radar(delay_d, doppler_au, &s_delay_d, &s_doppler_au).map_err(err)?;
        if !acc.is_empty() {
            a = a.with_observer_acceleration(acc).map_err(err)?;
        }
        if tx_states {
            let arr = tx.unwrap().extract::<PyReadonlyArray2<f64>>()?;
            let arr = arr.as_array();
            if arr.shape() != [n, 6] {
                return Err(PyValueError::new_err(format!("transmitter states have shape {:?}, expected ({}, 6)", arr.shape(), n)));
            }
            a = a.with_transmitter((0..n).map(|i| [arr[[i, 0]], arr[[i, 1]], arr[[i, 2]], arr[[i, 3]], arr[[i, 4]], arr[[i, 5]]]).collect()).map_err(err)?;
        }
        if let Some(sites) = tx_sites {
            a = a.with_transmitter_site(sites).map_err(err)?;
        }
    }
    Ok((a, frequency))
}

fn parse_nongrav(nongrav: Option<&Bound<'_, PyAny>>) -> PyResult<[bool; 3]> {
    let mut mask = [false; 3];
    let Some(obj) = nongrav else { return Ok(mask) };
    if obj.is_none() {
        return Ok(mask);
    }
    if let Ok(b) = obj.extract::<bool>() {
        mask[1] = b;
        return Ok(mask);
    }
    let names: Vec<String> = if let Ok(s) = obj.extract::<String>() {
        s.to_uppercase().split('A').filter(|p| !p.is_empty()).map(|p| format!("A{}", p.trim_matches(|c: char| c == ',' || c == ' '))).collect()
    } else {
        obj.extract::<Vec<String>>().map_err(|_| PyValueError::new_err("nongrav must be True, 'A2', 'A1A2A3' or a list like ['A1', 'A2']"))?
    };
    for name in names {
        match name.trim().to_uppercase().as_str() {
            "A1" => mask[0] = true,
            "A2" => mask[1] = true,
            "A3" => mask[2] = true,
            other => return Err(PyValueError::new_err(format!("unknown non-gravitational parameter '{}'", other))),
        }
    }
    Ok(mask)
}

#[allow(clippy::too_many_arguments)]
fn options(nongrav: Option<&Bound<'_, PyAny>>, nongrav_thresholds: Option<[f64; 3]>, gofr: Option<[f64; 5]>, max_iter: usize, epsilon: f64, chi2_threshold: f64, arc_gap: f64, iod: &str, engine: &str, robust: bool, outlier_sigma: f64, conv_frac: f64, per_arc: bool) -> PyResult<FitOptions> {
    let engine = match engine.to_lowercase().as_str() {
        "cartesian" => Engine::Cartesian,
        "bk_native" | "bk" => Engine::BkNative,
        other => return Err(PyValueError::new_err(format!("engine must be 'cartesian' or 'bk_native', not '{}'", other))),
    };
    let (bk_fallback, herget) = match iod.to_lowercase().as_str() {
        "auto" => (true, false),
        "gauss" => (false, false),
        "herget" => (false, true),
        other => return Err(PyValueError::new_err(format!("iod must be 'auto', 'gauss' or 'herget', not '{}'", other))),
    };
    if !(outlier_sigma > 0.0) {
        return Err(PyValueError::new_err("outlier_sigma must be positive"));
    }
    let auto = nongrav.and_then(|n| n.extract::<String>().ok()).map(|s| s.trim().eq_ignore_ascii_case("auto")).unwrap_or(false);
    let nongrav_auto = match (auto, nongrav_thresholds) {
        (true, Some([accept_reduced_chi2, delta_chi2_per_param, nsigma])) => Some(NongravAuto { accept_reduced_chi2, delta_chi2_per_param, nsigma }),
        (true, None) => Some(NongravAuto::default()),
        (false, Some(_)) => return Err(PyValueError::new_err("nongrav_thresholds applies only to nongrav='auto'")),
        (false, None) => None,
    };
    let fit_nongrav = if auto { [false; 3] } else { parse_nongrav(nongrav)? };
    if per_arc && (auto || !fit_nongrav.iter().any(|&b| b)) {
        return Err(PyValueError::new_err("per_arc needs an explicit choice of non-gravitational parameters, e.g. nongrav='A1A2A3'"));
    }
    if !(conv_frac >= 0.0) {
        return Err(PyValueError::new_err("conv_frac must be non-negative"));
    }
    Ok(FitOptions { fit_nongrav, nongrav_auto, gofr, max_iter, epsilon, chi2_threshold, arc_gap, bk_fallback, herget, engine, robust, outlier_sigma, conv_frac, nongrav_per_arc: per_arc, ..FitOptions::default() })
}

/// Barycentric J2000 (epoch, state, nongrav) of a SpaceRock.
fn rock_state(rock: &SpaceRock, kernel: &spacerocks::SpiceKernel) -> PyResult<(f64, [f64; 6], [f64; 3])> {
    let mut r = rock.clone();
    r.to_ssb(kernel).map_err(err)?;
    r.change_reference_plane("J2000").map_err(err)?;
    let s = [r.position.x, r.position.y, r.position.z, r.velocity.x, r.velocity.y, r.velocity.z];
    Ok((r.epoch.tdb().jd(), s, r.nongrav().unwrap_or([0.0; 3])))
}

/// The result of an orbit fit.
#[pyclass(name = "OrbitFit", module = "spacerocks.orbfit")]
pub struct PyOrbitFit {
    pub inner: OrbitFit,
    pub name: String,
    /// Transmit frequency per detection (Hz), to express Doppler residuals in Hz.
    pub frequency: Vec<f64>,
    /// Keys of the detections the fit covers (to update it later).
    pub keys: Vec<u64>,
    /// How `fit(..., prior=...)` handled it, or None.
    pub route: Option<UpdateRoute>,
}

/// Per-detection residuals in Python units: RA, Dec (rad), rates (rad/day), delay (s), Doppler
/// (Hz at the transmit frequency).
fn residual_array<'py>(py: Python<'py>, per: &[[f64; 6]], frequency: &[f64]) -> Bound<'py, PyArray2<f64>> {
    numpy::ndarray::Array2::from_shape_fn((per.len(), 6), |(i, j)| match j {
        4 => per[i][4] * 86400.0,
        5 => -per[i][5] * frequency.get(i).copied().unwrap_or(f64::NAN) / C_AU_DAY,
        _ => per[i][j],
    })
    .into_pyarray(py)
}

#[pymethods]
impl PyOrbitFit {
    /// The fitted orbit as a barycentric J2000 SpaceRock (with its non-gravitational parameters
    /// when they were fitted).
    #[getter]
    fn rock(&self) -> PyResult<PySpaceRock> {
        let f = &self.inner;
        let t = Time::new(f.epoch, "tdb", "jd").map_err(err)?;
        let s = f.state;
        let mut rock = SpaceRock::from_xyz(&self.name, s[0], s[1], s[2], s[3], s[4], s[5], t, "J2000", "SSB").map_err(err)?;
        if f.fit_nongrav.iter().any(|&b| b) || f.nongrav.iter().any(|&a| a != 0.0) {
            rock.set_nongrav(f.nongrav[0], f.nongrav[1], f.nongrav[2]);
        }
        Ok(PySpaceRock { inner: rock })
    }

    /// With `per_arc=True`, arc B's (A1, A2, A3): the detections after the epoch (AU/day^2);
    /// NaN otherwise. `nongrav` then holds arc A's.
    #[getter]
    fn nongrav_arc2<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        let f = &self.inner;
        (if f.per_arc { f.nongrav_arc2.to_vec() } else { vec![f64::NAN; 3] }).into_pyarray(py)
    }

    /// 1-sigma uncertainties of `nongrav_arc2` (NaN where not fitted).
    #[getter]
    fn nongrav_arc2_sigma<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.inner.nongrav_arc2_sigma().to_vec().into_pyarray(py)
    }

    /// Which detections (in input order) the final fit used: all of them unless `robust=True`
    /// rejected some; empty if there is no orbit.
    #[getter]
    fn used<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<bool>> {
        self.inner.used.clone().into_pyarray(py)
    }

    /// Barycentric J2000 state [x, y, z, vx, vy, vz] (AU, AU/day).
    #[getter]
    fn state<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.inner.state.to_vec().into_pyarray(py)
    }

    /// Epoch of the state (a TDB Time).
    #[getter]
    fn epoch(&self) -> PyResult<PyTime> {
        Ok(PyTime { inner: Time::new(self.inner.epoch, "tdb", "jd").map_err(err)? })
    }

    /// Covariance of the state and fitted non-gravitational parameters, (npar, npar).
    #[getter]
    fn covariance<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let n = self.inner.npar;
        let a = numpy::ndarray::Array2::from_shape_vec((n, n), self.inner.covariance.clone()).map_err(err)?;
        Ok(a.into_pyarray(py))
    }

    /// 6x6 covariance of the state.
    #[getter]
    fn state_covariance<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        let c = self.inner.state_covariance();
        numpy::ndarray::Array2::from_shape_fn((6, 6), |(i, j)| c[i][j]).into_pyarray(py)
    }

    /// Non-gravitational parameters [A1, A2, A3] (AU/day^2).
    #[getter]
    fn nongrav<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.inner.nongrav.to_vec().into_pyarray(py)
    }

    /// 1-sigma uncertainties of the fitted non-gravitational parameters (NaN if not fitted).
    #[getter]
    fn nongrav_sigma<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.inner.nongrav_sigma().to_vec().into_pyarray(py)
    }

    #[getter]
    fn chi2(&self) -> f64 {
        self.inner.chi2
    }

    #[getter]
    fn ndof(&self) -> i64 {
        self.inner.ndof
    }

    #[getter]
    fn niter(&self) -> usize {
        self.inner.niter
    }

    /// layup's flag: 0 converged, 1 did not converge, 2 reduced chi-square too large, 3 no
    /// initial orbit converged, 4 incremental fit failed, 5 no initial orbit, 6 non-grav
    /// parameter unconstrained, 7 prior covariance not positive definite and 8 update too large
    /// (sequential updates without a refit), 9 implausible hyperbolic excess speed.
    #[getter]
    fn flag(&self) -> i32 {
        self.inner.flag.code()
    }

    #[getter]
    fn status(&self) -> &'static str {
        self.inner.flag.description()
    }

    #[getter]
    fn converged(&self) -> bool {
        self.inner.converged()
    }

    /// Residuals (observed - computed) at the last iteration, shape (n, 6): RA (times cos Dec)
    /// and Dec in radians, RA and Dec rates in radians/day, radar delay in seconds and Doppler in
    /// Hz; NaN where a detection has no such measurement.
    #[getter]
    fn residuals<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        residual_array(py, &self.inner.residuals, &self.frequency)
    }

    /// With `fit(..., prior=...)`: how the prior was brought up to date: "skip" (same
    /// detections), "sequential" (detections added; only they were fitted, with the prior's
    /// covariance as a prior), "sequential_fallback" (that moved too far, so all were refitted),
    /// "full" (detections removed or changed; all refitted from the prior), or "cold" (no usable
    /// prior). None otherwise.
    #[getter]
    fn route(&self) -> Option<&'static str> {
        self.route.map(|r| r.name())
    }

    /// 64-bit keys of the detections this fit covers (a hash of each detection's inputs).
    #[getter]
    fn keys<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<u64>> {
        self.keys.clone().into_pyarray(py)
    }

    /// Order-independent fingerprint of the detections this fit covers (16 hex digits).
    #[getter]
    fn fingerprint(&self) -> String {
        format!("{:016x}", orbfit::fingerprint(&self.keys))
    }

    fn __repr__(&self) -> String {
        let f = &self.inner;
        let route = self.route.map(|r| format!(", route={}", r.name())).unwrap_or_default();
        format!("OrbitFit(flag={} [{}], chi2={:.6}, ndof={}, niter={}, epoch={:.8} TDB{})", f.flag.code(), f.flag.description(), f.chi2, f.ndof, f.niter, f.epoch, route)
    }
}

/// Fit an orbit to astrometry, sky-motion rates and radar (a port of layup's orbit fit).
///
/// ra, dec: radians (J2000/ICRF); NaN for detections without optical astrometry (radar).
/// epoch: Times, or Julian dates in `timescale` (radar: the receive time). observer: an
/// Observatory, an MPC code, a sequence of MPC codes or Observers, or an (n, 3), (n, 6) or
/// (n, 9) array of barycentric J2000 positions (AU), velocities (AU/day) and accelerations
/// (AU/day^2). sigma_ra, sigma_dec: on-sky uncertainties in radians (scalar or per detection;
/// default 1"). ra_rate, dec_rate: on-sky sky-motion rates (cos(Dec) dRA/dt and dDec/dt,
/// radians/day, as ADES raRate/decRate), NaN where not measured, with sigma_ra_rate,
/// sigma_dec_rate (default 24"/day). delay: round-trip radar delay (s); doppler: Doppler shift
/// (Hz) at the transmit `frequency` (Hz); with sigma_delay (s, default 1 us) and sigma_doppler
/// (Hz, default 1 Hz); NaN where not measured. transmitter: for bistatic radar, the MPC code(s)
/// of the transmitting antenna, or an (n, 6) array of its barycentric J2000 state at the
/// transmit time. Default: monostatic. For stations given by MPC code or Observatory the
/// transmitter's state is computed exactly at the model's transmit time; for observers given
/// as arrays or Observers it is extrapolated from the observer's velocity and acceleration. initial: a SpaceRock to start from instead of
/// Gauss's method. nongrav: non-gravitational parameters to fit ('A2', 'A1A2A3', ['A1', 'A2'],
/// or True for A2), or 'auto' to choose the model as layup does: kept gravity-only if its
/// reduced chi-square is at most 1.5, else the first of A2/A1/A3, then A1A2, then A1A2A3 that
/// lowers chi-square by more than 9 per parameter with every parameter above 3 sigma.
/// nongrav_thresholds: those three numbers, for 'auto'. gofr: Marsden g(r) as [alpha, m, n, k, r0] (default inverse square).
/// iod: 'auto' (Gauss's method, then a Bernstein-Khushalani seed if no Gauss candidate
/// converges), 'gauss' (Gauss only) or 'herget' (Herget's method on the primary arc, no
/// fallback), as in layup. engine: 'cartesian' (default) or 'bk_native', layup's engines for the
/// gravity-only fits of the pipeline; 'bk_native' fits in Bernstein-Khushalani parameters with a
/// bound-orbit prior on the line-of-sight velocity, which steadies short arcs of distant objects
/// (fits with non-gravitational parameters stay Cartesian; a fit from `initial` uses the engine). robust: reject outliers (detections more
/// than `outlier_sigma` from the orbit, re-evaluated every round), and if layup's pipeline finds
/// no orbit, start over from the short window of detections that gives one and widen it (not in
/// layup; off by default). `OrbitFit.used` marks the detections in the final fit. conv_frac:
/// layup's scaled convergence test; a positive value accepts a step once every parameter moves by
/// less than that fraction of its formal sigma (or 1e-12, whichever is larger). per_arc: fit the
/// non-gravitational parameters separately for the detections before the fit epoch (arc A,
/// `nongrav`) and after it (arc B, `nongrav_arc2`), with one state: layup's comet linkage.
/// Give an `initial` orbit whose epoch lies between the two apparitions.
///
/// prior: an OrbitFit from an earlier call, to bring up to date with the current detections as
/// layup's incremental fitting does (`OrbitFit.route` says how): unchanged detections return
/// the prior; added detections are fitted alone with the prior's covariance as a Gaussian prior
/// (a sequential update), falling back to a refit of all of them if the orbit moves more than
/// `max_update_sigma` prior sigmas; otherwise all are refitted from the prior.
#[pyfunction]
#[pyo3(signature = (ra, dec, epoch, observer, kernel, sigma_ra=None, sigma_dec=None, timescale="utc", initial=None, prior=None, max_update_sigma=4.0, nongrav=None, nongrav_thresholds=None, gofr=None, max_iter=100, epsilon=1e-9, chi2_threshold=10.0, arc_gap=90.0, iod="auto", engine="cartesian", robust=false, outlier_sigma=4.0, conv_frac=0.0, per_arc=false, name="rock", ra_rate=None, dec_rate=None, sigma_ra_rate=None, sigma_dec_rate=None, delay=None, doppler=None, sigma_delay=None, sigma_doppler=None, frequency=None, transmitter=None))]
#[allow(clippy::too_many_arguments)]
pub fn fit<'py>(
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
    prior: Option<PyRef<PyOrbitFit>>,
    max_update_sigma: f64,
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
    per_arc: bool,
    name: &str,
    ra_rate: Option<&Bound<'py, PyAny>>,
    dec_rate: Option<&Bound<'py, PyAny>>,
    sigma_ra_rate: Option<&Bound<'py, PyAny>>,
    sigma_dec_rate: Option<&Bound<'py, PyAny>>,
    delay: Option<&Bound<'py, PyAny>>,
    doppler: Option<&Bound<'py, PyAny>>,
    sigma_delay: Option<&Bound<'py, PyAny>>,
    sigma_doppler: Option<&Bound<'py, PyAny>>,
    frequency: Option<&Bound<'py, PyAny>>,
    transmitter: Option<&Bound<'py, PyAny>>,
) -> PyResult<PyOrbitFit> {
    let k = &kernel.inner;
    let extras = Extras { ra_rate, dec_rate, sigma_ra_rate, sigma_dec_rate, delay, doppler, sigma_delay, sigma_doppler, frequency, transmitter, want_velocity: false };
    let (a, frequency) = build_astrometry(ra, dec, epoch, observer, k, sigma_ra, sigma_dec, timescale, &extras)?;
    let mut opts = options(nongrav, nongrav_thresholds, gofr, max_iter, epsilon, chi2_threshold, arc_gap, iod, engine, robust, outlier_sigma, conv_frac, per_arc)?;
    opts.max_update_sigma = max_update_sigma;
    let keys = a.detection_keys();
    if let Some(p) = prior {
        if initial.is_some() {
            return Err(PyValueError::new_err("give initial or prior, not both"));
        }
        let pf = PriorFit { fit: p.inner.clone(), keys: p.keys.clone() };
        let (fit, route) = py.detach(|| orbfit::update_orbit(&a, Some(&pf), k, &opts));
        return Ok(PyOrbitFit { inner: fit, name: name.to_string(), frequency, keys, route: Some(route) });
    }
    let init = match initial {
        Some(r) => Some(rock_state(&r.inner, k)?),
        None => None,
    };
    let fit = py.detach(|| orbfit::determine_orbit(&a, init, k, &opts));
    Ok(PyOrbitFit { inner: fit, name: name.to_string(), frequency, keys, route: None })
}

/// Update `prior` (an OrbitFit) with new detections only: layup's `sequential_update` without
/// its refit. The detections here are fitted alone, from the prior's state at its epoch, with the
/// prior's covariance as a Gaussian prior on the parameters it fitted. The result reports the
/// chi-square of these detections only (ndof = their residual rows) and the posterior
/// covariance. flag 8 means the orbit moved more than `max_update_sigma` prior sigmas, so the
/// update can't be trusted: refit everything instead, or use `fit(..., prior=...)`, which does
/// that for you. Arguments as for `fit`.
#[pyfunction]
#[pyo3(signature = (prior, ra, dec, epoch, observer, kernel, sigma_ra=None, sigma_dec=None, timescale="utc", max_update_sigma=4.0, max_iter=100, epsilon=1e-9, chi2_threshold=10.0, name="rock", ra_rate=None, dec_rate=None, sigma_ra_rate=None, sigma_dec_rate=None, delay=None, doppler=None, sigma_delay=None, sigma_doppler=None, frequency=None, transmitter=None))]
#[allow(clippy::too_many_arguments)]
pub fn sequential_update<'py>(
    py: Python<'py>,
    prior: PyRef<PyOrbitFit>,
    ra: &Bound<'py, PyAny>,
    dec: &Bound<'py, PyAny>,
    epoch: &Bound<'py, PyAny>,
    observer: &Bound<'py, PyAny>,
    kernel: PyRef<PySpiceKernel>,
    sigma_ra: Option<&Bound<'py, PyAny>>,
    sigma_dec: Option<&Bound<'py, PyAny>>,
    timescale: &str,
    max_update_sigma: f64,
    max_iter: usize,
    epsilon: f64,
    chi2_threshold: f64,
    name: &str,
    ra_rate: Option<&Bound<'py, PyAny>>,
    dec_rate: Option<&Bound<'py, PyAny>>,
    sigma_ra_rate: Option<&Bound<'py, PyAny>>,
    sigma_dec_rate: Option<&Bound<'py, PyAny>>,
    delay: Option<&Bound<'py, PyAny>>,
    doppler: Option<&Bound<'py, PyAny>>,
    sigma_delay: Option<&Bound<'py, PyAny>>,
    sigma_doppler: Option<&Bound<'py, PyAny>>,
    frequency: Option<&Bound<'py, PyAny>>,
    transmitter: Option<&Bound<'py, PyAny>>,
) -> PyResult<PyOrbitFit> {
    let k = &kernel.inner;
    let extras = Extras { ra_rate, dec_rate, sigma_ra_rate, sigma_dec_rate, delay, doppler, sigma_delay, sigma_doppler, frequency, transmitter, want_velocity: false };
    let (a, frequency) = build_astrometry(ra, dec, epoch, observer, k, sigma_ra, sigma_dec, timescale, &extras)?;
    if !prior.inner.converged() {
        return Err(PyValueError::new_err("the prior fit did not converge"));
    }
    let opts = FitOptions { max_iter, epsilon, chi2_threshold, max_update_sigma, ..FitOptions::default() };
    let prior_fit = prior.inner.clone();
    let (fit, _) = py.detach(|| orbfit::sequential_update(&a, &prior_fit, None, k, &opts));
    let mut keys = prior.keys.clone();
    keys.extend(a.detection_keys());
    Ok(PyOrbitFit { inner: fit, name: name.to_string(), frequency, keys, route: None })
}

/// Fit many objects at once, in parallel. `ids` labels each detection with its object; the
/// other arguments are as for `fit`. Returns a dict of arrays, one row per object in order of
/// first appearance: id, flag, chi2, ndof, niter, epoch (TDB JD), state (m, 6), covariance
/// (m, 6, 6), nongrav (m, 3), nongrav_sigma (m, 3); and `used`, one bool per detection in input
/// order, true where the detection was in its object's final fit. Also, to bring the results up
/// to date later: fit_nongrav (m, 3) bool, parameter_covariance (a list of (npar, npar) arrays),
/// fingerprint (16 hex digits per object), keys (a list of uint64 arrays: the detections each
/// object's fit covers), and route.
///
/// prior: the dict from an earlier `fit_many`, to update those fits with the current detections
/// as layup's `incremental_orbitfit` does, per object (see `fit`'s `prior`): "skip" if its
/// detections are unchanged, "sequential" if detections were only added (fitting just those,
/// with a refit if the orbit moves more than `max_update_sigma` prior sigmas:
/// "sequential_fallback"), "full" if some were removed or changed, "cold" for objects without a
/// converged prior. `route` says which.
#[pyfunction]
#[pyo3(signature = (ids, ra, dec, epoch, observer, kernel, sigma_ra=None, sigma_dec=None, timescale="utc", prior=None, max_update_sigma=4.0, nongrav=None, nongrav_thresholds=None, gofr=None, max_iter=100, epsilon=1e-9, chi2_threshold=10.0, arc_gap=90.0, iod="auto", engine="cartesian", robust=false, outlier_sigma=4.0, conv_frac=0.0, per_arc=false, ra_rate=None, dec_rate=None, sigma_ra_rate=None, sigma_dec_rate=None, delay=None, doppler=None, sigma_delay=None, sigma_doppler=None, frequency=None, transmitter=None))]
#[allow(clippy::too_many_arguments)]
pub fn fit_many<'py>(
    py: Python<'py>,
    ids: &Bound<'py, PyAny>,
    ra: &Bound<'py, PyAny>,
    dec: &Bound<'py, PyAny>,
    epoch: &Bound<'py, PyAny>,
    observer: &Bound<'py, PyAny>,
    kernel: PyRef<PySpiceKernel>,
    sigma_ra: Option<&Bound<'py, PyAny>>,
    sigma_dec: Option<&Bound<'py, PyAny>>,
    timescale: &str,
    prior: Option<&Bound<'py, PyDict>>,
    max_update_sigma: f64,
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
    per_arc: bool,
    ra_rate: Option<&Bound<'py, PyAny>>,
    dec_rate: Option<&Bound<'py, PyAny>>,
    sigma_ra_rate: Option<&Bound<'py, PyAny>>,
    sigma_dec_rate: Option<&Bound<'py, PyAny>>,
    delay: Option<&Bound<'py, PyAny>>,
    doppler: Option<&Bound<'py, PyAny>>,
    sigma_delay: Option<&Bound<'py, PyAny>>,
    sigma_doppler: Option<&Bound<'py, PyAny>>,
    frequency: Option<&Bound<'py, PyAny>>,
    transmitter: Option<&Bound<'py, PyAny>>,
) -> PyResult<Bound<'py, PyDict>> {
    let k = &kernel.inner;
    let extras = Extras { ra_rate, dec_rate, sigma_ra_rate, sigma_dec_rate, delay, doppler, sigma_delay, sigma_doppler, frequency, transmitter, want_velocity: false };
    let labels: Vec<String> = ids
        .try_iter()?
        .map(|x| x.and_then(|x| x.str().map(|s| s.to_string())))
        .collect::<PyResult<_>>()?;
    let (all, _) = build_astrometry(ra, dec, epoch, observer, k, sigma_ra, sigma_dec, timescale, &extras)?;
    if labels.len() != all.len() {
        return Err(PyValueError::new_err(format!("{} ids for {} detections", labels.len(), all.len())));
    }
    let mut opts = options(nongrav, nongrav_thresholds, gofr, max_iter, epsilon, chi2_threshold, arc_gap, iod, engine, robust, outlier_sigma, conv_frac, per_arc)?;
    opts.max_update_sigma = max_update_sigma;
    let priors = match prior {
        Some(d) => Some(parse_priors(d)?),
        None => None,
    };

    let mut order: Vec<String> = Vec::new();
    let mut groups: HashMap<&str, Vec<usize>> = HashMap::new();
    for (i, l) in labels.iter().enumerate() {
        groups.entry(l.as_str()).or_insert_with(|| {
            order.push(l.clone());
            Vec::new()
        }).push(i);
    }
    let subsets: Vec<Astrometry> = order.iter().map(|l| all.subset(&groups[l.as_str()])).collect();
    let results: Vec<(OrbitFit, Option<UpdateRoute>)> = py.detach(|| {
        let serial = FitOptions { parallel: false, ..opts.clone() };
        subsets
            .par_iter()
            .zip(order.par_iter())
            .map(|(a, id)| match &priors {
                Some(p) => {
                    let (f, r) = orbfit::update_orbit(a, p.get(id.as_str()), k, &serial);
                    (f, Some(r))
                }
                None => (orbfit::determine_orbit(a, None, k, &serial), None),
            })
            .collect()
    });
    let routes: Vec<Option<UpdateRoute>> = results.iter().map(|r| r.1).collect();
    let fits: Vec<OrbitFit> = results.into_iter().map(|r| r.0).collect();

    let m = fits.len();
    let d = PyDict::new(py);
    d.set_item("id", order.clone())?;
    d.set_item("flag", fits.iter().map(|f| f.flag.code()).collect::<Vec<_>>().into_pyarray(py))?;
    d.set_item("chi2", fits.iter().map(|f| f.chi2).collect::<Vec<_>>().into_pyarray(py))?;
    d.set_item("ndof", fits.iter().map(|f| f.ndof).collect::<Vec<_>>().into_pyarray(py))?;
    d.set_item("niter", fits.iter().map(|f| f.niter as i64).collect::<Vec<_>>().into_pyarray(py))?;
    d.set_item("epoch", fits.iter().map(|f| f.epoch).collect::<Vec<_>>().into_pyarray(py))?;
    let state = numpy::ndarray::Array2::from_shape_fn((m, 6), |(i, j)| fits[i].state[j]);
    d.set_item("state", state.into_pyarray(py))?;
    let covs: Vec<[[f64; 6]; 6]> = fits.iter().map(|f| f.state_covariance()).collect();
    let cov = numpy::ndarray::Array3::from_shape_fn((m, 6, 6), |(i, j, l)| covs[i][j][l]);
    d.set_item("covariance", cov.into_pyarray(py))?;
    let ng = numpy::ndarray::Array2::from_shape_fn((m, 3), |(i, j)| fits[i].nongrav[j]);
    d.set_item("nongrav", ng.into_pyarray(py))?;
    let ngs: Vec<[f64; 3]> = fits.iter().map(|f| f.nongrav_sigma()).collect();
    let ngs = numpy::ndarray::Array2::from_shape_fn((m, 3), |(i, j)| ngs[i][j]);
    d.set_item("nongrav_sigma", ngs.into_pyarray(py))?;
    if opts.nongrav_per_arc {
        let ng2 = numpy::ndarray::Array2::from_shape_fn((m, 3), |(i, j)| if fits[i].per_arc { fits[i].nongrav_arc2[j] } else { f64::NAN });
        d.set_item("nongrav_arc2", ng2.into_pyarray(py))?;
        let s2: Vec<[f64; 3]> = fits.iter().map(|f| f.nongrav_arc2_sigma()).collect();
        let s2 = numpy::ndarray::Array2::from_shape_fn((m, 3), |(i, j)| s2[i][j]);
        d.set_item("nongrav_arc2_sigma", s2.into_pyarray(py))?;
    }
    // Per detection, in input order: whether it was in its object's final fit.
    let mut used = vec![false; all.len()];
    for (l, f) in order.iter().zip(&fits) {
        for (k, &i) in groups[l.as_str()].iter().enumerate() {
            used[i] = f.used.get(k).copied().unwrap_or(false);
        }
    }
    d.set_item("used", used.into_pyarray(py))?;
    let fitted = numpy::ndarray::Array2::from_shape_fn((m, 3), |(i, j)| fits[i].fit_nongrav[j]);
    d.set_item("fit_nongrav", fitted.into_pyarray(py))?;
    let pcov: Vec<Bound<'py, PyArray2<f64>>> = fits
        .iter()
        .map(|f| numpy::ndarray::Array2::from_shape_fn((f.npar, f.npar), |(i, j)| f.covariance.get(i * f.npar + j).copied().unwrap_or(f64::NAN)).into_pyarray(py))
        .collect();
    d.set_item("parameter_covariance", pcov)?;
    let keys: Vec<Vec<u64>> = subsets.iter().map(|a| a.detection_keys()).collect();
    d.set_item("fingerprint", keys.iter().map(|k| format!("{:016x}", orbfit::fingerprint(k))).collect::<Vec<_>>())?;
    d.set_item("keys", keys.into_iter().map(|k| k.into_pyarray(py)).collect::<Vec<_>>())?;
    d.set_item("route", routes.iter().map(|r| r.map(|r| r.name())).collect::<Vec<_>>())?;
    Ok(d)
}

/// Priors by object id, from an earlier `fit_many` result.
fn parse_priors(d: &Bound<'_, PyDict>) -> PyResult<HashMap<String, PriorFit>> {
    let get = |k: &str| d.get_item(k)?.ok_or_else(|| PyValueError::new_err(format!("prior is missing '{}' (pass the dict fit_many returned)", k)));
    let ids: Vec<String> = get("id")?.try_iter()?.map(|x| x.and_then(|x| x.str().map(|s| s.to_string()))).collect::<PyResult<_>>()?;
    let flag: Vec<i32> = get("flag")?.extract()?;
    let epoch: Vec<f64> = get("epoch")?.extract()?;
    let state: Vec<[f64; 6]> = get("state")?.extract::<numpy::PyReadonlyArray2<f64>>()?.as_array().rows().into_iter().map(|r| [r[0], r[1], r[2], r[3], r[4], r[5]]).collect();
    let nongrav: Vec<[f64; 3]> = get("nongrav")?.extract::<numpy::PyReadonlyArray2<f64>>()?.as_array().rows().into_iter().map(|r| [r[0], r[1], r[2]]).collect();
    let fitted: Vec<[bool; 3]> = get("fit_nongrav")?.extract::<numpy::PyReadonlyArray2<bool>>()?.as_array().rows().into_iter().map(|r| [r[0], r[1], r[2]]).collect();
    let pcov: Vec<numpy::PyReadonlyArray2<f64>> = get("parameter_covariance")?.extract()?;
    let keys: Vec<Vec<u64>> = get("keys")?.try_iter()?.map(|x| x.and_then(|x| x.extract::<Vec<u64>>())).collect::<PyResult<_>>()?;
    let chi2: Vec<f64> = get("chi2")?.extract()?;
    let ndof: Vec<i64> = get("ndof")?.extract()?;
    let niter: Vec<i64> = get("niter")?.extract()?;
    let arc2: Option<Vec<[f64; 3]>> = match d.get_item("nongrav_arc2")? {
        Some(x) => Some(x.extract::<numpy::PyReadonlyArray2<f64>>()?.as_array().rows().into_iter().map(|r| [r[0], r[1], r[2]]).collect()),
        None => None,
    };
    let m = ids.len();
    if [flag.len(), epoch.len(), state.len(), nongrav.len(), fitted.len(), pcov.len(), keys.len(), chi2.len(), ndof.len(), niter.len()].iter().any(|&l| l != m) {
        return Err(PyValueError::new_err("prior: every entry needs one row per object"));
    }
    let mut out = HashMap::new();
    for i in 0..m {
        let c = pcov[i].as_array();
        let npar = c.nrows();
        let per_arc = arc2.as_ref().map(|a| a[i].iter().all(|v| v.is_finite())).unwrap_or(false);
        let fit = OrbitFit {
            epoch: epoch[i],
            state: state[i],
            nongrav: nongrav[i],
            fit_nongrav: fitted[i],
            nongrav_arc2: arc2.as_ref().map(|a| a[i]).filter(|_| per_arc).unwrap_or([0.0; 3]),
            per_arc,
            covariance: c.iter().copied().collect(),
            npar,
            chi2: chi2[i],
            ndof: ndof[i],
            niter: niter[i].max(0) as usize,
            flag: FitFlag::from_code(flag[i]).unwrap_or(FitFlag::NotConverged),
            residuals: Vec::new(),
            used: Vec::new(),
        };
        out.insert(ids[i].clone(), PriorFit { fit, keys: keys[i].clone() });
    }
    Ok(out)
}

/// Residuals (observed - computed) of `rock`, shape (n, 6): RA (times cos Dec) and Dec in
/// radians, RA and Dec rates in radians/day, radar delay in seconds and Doppler in Hz; NaN
/// where not measured. Light-time corrected, with the ASSIST force model (and the rock's
/// non-gravitational parameters). Arguments as for `fit`.
#[pyfunction]
#[pyo3(signature = (rock, ra, dec, epoch, observer, kernel, timescale="utc", epsilon=1e-9, ra_rate=None, dec_rate=None, sigma_ra_rate=None, sigma_dec_rate=None, delay=None, doppler=None, sigma_delay=None, sigma_doppler=None, frequency=None, transmitter=None))]
#[allow(clippy::too_many_arguments)]
pub fn residuals<'py>(
    py: Python<'py>,
    rock: PyRef<PySpaceRock>,
    ra: &Bound<'py, PyAny>,
    dec: &Bound<'py, PyAny>,
    epoch: &Bound<'py, PyAny>,
    observer: &Bound<'py, PyAny>,
    kernel: PyRef<PySpiceKernel>,
    timescale: &str,
    epsilon: f64,
    ra_rate: Option<&Bound<'py, PyAny>>,
    dec_rate: Option<&Bound<'py, PyAny>>,
    sigma_ra_rate: Option<&Bound<'py, PyAny>>,
    sigma_dec_rate: Option<&Bound<'py, PyAny>>,
    delay: Option<&Bound<'py, PyAny>>,
    doppler: Option<&Bound<'py, PyAny>>,
    sigma_delay: Option<&Bound<'py, PyAny>>,
    sigma_doppler: Option<&Bound<'py, PyAny>>,
    frequency: Option<&Bound<'py, PyAny>>,
    transmitter: Option<&Bound<'py, PyAny>>,
) -> PyResult<Bound<'py, PyArray2<f64>>> {
    let k = &kernel.inner;
    let extras = Extras { ra_rate, dec_rate, sigma_ra_rate, sigma_dec_rate, delay, doppler, sigma_delay, sigma_doppler, frequency, transmitter, want_velocity: false };
    let (a, frequency) = build_astrometry(ra, dec, epoch, observer, k, None, None, timescale, &extras)?;
    let (e, s, ng) = rock_state(&rock.inner, k)?;
    let opts = FitOptions { epsilon, ..FitOptions::default() };
    let r = py.detach(|| orbfit::residuals(&a, e, &s, &ng, k, &opts, false)).map_err(err)?;
    Ok(residual_array(py, &r.per_detection(a.len()), &frequency))
}

/// Bernstein–Khushalani linear initial orbit (layup's `run_bk_iod`): a closed-form
/// straight-line fit in a tangent-plane frame, suited to short arcs of distant objects where
/// Gauss's method is ill-conditioned. Returns a barycentric J2000 SpaceRock at `epoch0` (a TDB
/// Julian date; default: the middle detection in time), or None. Arguments as for `fit`.
#[pyfunction]
#[pyo3(signature = (ra, dec, epoch, observer, kernel, sigma_ra=None, sigma_dec=None, timescale="utc", epoch0=None, name="rock"))]
#[allow(clippy::too_many_arguments)]
pub fn bk_iod(
    ra: &Bound<'_, PyAny>,
    dec: &Bound<'_, PyAny>,
    epoch: &Bound<'_, PyAny>,
    observer: &Bound<'_, PyAny>,
    kernel: PyRef<PySpiceKernel>,
    sigma_ra: Option<&Bound<'_, PyAny>>,
    sigma_dec: Option<&Bound<'_, PyAny>>,
    timescale: &str,
    epoch0: Option<f64>,
    name: &str,
) -> PyResult<Option<PySpaceRock>> {
    let (a, _) = build_astrometry(ra, dec, epoch, observer, &kernel.inner, sigma_ra, sigma_dec, timescale, &Extras::default())?;
    let a = a.subset(&a.optical());
    if a.is_empty() {
        return Ok(None);
    }
    let a = a.subset(&a.time_order());
    let e = epoch0.unwrap_or(a.epoch[a.len() / 2]);
    let Some(s) = orbfit::bk_iod(&a, e) else { return Ok(None) };
    let t = Time::new(e, "tdb", "jd").map_err(err)?;
    let rock = SpaceRock::from_xyz(name, s[0], s[1], s[2], s[3], s[4], s[5], t, "J2000", "SSB").map_err(err)?;
    Ok(Some(PySpaceRock { inner: rock }))
}

/// Herget's initial orbit (layup's `herget_iod`): the ranges to the first and last detections
/// (in time) are iterated, from 2, then 5, then 40 AU, until the orbit through them fits the
/// detections in between. Every detection passed is used, so give one arc (layup uses the
/// primary one). Returns a barycentric J2000 SpaceRock at the first detection (TDB), or None.
/// Arguments as for `fit`.
#[pyfunction]
#[pyo3(signature = (ra, dec, epoch, observer, kernel, sigma_ra=None, sigma_dec=None, timescale="utc", name="rock"))]
#[allow(clippy::too_many_arguments)]
pub fn herget_iod(
    py: Python<'_>,
    ra: &Bound<'_, PyAny>,
    dec: &Bound<'_, PyAny>,
    epoch: &Bound<'_, PyAny>,
    observer: &Bound<'_, PyAny>,
    kernel: PyRef<PySpiceKernel>,
    sigma_ra: Option<&Bound<'_, PyAny>>,
    sigma_dec: Option<&Bound<'_, PyAny>>,
    timescale: &str,
    name: &str,
) -> PyResult<Option<PySpaceRock>> {
    let k = &kernel.inner;
    let (a, _) = build_astrometry(ra, dec, epoch, observer, k, sigma_ra, sigma_dec, timescale, &Extras::default())?;
    let a = a.subset(&a.optical());
    let a = a.subset(&a.time_order());
    let idx: Vec<usize> = (0..a.len()).collect();
    let opts = FitOptions::default();
    let Some((e, s)) = py.detach(|| orbfit::herget_iod(&a, &idx, k, &opts)) else { return Ok(None) };
    let t = Time::new(e, "tdb", "jd").map_err(err)?;
    let rock = SpaceRock::from_xyz(name, s[0], s[1], s[2], s[3], s[4], s[5], t, "J2000", "SSB").map_err(err)?;
    Ok(Some(PySpaceRock { inner: rock }))
}

/// Predicted positions of a fitted orbit with their uncertainty (layup's `predict`), at `epoch`
/// (Times, or Julian dates in `timescale`) for `observer` (as for `fit`). The orbit keeps its
/// fitted non-gravitational parameters. Returns a dict: `epoch` (TDB JD), `ra`, `dec` (radians,
/// astrometric: light-time corrected), `delta` (AU), `covariance` (n, 2, 2) on the sky in
/// radians² along (RA·cos Dec, Dec), and the 1-sigma error ellipse: `sigma_major`,
/// `sigma_minor` (arcsec) and `pa` (degrees, North through East). The covariance is the fit's
/// mapped linearly: `B C B^T`, with B the partials of the direction including light time, and C
/// the fit's full covariance (so fitted A1-A3 contribute; layup maps the state covariance only).
#[pyfunction]
#[pyo3(signature = (fit, epoch, observer, kernel, timescale="utc", epsilon=1e-9))]
pub fn predict<'py>(py: Python<'py>, fit: PyRef<PyOrbitFit>, epoch: &Bound<'py, PyAny>, observer: &Bound<'py, PyAny>, kernel: PyRef<PySpiceKernel>, timescale: &str, epsilon: f64) -> PyResult<Bound<'py, PyDict>> {
    let k = &kernel.inner;
    let n = extract_times(epoch, timescale)?.len();
    let zeros = vec![0.0f64; n].into_pyarray(py).into_any();
    let (a, _) = build_astrometry(&zeros, &zeros, epoch, observer, k, None, None, timescale, &Extras::default())?;
    let f = fit.inner.clone();
    let opts = FitOptions { epsilon, ..FitOptions::default() };
    let p = py.detach(|| orbfit::predict_astrometry(&f, &a, k, &opts)).map_err(err)?;
    let d = PyDict::new(py);
    let as_rad = 180.0 / std::f64::consts::PI * 3600.0;
    d.set_item("epoch", p.iter().map(|x| x.epoch).collect::<Vec<_>>().into_pyarray(py))?;
    d.set_item("ra", p.iter().map(|x| x.ra).collect::<Vec<_>>().into_pyarray(py))?;
    d.set_item("dec", p.iter().map(|x| x.dec).collect::<Vec<_>>().into_pyarray(py))?;
    d.set_item("delta", p.iter().map(|x| x.delta).collect::<Vec<_>>().into_pyarray(py))?;
    let cov = numpy::ndarray::Array3::from_shape_fn((p.len(), 2, 2), |(i, u, v)| p[i].cov[u][v]);
    d.set_item("covariance", cov.into_pyarray(py))?;
    let ell: Vec<(f64, f64, f64)> = p.iter().map(|x| x.ellipse()).collect();
    d.set_item("sigma_major", ell.iter().map(|e| e.0 * as_rad).collect::<Vec<_>>().into_pyarray(py))?;
    d.set_item("sigma_minor", ell.iter().map(|e| e.1 * as_rad).collect::<Vec<_>>().into_pyarray(py))?;
    d.set_item("pa", ell.iter().map(|e| e.2).collect::<Vec<_>>().into_pyarray(py))?;
    Ok(d)
}

/// Original and future orbits of a long-period comet (layup's `comet`): the barycentric
/// osculating orbit (with the mass of the Sun and planets) at `reference_distance` AU from the
/// barycenter on the way in (original) and on the way out (future). `orbit` is an OrbitFit or a
/// SpaceRock (with its non-gravitational parameters; `gofr` as for `fit`, e.g. the water-ice
/// law). Returns {"original": ..., "future": ...}, each None if the orbit never gets that far,
/// else a dict: epoch (TDB JD), distance (AU), reached (False if the planetary ephemeris ended
/// first, so the elements are where it stopped: exact for two-body motion out there), inv_a
/// (1/AU; x 1e6 for the CODE catalogue's units), a, e, q (AU), inc (degrees, ecliptic).
#[pyfunction]
#[pyo3(signature = (orbit, kernel, reference_distance=250.0, gofr=None, epsilon=1e-9))]
pub fn comet_orbits<'py>(py: Python<'py>, orbit: &Bound<'py, PyAny>, kernel: PyRef<PySpiceKernel>, reference_distance: f64, gofr: Option<[f64; 5]>, epsilon: f64) -> PyResult<Bound<'py, PyDict>> {
    let k = &kernel.inner;
    let (epoch, state, nongrav) = if let Ok(f) = orbit.extract::<PyRef<PyOrbitFit>>() {
        (f.inner.epoch, f.inner.state, f.inner.nongrav)
    } else if let Ok(r) = orbit.extract::<PyRef<PySpaceRock>>() {
        rock_state(&r.inner, k)?
    } else {
        return Err(PyValueError::new_err("orbit must be an OrbitFit or a SpaceRock"));
    };
    let opts = FitOptions { epsilon, gofr, ..FitOptions::default() };
    let (orig, fut) = py.detach(|| {
        (
            orbfit::comet_orbit(epoch, &state, &nongrav, false, reference_distance, k, &opts),
            orbfit::comet_orbit(epoch, &state, &nongrav, true, reference_distance, k, &opts),
        )
    });
    let to_dict = |o: Option<orbfit::CometOrbit>| -> PyResult<Option<Bound<'py, PyDict>>> {
        let Some(o) = o else { return Ok(None) };
        let d = PyDict::new(py);
        d.set_item("epoch", o.epoch)?;
        d.set_item("distance", o.distance)?;
        d.set_item("reached", o.reached)?;
        d.set_item("inv_a", o.inv_a)?;
        d.set_item("a", o.a)?;
        d.set_item("e", o.e)?;
        d.set_item("q", o.q)?;
        d.set_item("inc", o.inc.to_degrees())?;
        Ok(Some(d))
    };
    let d = PyDict::new(py);
    d.set_item("original", to_dict(orig)?)?;
    d.set_item("future", to_dict(fut)?)?;
    Ok(d)
}
