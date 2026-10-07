//! Python bindings for `spacerocks::checker`: identify detections with known objects.

use std::path::PathBuf;

use numpy::{IntoPyArray, PyArray1, PyArray2, PyReadonlyArray1};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};

use spacerocks::batch::BatchOptions;
use spacerocks::checker::{self, mpc, Catalog, CheckOptions};

use crate::py_orbfit::fit::{build_astrometry, err, sigma_vec, Extras, PyOrbitFit};
use crate::py_spice::spicekernel::PySpiceKernel;
use crate::rockcollection::extract_times;

const ARCSEC: f64 = std::f64::consts::PI / (180.0 * 3600.0);

/// A catalog of orbits to check detections against (``spacerocks.checker.Catalog``).
///
/// Orbits are barycentric J2000 states at TDB epochs, with a covariance when one is known
/// (fits, MPC ``mpc_orb`` records) or else the MPC uncertainty parameter U (MPCORB).
#[pyclass(from_py_object, name = "Catalog", module = "spacerocks.checker")]
#[derive(Clone)]
pub struct PyCatalog {
    pub inner: Catalog,
}

fn json_text(obj: &Bound<'_, PyAny>) -> PyResult<String> {
    if let Ok(s) = obj.extract::<String>() {
        let p = PathBuf::from(&s);
        if !s.trim_start().starts_with(['{', '[']) && p.exists() {
            return std::fs::read_to_string(&p).map_err(err);
        }
        return Ok(s);
    }
    if let Ok(p) = obj.extract::<PathBuf>() {
        return std::fs::read_to_string(&p).map_err(err);
    }
    let json = obj.py().import("json")?;
    json.call_method1("dumps", (obj,))?.extract()
}

fn epoch_jd(epoch: &Bound<'_, PyAny>, timescale: &str) -> PyResult<f64> {
    if let Ok(jd) = epoch.extract::<f64>() {
        return Ok(spacerocks::Time::new(jd, timescale, "jd").map_err(err)?.tdb().jd());
    }
    let t = extract_times(epoch, timescale)?;
    if t.len() != 1 {
        return Err(PyValueError::new_err("give one epoch"));
    }
    Ok(t[0].tdb().jd())
}

#[pymethods]
impl PyCatalog {
    /// The MPC's catalog of minor-planet orbits (MPCORB, ~1.4 million orbits).
    ///
    /// ``path``: ``mpcorb_extended.json``, ``MPCORB.DAT``, or either gzipped. By default
    /// ``~/.spacerocks/mpc/mpcorb_extended.json.gz`` (or ``$SPACEROCKS_MPC_DIR``), downloaded from
    /// the MPC first if it is missing and ``download`` is set (``update=True`` downloads a fresh
    /// copy). The parsed orbits are cached next to the file (``<file>.srcat``), so later loads take
    /// seconds.
    #[staticmethod]
    #[pyo3(signature = (kernel, path=None, download=false, update=false))]
    fn mpcorb(py: Python<'_>, kernel: PyRef<PySpiceKernel>, path: Option<PathBuf>, download: bool, update: bool) -> PyResult<PyCatalog> {
        let k = &kernel.inner;
        let inner = py
            .detach(|| {
                let p = match path {
                    Some(p) => p,
                    None => mpc::mpcorb_path(download, update)?,
                };
                Catalog::mpcorb(Some(&p), false, k)
            })
            .map_err(err)?;
        Ok(PyCatalog { inner })
    }

    /// Orbits from fits (``OrbitFit``s, with their covariances): a dict ``{name: fit}`` or a
    /// list of fits with ``names``. ``h``: absolute magnitudes, for predicted magnitudes.
    #[staticmethod]
    #[pyo3(signature = (fits, names=None, h=None))]
    fn from_fits(fits: &Bound<'_, PyAny>, names: Option<Vec<String>>, h: Option<Vec<f64>>) -> PyResult<PyCatalog> {
        let (names, fits): (Vec<String>, Vec<spacerocks::orbfit::OrbitFit>) = if let Ok(d) = fits.cast::<PyDict>() {
            let mut n = Vec::new();
            let mut f = Vec::new();
            for (k, v) in d.iter() {
                n.push(k.extract::<String>()?);
                f.push(v.extract::<PyRef<PyOrbitFit>>()?.inner.clone());
            }
            (n, f)
        } else {
            let list: Vec<PyRef<PyOrbitFit>> = fits.extract()?;
            let names = match names {
                Some(n) => n,
                None => list.iter().map(|f| f.name.clone()).collect(),
            };
            (names, list.iter().map(|f| f.inner.clone()).collect())
        };
        let inner = Catalog::from_fits(&names, &fits, h.as_deref()).map_err(err)?;
        Ok(PyCatalog { inner })
    }

    /// Orbits from MPC ``mpc_orb`` records (the ``get-orb`` API's format, with the covariance of
    /// the Cartesian state): a dict, a list of them, JSON text, or a path to a JSON file.
    #[staticmethod]
    fn from_mpc_orb(orbs: &Bound<'_, PyAny>, kernel: PyRef<PySpiceKernel>) -> PyResult<PyCatalog> {
        let v: serde_json::Value = serde_json::from_str(&json_text(orbs)?).map_err(err)?;
        // A list of records, or one record (possibly wrapped in lists, as the API returns it).
        let records: Vec<serde_json::Value> = match &v {
            serde_json::Value::Array(a) if a.iter().any(|x| x.get("CAR").is_some() || x.get("mpc_orb").is_some()) => a.clone(),
            _ => vec![v.clone()],
        };
        let orbs: Vec<mpc::MpcOrb> = records.iter().map(mpc::parse_mpc_orb).collect::<Result<_, _>>().map_err(err)?;
        let inner = Catalog::from_mpc_orbs(&orbs, &kernel.inner).map_err(err)?;
        Ok(PyCatalog { inner })
    }

    /// Read a catalog written by ``save``.
    #[staticmethod]
    fn load(py: Python<'_>, path: PathBuf) -> PyResult<PyCatalog> {
        let inner = py.detach(|| Catalog::load(&path)).map_err(err)?;
        Ok(PyCatalog { inner })
    }

    /// Write the catalog (with its snapshot) to a binary file.
    fn save(&self, py: Python<'_>, path: PathBuf) -> PyResult<()> {
        py.detach(|| self.inner.save(&path)).map_err(err)
    }

    /// Integrate every orbit (N-body, ASSIST's force model) to ``epoch`` and keep the states as
    /// the catalog's snapshot. ``check`` moves objects by two-body motion from the closest
    /// reference state, so a snapshot within ``max_age`` days of the detections saves
    /// integrating the catalog on every call. Save the catalog to keep it.
    #[pyo3(signature = (epoch, kernel, timescale="tdb", chunk_size=64, method="nbody"))]
    fn snapshot(&mut self, py: Python<'_>, epoch: &Bound<'_, PyAny>, kernel: PyRef<PySpiceKernel>, timescale: &str, chunk_size: usize, method: &str) -> PyResult<()> {
        let t = epoch_jd(epoch, timescale)?;
        let opts = BatchOptions { chunk_size, method: spacerocks::batch::Method::from_str(method).map_err(err)?, ..BatchOptions::default() };
        let k = &kernel.inner;
        let cat = &mut self.inner;
        py.detach(|| cat.make_snapshot(t, k, &opts, 100_000));
        Ok(())
    }

    /// The orbits at ``indices`` (integers or a boolean mask), as a new catalog.
    fn select(&self, indices: &Bound<'_, PyAny>) -> PyResult<PyCatalog> {
        let idx: Vec<usize> = if let Ok(mask) = indices.extract::<PyReadonlyArray1<bool>>() {
            let m = mask.as_array();
            if m.len() != self.inner.len() {
                return Err(PyValueError::new_err(format!("mask has {} entries for {} orbits", m.len(), self.inner.len())));
            }
            m.iter().enumerate().filter(|(_, &b)| b).map(|(i, _)| i).collect()
        } else {
            let v: Vec<i64> = indices.extract()?;
            let n = self.inner.len() as i64;
            v.into_iter()
                .map(|i| {
                    let j = if i < 0 { i + n } else { i };
                    if (0..n).contains(&j) { Ok(j as usize) } else { Err(PyValueError::new_err(format!("index {} out of range", i))) }
                })
                .collect::<PyResult<_>>()?
        };
        Ok(PyCatalog { inner: self.inner.select(&idx) })
    }

    /// Append another catalog's orbits.
    fn extend(&mut self, other: PyRef<PyCatalog>) {
        self.inner.extend(&other.inner);
    }

    /// Index of the orbit named ``name`` (KeyError if absent).
    fn index(&self, name: &str) -> PyResult<usize> {
        self.inner.orbits.names.iter().position(|n| n == name).ok_or_else(|| pyo3::exceptions::PyKeyError::new_err(name.to_string()))
    }

    /// The orbit at index ``i`` as a barycentric J2000 SpaceRock at its epoch.
    fn rock(&self, i: usize) -> PyResult<crate::PySpaceRock> {
        if i >= self.inner.len() {
            return Err(PyValueError::new_err("index out of range"));
        }
        Ok(crate::PySpaceRock { inner: self.inner.rock(i).map_err(err)? })
    }

    #[getter]
    fn names(&self) -> Vec<String> {
        self.inner.orbits.names.clone()
    }

    /// TDB Julian dates of the orbits.
    #[getter]
    fn epoch<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.inner.orbits.epochs.clone().into_pyarray(py)
    }

    /// Barycentric J2000 states (n, 6), AU and AU/day.
    #[getter]
    fn states<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        let s = &self.inner.orbits.states;
        numpy::ndarray::Array2::from_shape_fn((s.len(), 6), |(i, j)| s[i][j]).into_pyarray(py)
    }

    #[getter]
    fn h<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.inner.h.clone().into_pyarray(py)
    }

    /// MPC uncertainty parameter U (NaN if unknown).
    #[getter]
    fn u<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.inner.u.clone().into_pyarray(py)
    }

    /// Whether each orbit has a covariance.
    #[getter]
    fn has_covariance<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<bool>> {
        self.inner.covariance.iter().map(|c| !c.is_empty()).collect::<Vec<_>>().into_pyarray(py)
    }

    /// TDB Julian date of the snapshot, or None.
    #[getter]
    fn snapshot_epoch(&self) -> Option<f64> {
        self.inner.has_snapshot().then_some(self.inner.snapshot_epoch)
    }

    fn __len__(&self) -> usize {
        self.inner.len()
    }

    fn __repr__(&self) -> String {
        let ncov = self.inner.covariance.iter().filter(|c| !c.is_empty()).count();
        let snap = if self.inner.has_snapshot() { format!(", snapshot at JD {} TDB", self.inner.snapshot_epoch) } else { String::new() };
        format!("Catalog({} orbits, {} with covariance{})", self.inner.len(), ncov, snap)
    }
}

/// Identify detections with known objects, like the MPC's MPChecker.
///
/// Every orbit in ``catalog`` is predicted at every detection: first by two-body motion from a
/// nearby reference state (the orbit's epoch, the catalog's snapshot, or a state integrated for
/// the purpose), then, for objects near a detection, with a full N-body prediction and its
/// uncertainty (the covariance mapped through the variational equations, or for MPCORB orbits
/// the along-track uncertainty of the U parameter). A detection is consistent with an object
/// when the Mahalanobis distance of their offset, with both uncertainties, is at most ``nsigma``.
///
/// ``ra``, ``dec``: radians (ICRF). ``epoch``: Julian dates in ``timescale`` or Time objects.
/// ``observer``: as for ``orbfit.fit`` (an Observatory, MPC code(s), Observers, or an (n, 3|6)
/// array of barycentric J2000 states). ``sigma_ra`` (of RA cos Dec), ``sigma_dec``: 1-sigma
/// uncertainties in radians (default 1"); ``correlation``: their correlation.
///
/// ``nsigma``: the consistency threshold (2-d Mahalanobis distance). ``radius`` (arcsec): also
/// report objects predicted this close even when inconsistent (0: consistent ones only).
/// ``max_age`` (days): the longest two-body step of the coarse search. ``max_uncertainty``
/// (arcsec): objects whose 1-sigma uncertainty is larger are only reported within ``radius``.
/// ``floor`` (arcsec): added in quadrature to every prediction's uncertainty. ``cross_track``:
/// for orbits without a covariance, the uncertainty across the motion as a fraction of the
/// along-track one.
///
/// Returns a dict of arrays with one value per reported (detection, object) pair, sorted by
/// detection, then consistent pairs first, then by likelihood (so ``pandas.DataFrame(result)``
/// makes a table): ``detection`` (index), ``object`` (index in the catalog), ``name``, ``ra``,
/// ``dec`` (predicted, radians), ``dra``, ``ddec`` (observed minus predicted, arcsec, along RA
/// cos Dec and Dec), ``separation`` (arcsec), ``distance`` (Mahalanobis), ``consistent``,
/// ``log_likelihood`` (of the offset, per arcsec²: ranks a precise orbit that fits above a vague
/// one), ``sigma_major``, ``sigma_minor`` (arcsec) and ``pa`` (degrees) of the prediction's
/// error ellipse and its covariance ``cov_ra``, ``cov_ra_dec``, ``cov_dec`` (arcsec², along RA
/// cos Dec and Dec, including ``floor``), ``ra_rate``, ``dec_rate`` (dRA/dt and dDec/dt,
/// radians/day), ``delta``, ``r_helio`` (AU), ``mag`` (predicted V) and ``from_covariance``.
/// Objects whose prediction fails are skipped with a warning.
#[pyfunction]
#[pyo3(signature = (catalog, ra, dec, epoch, observer, kernel, sigma_ra=None, sigma_dec=None, correlation=None, timescale="utc", nsigma=3.0, radius=60.0, max_age=30.0, max_uncertainty=600.0, floor=0.3, cross_track=0.1, epsilon=1e-9))]
#[allow(clippy::too_many_arguments)]
pub fn check<'py>(
    py: Python<'py>,
    catalog: PyRef<PyCatalog>,
    ra: &Bound<'py, PyAny>,
    dec: &Bound<'py, PyAny>,
    epoch: &Bound<'py, PyAny>,
    observer: &Bound<'py, PyAny>,
    kernel: PyRef<PySpiceKernel>,
    sigma_ra: Option<&Bound<'py, PyAny>>,
    sigma_dec: Option<&Bound<'py, PyAny>>,
    correlation: Option<&Bound<'py, PyAny>>,
    timescale: &str,
    nsigma: f64,
    radius: f64,
    max_age: f64,
    max_uncertainty: f64,
    floor: f64,
    cross_track: f64,
    epsilon: f64,
) -> PyResult<Bound<'py, PyDict>> {
    let k = &kernel.inner;
    let extras = Extras { want_velocity: true, ..Extras::default() };
    let (a, _) = build_astrometry(ra, dec, epoch, observer, k, sigma_ra, sigma_dec, timescale, &extras)?;
    let n = a.epoch.len();
    let corr = match correlation {
        Some(c) if !c.is_none() => sigma_vec(Some(c), n, "correlation")?.unwrap_or_default(),
        _ => Vec::new(),
    };
    let opts = CheckOptions { nsigma, radius: radius * ARCSEC, max_age, max_uncertainty: max_uncertainty * ARCSEC, floor: floor * ARCSEC, cross_track, epsilon, ..CheckOptions::default() };
    let cat = &catalog.inner;
    let out = py.detach(|| checker::check(cat, &a, &corr, k, &opts)).map_err(err)?;
    let m = &out.matches;
    let d = PyDict::new(py);
    let col = |f: &dyn Fn(&checker::Match) -> f64| m.iter().map(f).collect::<Vec<f64>>();
    d.set_item("detection", m.iter().map(|x| x.detection as i64).collect::<Vec<_>>().into_pyarray(py))?;
    d.set_item("object", m.iter().map(|x| x.object as i64).collect::<Vec<_>>().into_pyarray(py))?;
    d.set_item("name", PyList::new(py, m.iter().map(|x| cat.orbits.names[x.object].as_str()))?)?;
    d.set_item("ra", col(&|x| x.ra).into_pyarray(py))?;
    d.set_item("dec", col(&|x| x.dec).into_pyarray(py))?;
    d.set_item("dra", col(&|x| x.offset[0] / ARCSEC).into_pyarray(py))?;
    d.set_item("ddec", col(&|x| x.offset[1] / ARCSEC).into_pyarray(py))?;
    d.set_item("separation", col(&|x| x.separation / ARCSEC).into_pyarray(py))?;
    d.set_item("distance", col(&|x| x.distance).into_pyarray(py))?;
    d.set_item("consistent", m.iter().map(|x| x.consistent).collect::<Vec<_>>().into_pyarray(py))?;
    d.set_item("log_likelihood", col(&|x| x.log_likelihood).into_pyarray(py))?;
    let ell: Vec<(f64, f64, f64)> = m.iter().map(|x| x.ellipse()).collect();
    d.set_item("sigma_major", ell.iter().map(|e| e.0 / ARCSEC).collect::<Vec<_>>().into_pyarray(py))?;
    d.set_item("sigma_minor", ell.iter().map(|e| e.1 / ARCSEC).collect::<Vec<_>>().into_pyarray(py))?;
    d.set_item("pa", ell.iter().map(|e| e.2).collect::<Vec<_>>().into_pyarray(py))?;
    let a2 = ARCSEC * ARCSEC;
    d.set_item("cov_ra", col(&|x| x.covariance[0][0] / a2).into_pyarray(py))?;
    d.set_item("cov_ra_dec", col(&|x| x.covariance[0][1] / a2).into_pyarray(py))?;
    d.set_item("cov_dec", col(&|x| x.covariance[1][1] / a2).into_pyarray(py))?;
    d.set_item("ra_rate", col(&|x| x.ra_rate).into_pyarray(py))?;
    d.set_item("dec_rate", col(&|x| x.dec_rate).into_pyarray(py))?;
    d.set_item("delta", col(&|x| x.delta).into_pyarray(py))?;
    d.set_item("r_helio", col(&|x| x.r_helio).into_pyarray(py))?;
    d.set_item("mag", col(&|x| x.mag).into_pyarray(py))?;
    d.set_item("from_covariance", m.iter().map(|x| x.from_covariance).collect::<Vec<_>>().into_pyarray(py))?;
    if !out.failed.is_empty() {
        let list: Vec<String> = out.failed.iter().take(10).map(|(i, e)| format!("{}: {}", cat.orbits.names[*i], e)).collect();
        let msg = format!("{} object(s) could not be predicted and were skipped: {}{}", out.failed.len(), list.join("; "), if out.failed.len() > 10 { "; ..." } else { "" });
        py.import("warnings")?.call_method1("warn", (msg,))?;
    }
    Ok(d)
}

pub fn make_checker_submodule(py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    let submodule = PyModule::new(m.py(), "checker")?;
    submodule.add_class::<PyCatalog>()?;
    submodule.add_function(wrap_pyfunction!(check, submodule.clone())?)?;
    m.add_submodule(&submodule)?;
    py.import("sys")?.getattr("modules")?.set_item("spacerocks.checker", submodule.clone())?;
    submodule.setattr("__name__", "spacerocks.checker")?;
    Ok(())
}
