use pyo3::prelude::*;
use pyo3::exceptions::{PyIndexError, PyValueError};
use pyo3::{Python, PyResult};
use rayon::prelude::*;

use nalgebra::Vector3;
use numpy::{PyArray1, IntoPyArray, PyArray2};
use pyo3::types::PyDict;

use spacerocks::observing::{Apparent, Observer};
use spacerocks::batch::{self, BatchOptions, Method};
use spacerocks::state::{self, Elements, State};
use spacerocks::transforms::correct_for_ltt_vectors;
use spacerocks::{Origin, Population, ReferencePlane, Time};

use crate::py_time::time::PyTime;
use crate::PySpaceRock;
use crate::py_spice::spicekernel::PySpiceKernel;
use crate::py_observing::observer::PyObserver;
use crate::py_observing::observation::PyObservation;
use crate::py_observing::observatory::PyObservatory;

/// A collection of space rocks in one reference plane, about one origin.
///
/// The rocks are stored as a dense `(n, 6)` state array plus per-rock epochs, names and
/// properties (see `spacerocks::Population`). The reference plane and origin are shared: the
/// first rock added sets them unless they were given to the constructor, rocks in another
/// reference plane are rotated into the collection's, and adding a rock about another origin is
/// an error. Indexing builds a `SpaceRock` on demand.
#[pyclass(from_py_object)]
#[derive(Clone)]
pub struct RockCollection {
    pub inner: Population,
    /// Whether the reference plane and origin were set by the constructor (otherwise the first
    /// rock added sets them).
    fixed: bool,
}

impl RockCollection {
    pub fn from_population(inner: Population) -> Self {
        RockCollection { inner, fixed: true }
    }

    fn column<'py>(&self, py: Python<'py>, k: usize) -> Bound<'py, PyArray1<f64>> {
        let v: Vec<f64> = self.inner.states.iter().map(|s| s[k]).collect();
        v.into_pyarray(py)
    }

    fn mapped<'py>(&self, py: Python<'py>, f: fn(&State, f64) -> f64) -> Bound<'py, PyArray1<f64>> {
        let pop = &self.inner;
        py.detach(|| pop.map(f)).into_pyarray(py)
    }

    fn mapped_or_err<'py>(&self, py: Python<'py>, f: fn(&State, f64) -> Result<f64, Box<dyn std::error::Error>>) -> PyResult<Bound<'py, PyArray1<f64>>> {
        let pop = &self.inner;
        let v = py
            .detach(|| {
                pop.map(|s, mu| f(s, mu).map_err(|e| e.to_string()))
                    .into_iter()
                    .collect::<Result<Vec<f64>, String>>()
            })
            .map_err(PyValueError::new_err)?;
        Ok(v.into_pyarray(py))
    }
}

#[pymethods]
impl RockCollection {
    /// Creates a new, empty `RockCollection`. With `reference_plane` and `origin` given, rocks
    /// added later are rotated into that plane and must have that origin; otherwise the first
    /// rock added sets both.
    #[new]
    #[pyo3(signature = (reference_plane = None, origin = None))]
    pub fn new(reference_plane: Option<&str>, origin: Option<&str>) -> PyResult<Self> {
        let fixed = reference_plane.is_some() || origin.is_some();
        let plane = match reference_plane {
            Some(p) => ReferencePlane::from_str(p).map_err(PyValueError::new_err)?,
            None => ReferencePlane::J2000,
        };
        let origin = match origin {
            Some(o) => Origin::from_str(o).map_err(|e| PyValueError::new_err(e.to_string()))?,
            None => Origin::SSB,
        };
        Ok(RockCollection { inner: Population::new(plane, origin), fixed })
    }

    /// Add a rock (a copy of it).
    pub fn add(&mut self, rock: PyRef<PySpaceRock>) -> PyResult<()> {
        if self.inner.is_empty() && !self.fixed {
            self.inner.reference_plane = rock.inner.reference_plane.clone();
            self.inner.origin = rock.inner.origin.clone();
        }
        self.inner.push(rock.inner.clone()).map_err(|e| PyValueError::new_err(e.to_string()))
    }

    fn __getitem__(&self, index: isize) -> PyResult<PySpaceRock> {
        let n = self.inner.len() as isize;
        let i = if index < 0 { index + n } else { index };
        if i < 0 || i >= n {
            return Err(PyIndexError::new_err("Index out of range!"));
        }
        Ok(PySpaceRock { inner: self.inner.get(i as usize).unwrap() })
    }

    /// A new RockCollection with the rocks where `indices` (a boolean mask) is True.
    pub fn filter(&self, indices: Vec<bool>) -> PyResult<Self> {
        let inner = self.inner.filter(&indices).map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(RockCollection { inner, fixed: self.fixed })
    }

    pub fn observe(&mut self, py: Python<'_>, observer: PyRef<PyObserver>) -> PyResult<Vec<PyObservation>> {
        let o = &observer.inner;
        let pop = &self.inner;
        let observations = py
            .detach(|| {
                (0..pop.len())
                    .into_par_iter()
                    .map(|i| {
                        let mut rock = pop.get(i).unwrap();
                        rock.observe(o).map_err(|e| format!("{}: {}", rock.name, e))
                    })
                    .collect::<Result<Vec<_>, String>>()
            })
            .map_err(PyValueError::new_err)?;
        Ok(observations.into_iter().map(|obs| PyObservation { inner: obs }).collect())
    }

    /// Observe every rock from `observer` and return the results as NumPy arrays.
    ///
    /// Same computation as `observe`, but returns a dict of 1-D arrays (one entry per rock)
    /// instead of a list of Observation objects, which is much faster for large collections:
    /// `ra`, `dec` (radians), `ra_rate`, `dec_rate` (radians/day), `range` (AU),
    /// `range_rate` (AU/day), `r_helio` (AU), `phase`, `elong` (radians) and `mag`
    /// (NaN for rocks without an absolute magnitude).
    pub fn observe_arrays<'py>(&self, py: Python<'py>, observer: PyRef<PyObserver>) -> PyResult<Bound<'py, PyDict>> {
        let o = &observer.inner;
        let pop = &self.inner;
        let apps = py.detach(|| pop.apparent(o).map_err(|e| e.to_string())).map_err(PyValueError::new_err)?;
        apparent_dict(py, &apps, &[apps.len()])
    }

    /// Light-time corrected RA and Dec (radians) of every rock seen by `observer`, as an
    /// `(n, 2)` array. The collection and the observer must both be in J2000.
    pub fn calc_radec<'py>(
        &self,
        py: Python<'py>,
        observer: PyRef<PyObserver>,
    ) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let o = &observer.inner;

        if o.reference_plane != ReferencePlane::J2000 {
            return Err(PyValueError::new_err("Observer frame is not J2000. Cannot calculate RA/Dec."));
        }
        if self.inner.reference_plane != ReferencePlane::J2000 {
            return Err(PyValueError::new_err("Collection frame is not J2000. Cannot calculate RA/Dec."));
        }

        let n = self.inner.len();
        let mut data = vec![0.0f64; 2 * n];
        let states = &self.inner.states;
        let obs_vel = o.velocity.unwrap_or_else(Vector3::zeros);
        py.detach(|| {
            data.par_chunks_mut(2).zip(states.par_iter()).for_each(|(row, s)| {
                let (r, v) = state::pv(s);
                let (p, _) = correct_for_ltt_vectors(&r, &v, &o.position, &obs_vel);
                let mut ra = p.y.atan2(p.x);
                if ra < 0.0 {
                    ra += 2.0 * std::f64::consts::PI;
                }
                row[0] = ra;
                row[1] = (p.z / p.norm()).asin();
            })
        });

        let arr = numpy::ndarray::Array2::from_shape_vec((n, 2), data).unwrap();
        Ok(arr.into_pyarray(py))
    }

    /// Move every rock to `epoch` along its Keplerian orbit about the collection's origin.
    pub fn analytic_propagate(&mut self, py: Python<'_>, epoch: PyRef<PyTime>) -> PyResult<()> {
        let ep = &epoch.inner;
        let pop = &mut self.inner;
        py.detach(|| pop.analytic_propagate(ep).map_err(|e| format!("Failed to propagate rock: {}", e)))
            .map_err(PyValueError::new_err)
    }

    /// Rotate every rock into `reference_plane`.
    pub fn change_reference_plane(&mut self, reference_plane: &str) -> PyResult<()> {
        let plane = ReferencePlane::from_str(reference_plane).map_err(PyValueError::new_err)?;
        self.inner.change_reference_plane(&plane).map_err(|e| PyValueError::new_err(e.to_string()))
    }

    /// Propagate every rock to `epoch`.
    ///
    /// `method="nbody"` (default) integrates with IAS15 in the field of the Sun, planets, Moon,
    /// Pluto and 16 massive asteroids from `kernel`, with ASSIST's force model (relativity,
    /// Earth and solar harmonics, and non-gravitational forces for rocks with `set_nongrav`); `method="twobody"` uses Keplerian motion
    /// (then `kernel` is not used). Rocks are integrated together in groups of up to
    /// `chunk_size` (grouped by epoch and perihelion distance), with groups run in parallel;
    /// `chunk_size=1` integrates every rock on its own, exactly like `SpaceRock.propagate`.
    /// The collection keeps its reference plane and origin (SUN or SSB).
    #[pyo3(signature = (epoch, kernel, method = "nbody", chunk_size = 64))]
    pub fn propagate(&mut self, py: Python<'_>, epoch: PyRef<PyTime>, kernel: PyRef<PySpiceKernel>, method: &str, chunk_size: usize) -> PyResult<()> {
        let opts = BatchOptions { method: Method::from_str(method).map_err(PyValueError::new_err)?, chunk_size, ..Default::default() };
        let ep = epoch.inner.clone();
        let k = &kernel.inner;
        let pop = &mut self.inner;
        py.detach(|| {
            let mut rocks = pop.to_rocks();
            batch::propagate_batch(&mut rocks, &ep, k, &opts).map_err(|e| e.to_string())?;
            pop.update_from_rocks(&rocks).map_err(|e| e.to_string())
        })
        .map_err(PyValueError::new_err)
    }

    /// Ephemerides of every rock at many epochs.
    ///
    /// `epochs` is a list of `Time` objects or an array of Julian dates in `timescale`.
    /// `observer` is an `Observatory` (its position is computed at each epoch) or a list of
    /// `Observer` objects, one per epoch (then `epochs` may be None).
    ///
    /// Each rock is integrated once through all epochs (`method="nbody"`, the same force model
    /// as `propagate`) or moved on its Keplerian orbit (`method="twobody"`); the collection
    /// itself is not modified. Returns a dict of arrays of shape `(len(self), n_epochs)`:
    /// `ra`, `dec` (J2000 equator, radians), `ra_rate`, `dec_rate` (radians/day), `range`
    /// (AU), `range_rate` (AU/day), `r_helio` (AU), `phase`, `elong` (radians), `mag`
    /// (NaN without an absolute magnitude), and `epoch` (TDB Julian dates, shape
    /// `(n_epochs,)`). With `return_states=True` it also holds `states`, the barycentric J2000
    /// state of each rock at each epoch, shape `(len(self), n_epochs, 6)`.
    #[pyo3(signature = (epochs, observer, kernel, method = "nbody", timescale = "utc", return_states = false, chunk_size = 64))]
    pub fn ephemeris<'py>(
        &self,
        py: Python<'py>,
        epochs: Option<&Bound<'py, PyAny>>,
        observer: &Bound<'py, PyAny>,
        kernel: PyRef<PySpiceKernel>,
        method: &str,
        timescale: &str,
        return_states: bool,
        chunk_size: usize,
    ) -> PyResult<Bound<'py, PyDict>> {
        let k = &kernel.inner;
        let observers: Vec<Observer> = if let Ok(obs) = observer.cast::<PyObservatory>() {
            let obs = obs.borrow().inner.clone();
            let epochs = epochs.ok_or_else(|| PyValueError::new_err("epochs are required with an Observatory"))?;
            let times = extract_times(epochs, timescale)?;
            py.detach(|| {
                times
                    .par_iter()
                    .map(|t| obs.at(t, "J2000", "SSB", k).map_err(|e| e.to_string()))
                    .collect::<Result<Vec<_>, String>>()
            })
            .map_err(PyValueError::new_err)?
        } else {
            let list: Vec<PyObserver> = observer
                .extract()
                .map_err(|_| PyValueError::new_err("observer must be an Observatory or a list of Observers"))?;
            if let Some(e) = epochs {
                if !e.is_none() {
                    let times = extract_times(e, timescale)?;
                    let consistent = times.len() == list.len()
                        && times.iter().zip(&list).all(|(t, o)| (t.tdb().jd() - o.inner.epoch.tdb().jd()).abs() < 1e-8);
                    if !consistent {
                        return Err(PyValueError::new_err(
                            "epochs do not match the observers' epochs (pass epochs=None with a list of Observers)",
                        ));
                    }
                }
            }
            list.into_iter().map(|o| o.inner).collect()
        };

        let opts = BatchOptions {
            method: Method::from_str(method).map_err(PyValueError::new_err)?,
            chunk_size,
            with_states: return_states,
            ..Default::default()
        };
        let pop = &self.inner;
        let eph = py
            .detach(|| batch::ephemeris(&pop.to_rocks(), &observers, k, &opts).map_err(|e| e.to_string()))
            .map_err(PyValueError::new_err)?;

        let d = apparent_dict(py, &eph.apparent, &[eph.n_rocks, eph.n_epochs])?;
        d.set_item("epoch", eph.epochs.clone().into_pyarray(py))?;
        if let Some(states) = eph.states {
            let flat: Vec<f64> = states.iter().flat_map(|s| s.iter().copied()).collect();
            let arr = numpy::ndarray::Array3::from_shape_vec((eph.n_rocks, eph.n_epochs, 6), flat)
                .map_err(|e| PyValueError::new_err(e.to_string()))?;
            d.set_item("states", arr.into_pyarray(py))?;
        }
        Ok(d)
    }

    /// The collection's reference plane (shared by every rock).
    #[getter]
    pub fn reference_plane(&self) -> String {
        self.inner.reference_plane.to_string()
    }

    /// The collection's origin (shared by every rock).
    #[getter]
    pub fn origin(&self) -> String {
        self.inner.origin.to_string()
    }

    fn __len__(&self) -> usize {
        self.inner.len()
    }

    pub fn len(&self) -> usize {
        self.inner.len()
    }

    pub fn __repr__(&self) -> String {
        format!("RockCollection: {} rocks ({}, {})", self.inner.len(), self.inner.reference_plane, self.inner.origin)
    }

    /// The states of all rocks as an `(n, 6)` array of `x, y, z` (AU) and `vx, vy, vz`
    /// (AU/day). A copy: changing it does not change the collection.
    #[getter]
    pub fn states<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        let n = self.inner.len();
        let flat: Vec<f64> = self.inner.states.as_flattened().to_vec();
        numpy::ndarray::Array2::from_shape_vec((n, 6), flat).unwrap().into_pyarray(py)
    }

    #[getter]
    pub fn x<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.column(py, 0)
    }

    #[getter]
    pub fn y<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.column(py, 1)
    }

    #[getter]
    pub fn z<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.column(py, 2)
    }

    #[getter]
    pub fn vx<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.column(py, 3)
    }

    #[getter]
    pub fn vy<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.column(py, 4)
    }

    #[getter]
    pub fn vz<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.column(py, 5)
    }

    pub fn r<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.mapped(py, |s, _| state::position(s).norm())
    }

    /// Names of the rocks.
    #[getter]
    pub fn name(&self) -> Vec<String> {
        self.inner.names.clone()
    }

    /// All osculating elements of every rock, computed in one pass per rock, as a dict of
    /// arrays: `a`, `e`, `q` (AU), `inc`, `node`, `arg`, `true_anomaly`, `conic_anomaly` and
    /// `mean_anomaly` (radians; anomalies that cannot be computed are NaN).
    pub fn elements<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let pop = &self.inner;
        let fields: [(&str, fn(&Elements) -> f64); 9] = [
            ("a", |e| e.a),
            ("e", |e| e.e),
            ("q", |e| e.q),
            ("inc", |e| e.inc),
            ("node", |e| e.node),
            ("arg", |e| e.arg),
            ("true_anomaly", |e| e.true_anomaly),
            ("conic_anomaly", |e| e.conic_anomaly),
            ("mean_anomaly", |e| e.mean_anomaly),
        ];
        let columns: Vec<Vec<f64>> = py.detach(|| {
            let els = pop.elements();
            fields.par_iter().map(|(_, f)| els.iter().map(f).collect()).collect()
        });
        let d = PyDict::new(py);
        for ((name, _), v) in fields.iter().zip(columns) {
            d.set_item(*name, v.into_pyarray(py))?;
        }
        Ok(d)
    }

    pub fn a<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.mapped(py, state::semi_major_axis)
    }

    pub fn q<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.mapped(py, state::perihelion)
    }

    pub fn e<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.mapped(py, state::eccentricity)
    }

    pub fn inc<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.mapped(py, |s, _| state::inclination(s))
    }

    pub fn node<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.mapped(py, |s, _| state::node(s))
    }

    pub fn arg<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.mapped(py, state::argument_of_perihelion)
    }

    pub fn true_anomaly<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.mapped(py, state::true_anomaly)
    }

    pub fn mean_anomaly<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray1<f64>>> {
        self.mapped_or_err(py, state::mean_anomaly)
    }

    pub fn conic_anomaly<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray1<f64>>> {
        self.mapped_or_err(py, state::conic_anomaly)
    }

    #[getter]
    pub fn epoch(&self) -> Vec<PyTime> {
        self.inner.epochs.iter().map(|t| PyTime { inner: t.clone() }).collect()
    }

    pub fn get(&self, name: &str) -> PyResult<PySpaceRock> {
        match self.inner.index_of(name) {
            Some(i) => Ok(PySpaceRock { inner: self.inner.get(i).unwrap() }),
            None => Err(PyValueError::new_err(format!("No rock found with name: {}", name))),
        }
    }
}

/// Pack a slice of `Apparent` values into a dict of NumPy arrays with the given shape.
pub fn apparent_dict<'py>(py: Python<'py>, apps: &[Apparent], shape: &[usize]) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new(py);
    let fields: [(&str, fn(&Apparent) -> f64); 10] = [
        ("ra", |a| a.ra),
        ("dec", |a| a.dec),
        ("ra_rate", |a| a.ra_rate),
        ("dec_rate", |a| a.dec_rate),
        ("range", |a| a.range),
        ("range_rate", |a| a.range_rate),
        ("r_helio", |a| a.r_helio),
        ("phase", |a| a.phase),
        ("elong", |a| a.elong),
        ("mag", |a| a.mag),
    ];
    for (name, f) in fields.iter() {
        let v: Vec<f64> = apps.iter().map(f).collect();
        let arr = numpy::ndarray::ArrayD::from_shape_vec(shape.to_vec(), v)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        d.set_item(*name, arr.into_pyarray(py))?;
    }
    Ok(d)
}

/// A list of `Time` objects, or an array of Julian dates in `timescale`.
pub(crate) fn extract_times(obj: &Bound<'_, PyAny>, timescale: &str) -> PyResult<Vec<Time>> {
    if let Ok(list) = obj.extract::<Vec<PyTime>>() {
        return Ok(list.into_iter().map(|t| t.inner).collect());
    }
    if let Ok(t) = obj.extract::<PyTime>() {
        return Ok(vec![t.inner]);
    }
    let jds: Vec<f64> = obj
        .extract()
        .map_err(|_| PyValueError::new_err("epochs must be a list of Time objects or an array of Julian dates"))?;
    jds.into_iter()
        .map(|jd| Time::new(jd, timescale, "jd").map_err(|e| PyValueError::new_err(e.to_string())))
        .collect()
}
