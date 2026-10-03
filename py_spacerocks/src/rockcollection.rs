use pyo3::prelude::*;
use pyo3::exceptions::{PyIndexError,PyValueError};
use pyo3::{PyErr, Python, PyResult};
// use pyo3::types::{PyList, PyType, IntoPyDict};
use rayon::prelude::*;
// use pyo3::types::PySequence;

use spacerocks::spacerock::SpaceRock;
// use spacerocks::Time;
use spacerocks::ReferencePlane;

use crate::py_time::time::PyTime;
use crate::PySpaceRock;
use crate::py_spice::spicekernel::PySpiceKernel;
use crate::py_observing::observer::PyObserver;
use crate::py_observing::observation::PyObservation;

use numpy::{PyArray1, IntoPyArray, PyArray2};
use pyo3::types::PyDict;
use spacerocks::observing::{Apparent, Observer};
use spacerocks::batch::{self, BatchOptions, Method};
use spacerocks::Time;
use crate::py_observing::observatory::PyObservatory;

// use std::fs;
// use serde::{Serialize, Deserialize};
// use nalgebra::Vector3;
// use arrow::array::{Float64Array, StringArray};
// use crate::mpc::MPCHandler;
// use ndarray;


// use pyo3::impl_::pymethods::AsyncIterBaseKind;

// use numpy::{PyArray1, IntoPyArray, PyArray};

// pub fn create_mixed_array<T: pyo3::IntoPyObject>(data: Vec<Option<T>>, py: Python) -> pyo3::Bound<'_, PyArray<pyo3::Py<PyAny>, numpy::ndarray::Dim<[usize; 1]>>> {
//     let numpy_array: Vec<_> = data.into_iter()
//             .map(|opt| match opt {
//                     Some(value) => value.to_object(py),
//                     None => py.None(),
//                 }
//             ).collect();
//     numpy_array.into_pyarray(py).to_owned()
// }

/// Represents a collection of space rocks.
///
/// This struct is used to manage and manipulate a collection of 
/// `SpaceRock` objects, including operations such as filtering, 
/// observing, and converting formats.
#[pyclass(from_py_object)]
#[derive(Clone)]
pub struct RockCollection {
    /// A vector holding all `SpaceRock` instances.
    pub rocks: Vec<SpaceRock>,
}

#[pymethods]
impl RockCollection {
    /// Creates a new, empty `RockCollection`.
    
    #[new]
    pub fn new() -> Self {
        RockCollection { rocks: Vec::new() }
    }

    /// Constructs a `RockCollection` from MPC data.
    ///
    /// This method fetches and reads data from the Minor Planet Center (MPC)
    /// and constructs a `RockCollection` from the data.
    ///
    /// # Arguments
    /// * `mpc_path` - The path to the directory where the MPC data will be stored.
    /// * `catalog` - The name of the MPC catalog to fetch (i.e, mpcorb_extended).
    /// * `download_data` - A boolean flag indicating whether to download the data if it is not already present.
    ///
    /// # Returns
    /// A `RockCollection` instance.
    ///
    /// # Example
    /// ```python
    /// from spacerocks import RockCollection
    ///
    /// rocks = RockCollection.from_mpc("data/mpc", "mpcorb_extended", download_data=True)
    /// ```
    // #[staticmethod]
    // #[pyo3(signature = (catalog, download_data=false, mpc_path=None, orbit_type=None))]
    // pub fn from_mpc(catalog: String, download_data: bool, mpc_path: Option<PathBuf>, orbit_type: Option<String>) -> PyResult<Self> {
    //     let default_path = home_dir()
    //         .unwrap_or_default()
    //         .join(".spacerocks")
    //         .join("mpc");
        
    //     let final_path = mpc_path.unwrap_or(default_path);

    //     MPCHandler::create_rock_collection(
    //         final_path,
    //         catalog,
    //         download_data,
    //         orbit_type 
    //     )
    // }

    

    // #[classmethod]
    // pub fn random(_cls: &PyType, n: usize) -> Self {
    //     let rocks: Vec<SpaceRock> = (0..n).into_par_iter().map(|_| SpaceRock::random()).collect();
    //     RockCollection { rocks: rocks }
    // }


    pub fn add(&mut self, rock: PyRef<PySpaceRock>) {
        self.rocks.push(rock.inner.clone());
    }


    fn __getitem__(&self, index: usize) -> PyResult<PySpaceRock> {
        if index < self.rocks.len() {
            Ok(PySpaceRock { inner: self.rocks[index].clone() })
        } else {
            Err(PyIndexError::new_err("Index out of range!"))
        }
    }


    // function to filter rocks by a boolean array, and then return a new RockCollection of clones of the rocks that are True
    pub fn filter(&self, indices: Vec<bool>) -> PyResult<Self> {
        if indices.len() != self.rocks.len() {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                "Mask length must match the number of rocks.",
            ));
        }

        let filtered_rocks = self
            .rocks
            .iter()
            .zip(indices.iter())
            .filter_map(|(rock, &keep)| if keep { Some(rock.clone()) } else { None })
            .collect();

        Ok(Self {
            rocks: filtered_rocks,
        })
    }


    // pub fn calculate_orbit(&mut self) {
    //     self.rocks.par_iter_mut().for_each(|rock| rock.calculate_orbit());
    // }

    pub fn observe(&mut self, py: Python<'_>, observer: PyRef<PyObserver>) -> PyResult<Vec<PyObservation>> {
        let o = observer.inner.clone();

        // if o.reference_plane != ReferencePlane::J2000 {
        //     return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("Observer frame is not J2000. Cannot observe rocks.")));
        // }

        let rocks = &mut self.rocks;
        let observations = py
            .detach(|| {
                rocks
                    .par_iter_mut()
                    .map(|rock| rock.observe(&o).map_err(|e| format!("{}: {}", rock.name, e)))
                    .collect::<Result<Vec<_>, String>>()
            })
            .map_err(PyValueError::new_err)?;
        let py_observations: Vec<_> = observations.into_iter().map(|obs| PyObservation { inner: obs }).collect();
        Ok(py_observations)
           
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
        let rocks = &self.rocks;
        let apps = py
            .detach(|| {
                rocks
                    .par_iter()
                    .map(|rock| rock.apparent(o).map_err(|e| format!("{}: {}", rock.name, e)))
                    .collect::<Result<Vec<_>, String>>()
            })
            .map_err(PyValueError::new_err)?;
        apparent_dict(py, &apps, &[apps.len()])
    }

    pub fn calc_radec<'py>(
        &self,
        py: Python<'py>,
        observer: PyRef<PyObserver>,
    ) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let o = &observer.inner;

        if o.reference_plane != ReferencePlane::J2000 {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                "Observer frame is not J2000. Cannot calculate RA/Dec.",
            ));
        }

        let n = self.rocks.len();
        let mut data = vec![0.0f64; 2 * n];
        let rocks = &self.rocks;
        py.detach(|| {
            data.par_chunks_mut(2)
                .zip(rocks.par_iter())
                .try_for_each(|(row, rock)| {
                    let (ra, dec) = rock.calc_radec(o).map_err(|e| format!("{}: {}", rock.name, e))?;
                    row[0] = ra;
                    row[1] = dec;
                    Ok::<(), String>(())
                })
        })
        .map_err(PyValueError::new_err)?;

        let arr = numpy::ndarray::Array2::from_shape_vec((n, 2), data).unwrap();
        Ok(arr.into_pyarray(py))
    }


    // pub fn calc_radec_no_light_time<'py>(
    //     &self,
    //     py: Python<'py>,
    //     observer: PyRef<PyObserver>,
    // ) -> PyResult<Bound<'py, PyArray2<f64>>> {
    //     let o = &observer.inner;

    //     if o.reference_plane != ReferencePlane::J2000 {
    //         return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
    //             "Observer frame is not J2000. Cannot calculate RA/Dec.",
    //         ));
    //     }

    //     let n = self.rocks.len();
    //     let mut data = vec![0.0f64; 2 * n];

    //     data.par_chunks_mut(2)
    //         .zip(self.rocks.par_iter())
    //         .for_each(|(row, rock)| {
    //             let (ra, dec) = rock.calc_radec_no_light_time(o).unwrap();
    //             row[0] = ra;
    //             row[1] = dec;
    //         });

    //     let arr = ndarray::Array2::from_shape_vec((n, 2), data).unwrap();
    //     Ok(arr.into_pyarray(py))
    // }




    // pub fn calc_radec_no_light_time(&self, observer: PyRef<PyObserver>) -> PyResult<Vec<(f64, f64)>> {
    //     let o = &observer.inner;

    //     if o.reference_plane != ReferencePlane::J2000 {
    //         return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
    //             "Observer frame is not J2000. Cannot calculate RA/Dec."
    //         ));
    //     }

    //     let radec_values: Vec<_> = self
    //         .rocks
    //         .par_iter()
    //         .map(|rock| rock.calc_radec_no_light_time(o).unwrap())
    //         .collect();

    //     Ok(radec_values)
    // }

    // pub fn calc_radec_no_light_time(&mut self, observer: PyRef<PyObserver>) -> PyResult<Vec<(f64, f64)>> {
    //     let o = observer.inner.clone();

    //     if o.reference_plane != ReferencePlane::J2000 {
    //         return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("Observer frame is not J2000. Cannot calculate RA/Dec.")));
    //     }

    //     let radec_values: Vec<_> = self.rocks.par_iter_mut().map(|rock| rock.calc_radec_no_light_time(&o).unwrap()).collect();   
    //     Ok(radec_values)
           
    // }

    pub fn analytic_propagate(&mut self, epoch: PyRef<PyTime>) -> PyResult<()> {
        let ep = &epoch.inner;
    
        if let Some(error) = self.rocks
            .par_iter_mut()
            .filter_map(|rock| {
                match rock.analytic_propagate(ep) {
                    Err(e) => Some(format!("Failed to propagate rock: {}", e)),
                    Ok(_) => None
                }
            })
            .find_first(|_| true) {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(error));
        }
        
        Ok(())
    }

    pub fn change_reference_plane(&mut self, reference_plane: &str) -> PyResult<()> {
        ReferencePlane::from_str(reference_plane).map_err(PyValueError::new_err)?;
        self.rocks
            .par_iter_mut()
            .try_for_each(|rock| rock.change_reference_plane(reference_plane).map_err(|e| e.to_string()))
            .map_err(PyValueError::new_err)
    }

    /// Propagate every rock to `epoch`.
    ///
    /// `method="nbody"` (default) integrates with IAS15 in the field of the Sun, planets, Moon,
    /// Pluto and 16 massive asteroids from `kernel`, with ASSIST's force model (relativity,
    /// Earth and solar harmonics, and non-gravitational forces for rocks with `set_nongrav`); `method="twobody"` uses Keplerian motion
    /// (then `kernel` is not used). Rocks are integrated together in groups of up to
    /// `chunk_size` (grouped by epoch and perihelion distance), with groups run in parallel;
    /// `chunk_size=1` integrates every rock on its own, exactly like `SpaceRock.propagate`.
    /// Each rock keeps its reference plane and origin (SUN or SSB).
    #[pyo3(signature = (epoch, kernel, method = "nbody", chunk_size = 64))]
    pub fn propagate(&mut self, py: Python<'_>, epoch: PyRef<PyTime>, kernel: PyRef<PySpiceKernel>, method: &str, chunk_size: usize) -> PyResult<()> {
        let opts = BatchOptions { method: Method::from_str(method).map_err(PyValueError::new_err)?, chunk_size, ..Default::default() };
        let ep = epoch.inner.clone();
        let k = &kernel.inner;
        let rocks = &mut self.rocks;
        py.detach(|| batch::propagate_batch(rocks, &ep, k, &opts).map_err(|e| e.to_string()))
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
        let rocks = &self.rocks;
        let eph = py
            .detach(|| batch::ephemeris(rocks, &observers, k, &opts).map_err(|e| e.to_string()))
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

    #[getter]
    pub fn reference_plane(&self) -> Vec<String> {
        let reference_planes = self.rocks.par_iter().map(|rock| rock.reference_plane.to_string()).collect::<Vec<String>>();
        reference_planes
    }

    
    fn __len__(&self) -> usize {
        self.rocks.len()
    }

    pub fn len(&self) -> usize {
        self.rocks.len()
    }

    pub fn __repr__(&self) -> String {
        format!("RockCollection: {} rocks", self.rocks.len())
    }

    #[getter]
    pub fn x(&self, py: Python) -> Py<PyArray1<f64>> {
        let x: Vec<f64> = self.rocks.par_iter().map(|rock| rock.position[0]).collect();
        x.into_pyarray(py).to_owned().into()
    }

    #[getter]
    pub fn y(&self, py: Python) -> Py<PyArray1<f64>> {
        let y: Vec<f64> = self.rocks.par_iter().map(|rock| rock.position[1]).collect();
        y.into_pyarray(py).to_owned().into()
    }

    #[getter]
    pub fn z(&self, py: Python) -> Py<PyArray1<f64>> {
        let z: Vec<f64> = self.rocks.par_iter().map(|rock| rock.position[2]).collect();
        z.into_pyarray(py).to_owned().into()
    }

    #[getter]
    pub fn vx(&self, py: Python) -> Py<PyArray1<f64>> {
        let vx: Vec<f64> = self.rocks.par_iter().map(|rock| rock.velocity[0]).collect();
        vx.into_pyarray(py).to_owned().into()
    }

    #[getter]
    pub fn vy(&self, py: Python) -> Py<PyArray1<f64>> {
        let vy: Vec<f64> = self.rocks.par_iter().map(|rock| rock.velocity[1]).collect();
        vy.into_pyarray(py).to_owned().into()
    }

    #[getter]
    pub fn vz(&self, py: Python) -> Py<PyArray1<f64>> {
        let vz: Vec<f64> = self.rocks.par_iter().map(|rock| rock.velocity[2]).collect();
        vz.into_pyarray(py).to_owned().into()
    }

    // #[getter]
    pub fn r(&self, py: Python) -> Py<PyArray1<f64>> {
        let r: Vec<f64> = self.rocks.par_iter().map(|rock| rock.r()).collect();
        r.into_pyarray(py).to_owned().into()
    }

    // #[getter]
    // pub fn name(&self) -> Vec<String> {
    //     self.rocks.par_iter().map(|rock| rock.name.clone()).collect()
    // }

    // #[getter] 
    // pub fn name(&self, py: Python) -> PyResult<Py<PyArray1<PyObject>>> {
    //     let names: Vec<Option<String>> = self.rocks.par_iter().map(|rock| Some((*rock.name).clone())).collect();
    //     create_mixed_array(names, py)
    // }

    pub fn a(&self, py: Python) -> Py<PyArray1<f64>> {
        let a_values: Vec<f64> = self.rocks.par_iter().map(|rock| rock.a()).collect();
        a_values.into_pyarray(py).to_owned().into()
    }

    pub fn q(&self, py: Python) -> Py<PyArray1<f64>> {
        let q_values: Vec<f64> = self.rocks.par_iter().map(|rock| rock.q()).collect();
        q_values.into_pyarray(py).to_owned().into()
    }

    pub fn e(&self, py: Python) -> Py<PyArray1<f64>> {
        let e_values: Vec<f64> = self.rocks.par_iter().map(|rock| rock.e()).collect();
        e_values.into_pyarray(py).to_owned().into()
    }

    pub fn inc(&self, py: Python) -> Py<PyArray1<f64>> {
        let inc_values: Vec<f64> = self.rocks.par_iter().map(|rock| rock.inc()).collect();
        inc_values.into_pyarray(py).to_owned().into()
    }

    pub fn node(&self, py: Python) -> Py<PyArray1<f64>> {
        let node_values: Vec<f64> = self.rocks.par_iter().map(|rock| rock.node()).collect();
        node_values.into_pyarray(py).to_owned().into()
    }

    pub fn arg(&self, py: Python) -> Py<PyArray1<f64>> {
        let arg_values: Vec<f64> = self.rocks.par_iter().map(|rock| rock.arg()).collect();
        arg_values.into_pyarray(py).to_owned().into()
    }

    pub fn true_anomaly(&self, py: Python) -> Py<PyArray1<f64>> {
        let true_anomaly_values: Vec<f64> = self.rocks.par_iter().map(|rock| rock.true_anomaly()).collect();
        true_anomaly_values.into_pyarray(py).to_owned().into()
    }

    pub fn mean_anomaly(&self, py: Python) -> Py<PyArray1<f64>> {
        let mean_anomaly_values: Vec<f64> = self.rocks.par_iter().map(|rock| rock.mean_anomaly()).collect();
        mean_anomaly_values.into_pyarray(py).to_owned().into()
    }

    pub fn conic_anomaly(&self, py: Python) -> Py<PyArray1<f64>> {
        let conic_anomaly_values: Vec<f64> = self.rocks.par_iter().map(|rock| rock.conic_anomaly()).collect();
        conic_anomaly_values.into_pyarray(py).to_owned().into()
    }


    #[getter]
    pub fn epoch(&self) -> Vec<PyTime> {
        self.rocks.par_iter().map(|rock| PyTime { inner: rock.epoch.clone() }).collect()
    }

    pub fn get(&self, name: &str) -> PyResult<PySpaceRock> {
        let rock = self.rocks.iter().find(|rock| rock.name == name);
        match rock {
            Some(rock) => Ok(PySpaceRock { inner: rock.clone() }),
            None => Err(PyErr::new::<PyValueError, _>(
                format!("No rock found with name: {}", name)
            ))
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
