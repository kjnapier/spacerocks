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
use crate::py_observing::observer::PyObserver;
use crate::py_observing::observation::PyObservation;

use numpy::{PyArray1, IntoPyArray};

use std::path::{PathBuf};
// use std::fs;
// use serde::{Serialize, Deserialize};
// use nalgebra::Vector3;
use serde_json;
// use arrow::array::{Float64Array, StringArray};
// use crate::mpc::MPCHandler;
use dirs::home_dir;


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
#[pyclass]
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

    pub fn observe(&mut self, observer: PyRef<PyObserver>) -> PyResult<Vec<PyObservation>> {
        let o = observer.inner.clone();

        if o.reference_plane != ReferencePlane::J2000 {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("Observer frame is not J2000. Cannot observe rocks.")));
        }

        let observations: Vec<_> = self.rocks.par_iter_mut().map(|rock| rock.observe(&o).unwrap()).collect();   
        let py_observations: Vec<_> = observations.into_iter().map(|obs| PyObservation { inner: obs }).collect();
        Ok(py_observations)
           
    }

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

    pub fn change_reference_plane(&mut self, reference_plane: &str) {
        self.rocks.par_iter_mut().for_each(|rock| rock.change_reference_plane(reference_plane).expect("Failed to change frame"));
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