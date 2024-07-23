use pyo3::prelude::*;
use pyo3::types::PyType;
use pyo3::exceptions::PyIndexError;
use pyo3::types::PyBool;
use pyo3::types::PyList;
use pyo3::exceptions::PyValueError;
use numpy::ToPyArray;

use rayon::prelude::*;

use spacerocks::spacerock::SpaceRock;
use spacerocks::spacerock::CoordinateFrame;
use spacerocks::Detection;

use crate::py_time::time::PyTime;
use crate::py_spacerock::spacerock::PySpaceRock;
use crate::py_observing::observer::PyObserver;
use crate::py_observing::detectioncatalog::DetectionCatalog;

use std::collections::HashMap;
use std::sync::Arc;

use numpy::{PyArray1, IntoPyArray};

pub fn create_mixed_array<T: pyo3::ToPyObject>(data: Vec<Option<T>>, py: Python) -> PyResult<Py<PyArray1<PyObject>>> {

    let numpy_array: Vec<_> = data.into_iter()
            .map(|opt| match opt {
                    Some(value) => value.to_object(py),
                    None => py.None(),
                }
            ).collect();

    Ok(numpy_array.into_pyarray(py).to_owned())
    
}


#[pyclass]
pub struct RockCollection {
    pub rocks: Vec<SpaceRock>,
    pub name_hash_map: HashMap<String, usize>,
}

#[pymethods]
impl RockCollection {
    #[new]
    pub fn new() -> Self {
        RockCollection { rocks: Vec::new(), name_hash_map: HashMap::new() }
    }

    #[classmethod]
    pub fn random(_cls: &PyType, n: usize) -> Self {
        let rocks: Vec<SpaceRock> = (0..n).into_par_iter().map(|_| SpaceRock::random()).collect();
        let mut name_hash_map = HashMap::new();
        for (i, rock) in rocks.iter().enumerate() {
            // name_hash_map.insert((*rock.name).clone(), i);
            name_hash_map.insert(rock.name.to_string(), i);
        }
        RockCollection { rocks: rocks, name_hash_map: name_hash_map }
    }

    pub fn add(&mut self, rock: PyRef<PySpaceRock>) -> Result<(), PyErr> {
        // if the name is already in the hashmap, return an error
        if self.name_hash_map.contains_key(&rock.inner.name.to_string()) {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("There is already a rock with name {} in the RockCollection", rock.inner.name)));
        }
        // self.name_hash_map.insert((*rock.inner.name).clone(), self.rocks.len());
        self.name_hash_map.insert(rock.inner.name.to_string(), self.rocks.len());
        self.rocks.push(rock.inner.clone());
        Ok(())
    }


    pub fn get_by_name(&self, name: &str) -> PyResult<PySpaceRock> {
        let index = self.name_hash_map.get(name);
        if index.is_none() {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("Rock with name {} not found", name)));
        }
        let index = index.unwrap();
        Ok(PySpaceRock { inner: self.rocks[*index].clone() })
    }

    fn __getitem__(&self, index: usize) -> PyResult<PySpaceRock> {
        if index < self.rocks.len() {
            Ok(PySpaceRock { inner: self.rocks[index].clone() })
        } else {
            Err(PyIndexError::new_err("Index out of range!"))
        }
    }

    // function to filter rocks by a boolean array, and then return a new RockCollection of clones of the rocks that are True
    pub fn filter(&self, mask: &PyArray1<bool>) -> PyResult<RockCollection> {
        let mask = unsafe {mask.as_array()};
        let mut new_rocks = Vec::new();
        let mut new_name_hash_map = HashMap::new();
        for (i, rock) in self.rocks.iter().enumerate() {
            if mask[i] {
                new_rocks.push(rock.clone());
                new_name_hash_map.insert(rock.name.to_string(), new_rocks.len()-1);
            }
        }
        Ok(RockCollection { rocks: new_rocks, name_hash_map: new_name_hash_map })
    }


    pub fn calculate_orbit(&mut self) {
        self.rocks.par_iter_mut().for_each(|rock| rock.calculate_orbit());
    }

    pub fn observe(&mut self, observer: PyRef<PyObserver>) -> PyResult<DetectionCatalog> {
        let o = observer.inner.clone();

        if o.frame != CoordinateFrame::J2000 {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("Observer frame is not J2000. Cannot observe rocks.")));
        }

        // let observations: Vec<Result<PyDetection, PyErr>> = self.rocks.par_iter().map(|rock| {
        //     match rock.observe(&observer.inner) {
        //             Ok(detection) => detection,
        //             Err(e) => Err(PyValueError::new_err(format!("Error observing: {}", e))),
        //     }
        // }).collect();

        // let mut results = Vec::new();
        // let mut errors = Vec::new();

        // for observation in observations {
        //     match observation {
        //         Ok(detection) => results.push(detection),
        //         Err(err) => errors.push(err.to_string()),
        //     }
        // }

        // if !errors.is_empty() {
        //     return Err(PyValueError::new_err(format!("Errors occurred: {}", errors.join(", "))));
        // }

        // Ok(DetectionCatalog { observations: results })

        let observations: Vec<_> = self.rocks.par_iter_mut().map(|rock| rock.observe(&o).unwrap()).collect();   
        Ok(DetectionCatalog { observations: observations })     
    }

    pub fn analytic_propagate(&mut self, t: PyRef<PyTime>) {
        let ep = &t.inner;
        self.rocks.par_iter_mut().for_each(|rock| rock.analytic_propagate(ep));
    }

    pub fn change_frame(&mut self, frame: &str) {
        self.rocks.par_iter_mut().for_each(|rock| rock.change_frame(frame).expect("Failed to change frame"));
    }

    #[getter]
    pub fn frame(&self) -> Vec<String> {
        let frames = self.rocks.par_iter().map(|rock| rock.frame.to_string()).collect::<Vec<String>>();
        frames
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
        x.into_pyarray(py).to_owned()
    }

    #[getter]
    pub fn y(&self, py: Python) -> Py<PyArray1<f64>> {
        let y: Vec<f64> = self.rocks.par_iter().map(|rock| rock.position[1]).collect();
        y.into_pyarray(py).to_owned()
    }

    #[getter]
    pub fn z(&self, py: Python) -> Py<PyArray1<f64>> {
        let z: Vec<f64> = self.rocks.par_iter().map(|rock| rock.position[2]).collect();
        z.into_pyarray(py).to_owned()
    }

    #[getter]
    pub fn vx(&self, py: Python) -> Py<PyArray1<f64>> {
        let vx: Vec<f64> = self.rocks.par_iter().map(|rock| rock.velocity[0]).collect();
        vx.into_pyarray(py).to_owned()
    }

    #[getter]
    pub fn vy(&self, py: Python) -> Py<PyArray1<f64>> {
        let vy: Vec<f64> = self.rocks.par_iter().map(|rock| rock.velocity[1]).collect();
        vy.into_pyarray(py).to_owned()
    }

    #[getter]
    pub fn vz(&self, py: Python) -> Py<PyArray1<f64>> {
        let vz: Vec<f64> = self.rocks.par_iter().map(|rock| rock.velocity[2]).collect();
        vz.into_pyarray(py).to_owned()
    }

    #[getter]
    pub fn r(&self, py: Python) -> Py<PyArray1<f64>> {
        let r: Vec<f64> = self.rocks.par_iter().map(|rock| rock.r()).collect();
        r.into_pyarray(py).to_owned()
    }

    // #[getter]
    // pub fn name(&self) -> Vec<String> {
    //     self.rocks.par_iter().map(|rock| rock.name.clone()).collect()
    // }
    #[getter] 
    pub fn name(&self, py: Python) -> PyResult<Py<PyArray1<PyObject>>> {
        let names: Vec<Option<String>> = self.rocks.par_iter().map(|rock| Some((*rock.name).clone())).collect();
        create_mixed_array(names, py)
    }


    #[getter]
    pub fn node(&self, py: Python) -> Py<PyArray1<f64>> {
        // self.rocks.par_iter().map(|rock| rock.orbit.as_ref().unwrap().node).collect()
        let nodes: Vec<f64> = self.rocks.par_iter().map(|rock| rock.orbit.as_ref().unwrap().node).collect();
        nodes.into_pyarray(py).to_owned()
    }

    #[getter]
    pub fn inc(&self, py: Python) -> Py<PyArray1<f64>> {
        // self.rocks.par_iter().map(|rock| rock.orbit.as_ref().unwrap().inc).collect()
        let incs: Vec<f64> = self.rocks.par_iter().map(|rock| rock.orbit.as_ref().unwrap().inc).collect();
        incs.into_pyarray(py).to_owned()
    }

    #[getter]
    pub fn e(&self, py: Python) -> Py<PyArray1<f64>> {
        // self.rocks.par_iter().map(|rock| rock.orbit.as_ref().unwrap().e).collect()
        let es: Vec<f64> = self.rocks.par_iter().map(|rock| rock.orbit.as_ref().unwrap().e).collect();
        es.into_pyarray(py).to_owned()
    }

    #[getter]
    pub fn a(&self, py: Python) -> Py<PyArray1<f64>> {
        // self.rocks.par_iter().map(|rock| rock.orbit.as_ref().unwrap().a).collect()
        let a_values: Vec<f64> = self.rocks.par_iter().map(|rock| rock.orbit.as_ref().unwrap().a).collect();
        a_values.into_pyarray(py).to_owned()
    }

    #[getter]
    pub fn arg(&self, py: Python) -> Py<PyArray1<f64>> {
        // self.rocks.par_iter().map(|rock| rock.orbit.as_ref().unwrap().arg).collect()
        let args: Vec<f64> = self.rocks.par_iter().map(|rock| rock.orbit.as_ref().unwrap().arg).collect();
        args.into_pyarray(py).to_owned()
    }

    #[getter]
    pub fn varpi(&self, py: Python) -> Py<PyArray1<f64>> {
        // self.rocks.par_iter().map(|rock| rock.orbit.as_ref().unwrap().varpi()).collect()
        let varpis: Vec<f64> = self.rocks.par_iter().map(|rock| rock.orbit.as_ref().unwrap().varpi()).collect();
        varpis.into_pyarray(py).to_owned()
    }

    // #[getter]
    // pub fn epoch(&self) -> Vec<f64> {
    //     // self.rocks.par_iter().map(|rock| rock.epoch).collect()
    //     let epochs: Vec<f64> = self.rocks.par_iter().map(|rock| rock.epoch).collect();

    // }

}
