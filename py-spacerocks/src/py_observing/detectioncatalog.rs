use pyo3::prelude::*;

use rayon::prelude::*;

use spacerocks::Detection;

use pyo3::exceptions::PyIndexError;

use numpy::{PyArray1, IntoPyArray};

use crate::py_observing::detection::PyDetection;

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
pub struct DetectionCatalog {
   pub observations: Vec<Detection>,
}

#[pymethods]
impl DetectionCatalog {

    #[new]
    pub fn new() -> Self {
        DetectionCatalog { observations: Vec::new() }
    }

    pub fn add(&mut self, observation: PyRef<PyDetection>) -> Result<(), PyErr> {
        self.observations.push(observation.inner.clone());
        Ok(())
    }

    // make indexable
    fn __getitem__(&self, index: usize) -> PyResult<PyDetection> {
        if index < self.observations.len() {
            Ok(PyDetection { inner: self.observations[index].clone() })
        } else {
            Err(PyIndexError::new_err("Index out of range!"))
        }
    }

    pub fn len(&self) -> usize {
        self.observations.len()
    }

    pub fn __repr__(&self) -> String {
        format!("DetectionCatalog: {} detections", self.observations.len())
    }

    #[getter]
    pub fn ra(&self, py: Python) -> Py<PyArray1<f64>> {
        let ra: Vec<f64> = self.observations.par_iter().map(|obs| obs.ra).collect();
        ra.into_pyarray(py).to_owned()
    }

    #[getter]
    pub fn dec(&self, py: Python) -> Py<PyArray1<f64>> {
        let dec: Vec<f64> = self.observations.par_iter().map(|obs| obs.dec).collect();
        dec.into_pyarray(py).to_owned()
    }

    #[getter]
    pub fn ra_rate(&self, py: Python) -> PyResult<Py<PyArray1<PyObject>>> {
        let ra_rates: Vec<Option<f64>> = self.observations.par_iter().map(|obs| obs.ra_rate).collect();
        create_mixed_array(ra_rates, py)
    }

    #[getter]
    pub fn dec_rate(&self, py: Python) -> PyResult<Py<PyArray1<PyObject>>> {
        let dec_rates: Vec<Option<f64>> = self.observations.par_iter().map(|obs| obs.dec_rate).collect();
        create_mixed_array(dec_rates, py)
    }

    #[getter]
    pub fn rho(&self, py: Python) -> PyResult<Py<PyArray1<PyObject>>> {
        let rhos: Vec<Option<f64>> = self.observations.par_iter().map(|obs| obs.rho).collect();
        create_mixed_array(rhos, py)
    }

    #[getter]
    pub fn rho_rate(&self, py: Python) -> PyResult<Py<PyArray1<PyObject>>> {
        let rho_rates: Vec<Option<f64>> = self.observations.par_iter().map(|obs| obs.rho_rate).collect();
        create_mixed_array(rho_rates, py)
    }

    #[getter]
    pub fn mag(&self, py: Python) -> PyResult<Py<PyArray1<PyObject>>> {
        let mags: Vec<Option<f64>> = self.observations.par_iter().map(|obs| obs.mag).collect();
        create_mixed_array(mags, py)
    }

    #[getter]
    pub fn epoch(&self, py: Python) -> Py<PyArray1<f64>> {
        let epochs: Vec<f64> = self.observations.par_iter().map(|obs| obs.epoch.epoch).collect();
        epochs.into_pyarray(py).to_owned()
    }

    pub fn calc_altaz(&self, py: Python) -> PyResult<(Py<PyArray1<f64>>, Py<PyArray1<f64>>)> {
        
        // process the altaz in parallel, but hanfle errors
        let altaz: Vec<_> = self.observations.par_iter().map(|obs| obs.calc_altaz().unwrap()).collect();
        // let mut alts = altaz.par_iter().map(|(alt, az)| *alt).collect::<Vec<f64>>();
        // let mut azs = altaz.par_iter().map(|(alt, az)| *az).collect::<Vec<f64>>();
        let alts: Vec<_> = altaz.par_iter().map(|(alt, _az)| *alt).collect();
        let azs: Vec<_> = altaz.par_iter().map(|(_alt, az)| *az).collect();
        
        Ok((alts.into_pyarray(py).to_owned(), azs.into_pyarray(py).to_owned()))
    }
    


}
