use pyo3::prelude::*;
use pyo3::types::PyType;
use pyo3::exceptions::PyValueError;

use nalgebra::Vector3;
use spacerocks::observing::Detection;
use crate::observing::observer::PyObserver;

#[pyclass]
#[pyo3(name = "Detection")]
pub struct PyDetection {
    pub inner: Detection,
}

#[pymethods]
impl PyDetection {

    fn __repr__(&self) -> String {
        format!("Detection: {} at time: {:?}", self.inner.name, self.inner.epoch)
    }

    #[getter]
    fn ra(&self) -> f64 {
        self.inner.ra
    }

    #[getter]
    fn dec(&self) -> f64 {
        self.inner.dec
    }

    #[getter]
    fn ra_rate(&self) -> Option<f64> {
        self.inner.ra_rate
    }

    #[getter]
    fn dec_rate(&self) -> Option<f64> {
        self.inner.dec_rate
    }

    #[getter]
    fn rho(&self) -> Option<f64> {
        self.inner.rho
    }

    #[getter]
    fn rho_rate(&self) -> Option<f64> {
        self.inner.rho_rate
    }

    #[getter]
    fn mag(&self) -> Option<f64> {
        self.inner.mag
    }

    #[getter]
    fn filter(&self) -> Option<String> {
        let f = match &self.inner.filter {
            Some(f) => f.clone(),
            None => return None,
        };
        Some(f)
    }    

    #[getter]
    fn observer(&self) -> PyResult<PyObserver> {
        Ok(PyObserver { inner: self.inner.observer.clone() })
    }

    fn calc_altaz(&self) -> PyResult<(f64, f64)> {
        let altaz = self.inner.calc_altaz();
        let altaz = match self.inner.calc_altaz() {
            Ok(altaz) => altaz,
            Err(e) => return Err(PyValueError::new_err("Could not calculate altaz")),
        };
        Ok(altaz)
    }

    
}