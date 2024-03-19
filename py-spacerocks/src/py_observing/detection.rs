use pyo3::prelude::*;
use pyo3::types::PyType;
use pyo3::exceptions::PyValueError;

use spacerocks::observing::Detection;
use crate::py_observing::observer::PyObserver;
use crate::py_time::time::PyTime;

#[pyclass]
#[pyo3(name = "Detection")]
#[derive(Clone)]
pub struct PyDetection {
    pub inner: Detection,
}

#[pymethods]
impl PyDetection {

    #[new]
    pub fn new(name: String, ra: f64, dec: f64, epoch: PyRef<PyTime>, observer: PyRef<PyObserver>, mag: Option<f64>, 
               filter: Option<String>, ra_uncertainty: Option<f64>, dec_uncertainty: Option<f64>) -> Self {
        let inner = Detection {
            name: name.into(),
            epoch: epoch.inner.clone(),
            observer: observer.inner.clone(),
            ra,
            dec,
            ra_rate: None,
            dec_rate: None,
            rho: None,
            rho_rate: None,
            mag,
            filter,
            ra_uncertainty: ra_uncertainty,
            dec_uncertainty: dec_uncertainty,
            ra_rate_uncertainty: None,
            dec_rate_uncertainty: None,
            rho_uncertainty: None,
            rho_rate_uncertainty: None,
            mag_uncertainty: None,
        };

        PyDetection { inner }
    }


    #[classmethod]
    pub fn streak(_cls: &PyType, name: String, ra: f64, dec: f64, ra_rate: f64, dec_rate: f64, epoch: PyRef<PyTime>, observer: PyRef<PyObserver>, mag: Option<f64>, filter: Option<String>) -> Self {
        let inner = Detection {
            name: name.into(),
            epoch: epoch.inner.clone(),
            observer: observer.inner.clone(),
            ra,
            dec,
            ra_rate: Some(ra_rate),
            dec_rate: Some(dec_rate),
            rho: None,
            rho_rate: None,
            mag,
            filter,
            ra_uncertainty: None,
            dec_uncertainty: None,
            ra_rate_uncertainty: None,
            dec_rate_uncertainty: None,
            rho_uncertainty: None,
            rho_rate_uncertainty: None,
            mag_uncertainty: None,
        };

        PyDetection { inner }
    }

    fn __repr__(&self) -> String {
        format!("Detection: {} at time: {:?}", self.inner.name, self.inner.epoch)
    }

    #[getter]
    fn name(&self) -> &str {
        &self.inner.name
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

    #[getter]
    fn epoch(&self) -> PyResult<PyTime> {
        Ok(PyTime { inner: self.inner.epoch.clone() })
    }

    #[getter]
    fn ra_uncertainty(&self) -> Option<f64> {
        self.inner.ra_uncertainty
    }

    #[getter]
    fn dec_uncertainty(&self) -> Option<f64> {
        self.inner.dec_uncertainty
    }

    #[getter]
    fn ra_rate_uncertainty(&self) -> Option<f64> {
        self.inner.ra_rate_uncertainty
    }

    #[getter]
    fn dec_rate_uncertainty(&self) -> Option<f64> {
        self.inner.dec_rate_uncertainty
    }

    #[getter]
    fn rho_uncertainty(&self) -> Option<f64> {
        self.inner.rho_uncertainty
    }

    #[getter]
    fn rho_rate_uncertainty(&self) -> Option<f64> {
        self.inner.rho_rate_uncertainty
    }

    fn calc_altaz(&self) -> PyResult<(f64, f64)> {
        // let altaz = self.inner.calc_altaz();
        let altaz = match self.inner.calc_altaz() {
            Ok(altaz) => altaz,
            // Err(e) => return Err(PyValueError::new_err("Could not calculate altaz")),
            Err(e) => return Err(PyValueError::new_err(e.to_string())),
        };
        Ok(altaz)
    }

    
}