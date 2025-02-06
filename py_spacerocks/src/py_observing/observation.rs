use pyo3::prelude::*;
use pyo3::types::PyType;

use spacerocks::Observation;

use crate::py_observing::observer::PyObserver;
use crate::py_time::time::PyTime;
// use crate::py_spacerock::origin::PyOrigin;

use nalgebra::DMatrix;

#[pyclass]
#[pyo3(name = "Observation")]
pub struct PyObservation {
    pub inner: Observation,
}

#[pymethods]
impl PyObservation {

    #[classmethod]
    #[pyo3(signature = (epoch, ra, dec, observer, covariance=None, mag=None, mag_err=None))]
    fn from_astrometry(_cls: Py<PyType>, epoch: PyTime, ra: f64, dec: f64, observer: PyObserver, covariance: Option<[[f64; 2]; 2]>, mag: Option<f64>, mag_err: Option<f64>) -> PyResult<PyObservation> {
        let obs = Observation::from_astrometry(epoch.inner, ra, dec, observer.inner, covariance, mag, mag_err);
        if obs.is_err() {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("Failed to create Observation from Astrometry")));
        }
        Ok(PyObservation { inner: obs.unwrap() })
    }

    #[classmethod]
    #[pyo3(signature = (epoch, ra, dec, ra_rate, dec_rate, observer, covariance=None, mag=None, mag_err=None))]
    fn from_streak(_cls: Py<PyType>, epoch: PyTime, ra: f64, dec: f64, ra_rate: f64, dec_rate: f64, observer: PyObserver, covariance: Option<[[f64; 4]; 4]>, mag: Option<f64>, mag_err: Option<f64>) -> PyResult<PyObservation> {
        let obs = Observation::from_streak(epoch.inner, ra, dec, ra_rate, dec_rate, observer.inner, covariance, mag, mag_err);
        if obs.is_err() {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("Failed to create Observation from Streak")));
        }
        Ok(PyObservation { inner: obs.unwrap() })
    }

    // repr
    fn __repr__(&self) -> PyResult<String> {
        let mut components = Vec::new();
        
        components.push("Observation:".to_string());
        
        // Format positional components with 4 decimal places
        components.push(format!("(RA={:.4}°, Dec={:.4}°)", 
            self.inner.ra().to_degrees(), 
            self.inner.dec().to_degrees()));

        // Add motion
        if let (Some(ra_rate), Some(dec_rate)) = (self.inner.ra_rate(), self.inner.dec_rate()) {
            components.push(format!("ra/dec rates = ({:.4}°/day, {:.4}°/day)", 
                ra_rate.to_degrees(), 
                dec_rate.to_degrees()));
        }

        // Range
        if let (Some(r), Some(r_rate)) = (self.inner.range(), self.inner.range_rate()) {
            components.push(format!("distance: {:.2}AU (rate: {:.4}AU/day)", 
                r, r_rate));
        }

        // Magnitude
        if let Some(m) = self.inner.mag() {
            components.push(format!("mag: {:.2}", m));
        }

        // Add epoch (always present)
        components.push(format!("epoch: {}", self.inner.epoch));

        // Join all components with newlines
        Ok(components.join("\n"))
    }

    #[getter]
    fn ra(&self) -> f64 {
        self.inner.ra()
    }

    #[getter]
    fn dec(&self) -> f64 {
        self.inner.dec()
    }

    #[getter]
    fn ra_rate(&self) -> Option<f64> {
        self.inner.ra_rate()
    }

    #[getter]
    fn dec_rate(&self) -> Option<f64> {
        self.inner.dec_rate()
    }

    #[getter]
    fn range(&self) -> Option<f64> {
        self.inner.range()
    }

    #[getter]
    fn range_rate(&self) -> Option<f64> {
        self.inner.range_rate()
    }

    #[getter]
    fn mag(&self) -> Option<f64> {
        self.inner.mag()
    }

    #[getter]
    fn epoch(&self) -> PyTime {
        PyTime { inner: self.inner.epoch.clone() }
    }

    #[getter]
    fn observer(&self) -> PyObserver {
        PyObserver { inner: self.inner.observer.clone() }
    }

    #[getter]
    fn covariance(&self) -> Option<Vec<Vec<f64>>> {
        match self.inner.covariance() {
            Some(cov) => Some(cov.data.as_vec().chunks(cov.ncols()).map(|x| x.to_vec()).collect()),
            None => None,
        }
    }

    fn set_covariance(&mut self, covariance: Vec<Vec<f64>>) {
        // ingest as a DMatrix
        let covariance = DMatrix::from_vec(covariance.len(), covariance[0].len(), covariance.iter().flatten().cloned().collect());
        self.inner.set_covariance(covariance);
    }

    fn set_mag(&mut self, mag: f64) {
        self.inner.set_mag(mag);
    }

    fn set_mag_err(&mut self, mag_err: f64) {
        self.inner.set_mag_err(mag_err);
    }


    

}