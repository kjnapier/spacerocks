use pyo3::prelude::*;
use pyo3::types::PyType;

use spacerocks::Observation;

use crate::py_observing::observer::PyObserver;
use crate::py_time::time::PyTime;
// use crate::py_spacerock::origin::PyOrigin;

#[pyclass]
#[pyo3(name = "Observation")]
pub struct PyObservation {
    pub inner: Observation,
}

#[pymethods]
impl PyObservation {

    #[classmethod]
    #[pyo3(signature = (epoch, ra, dec, observer, mag = None))]
    fn from_astrometry(_cls: Py<PyType>, epoch: PyTime, ra: f64, dec: f64, observer: PyObserver, mag: Option<f64>) -> PyResult<PyObservation> {
        Ok(PyObservation { inner: Observation::from_astrometry(epoch.inner, ra, dec, mag, observer.inner) })
    }

    #[classmethod]
    #[pyo3(signature = (epoch, ra, dec, ra_rate, dec_rate, observer, mag = None))]
    fn from_streak(_cls: Py<PyType>, epoch: PyTime, ra: f64, dec: f64, ra_rate: f64, dec_rate: f64, observer: PyObserver, mag: Option<f64>) -> PyResult<PyObservation> {
        Ok(PyObservation { inner: Observation::from_streak(epoch.inner, ra, dec, ra_rate, dec_rate, mag, observer.inner) })
    }

    // repr
    fn __repr__(&self) -> PyResult<String> {
        // Ok(format!("Observation: {:?}", self.inner))
        // println!("Observation:");
        // println!("  ra: {:?}", self.inner.ra());
        // println!("  dec: {:?}", self.inner.dec());
        // println!("  ra_rate: {:?}", self.inner.ra_rate());
        // println!("  dec_rate: {:?}", self.inner.dec_rate());
        // println!("  range: {:?}", self.inner.range());
        // println!("  range_rate: {:?}", self.inner.range_rate());
        // println!("  epoch: {:?}", self.inner.epoch);
        // println!("  observer: {:?}", self.inner.observer);
        // Ok("".to_string())

        // Format positional components with 4 decimal places
        let position = format!("(RA={:.4}°, Dec={:.4}°)", 
        self.inner.ra().to_degrees(), self.inner.dec().to_degrees());

        // Format motion components if they exist
        let motion = match (self.inner.ra_rate(), self.inner.dec_rate()) {
            (Some(ra_rate), Some(dec_rate)) => 
                format!("ra/dec rates = ({:.4}°/day, {:.4}°/day)", ra_rate.to_degrees(), dec_rate.to_degrees()),
            _ => String::new()
        };

        // Format range components if they exist
        let range = match (self.inner.range(), self.inner.range_rate()) {
            (Some(r), Some(r_rate)) => 
                format!("distance = {:.2}AU (rate = {:.4}AU/day)", r, r_rate),
            _ => String::new()
        };

        // Format magnitude if it exists
        let mag = self.inner.mag()
            .map_or(String::new(), |m| format!("mag = {:.1}", m));

        Ok(format!("Observation:\n{}\n{}\n{}\n{}\nepoch: {}", 
            position,
            motion,
            range,
            mag, 
            self.inner.epoch.iso()))

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


    

}