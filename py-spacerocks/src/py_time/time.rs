use pyo3::prelude::*;
use pyo3::types::PyType;

use spacerocks::time::Time;

#[pyclass]
#[pyo3(name = "Time")]
#[derive(Clone)]
pub struct PyTime {
    pub inner: Time,
}

#[pymethods]
impl PyTime {

    #[new]
    fn new(epoch: f64, timescale: &str, format: &str) -> Self {
        PyTime { inner: Time::new(epoch, timescale, format) }
    }

    #[classmethod]
    fn now(_cls: &PyType) -> PyResult<Self> {
        Ok(PyTime { inner: Time::now() })
    }

    fn to_utc(&mut self) {
        self.inner.to_utc();
    }

    fn to_tdb(&mut self) {
        self.inner.to_tdb();
    }

    fn __repr__(&self) -> String {
        format!("Time: {} {:?} {:?}", self.inner.epoch, self.inner.timescale, self.inner.format)
    }

    #[getter]
    fn epoch(&self) -> f64 {
        self.inner.epoch
    }

    #[getter]
    fn timescale(&self) -> &str {
        // &self.inner.timescale
        match self.inner.timescale {
            spacerocks::time::TimeScale::TDB => "TDB",
            spacerocks::time::TimeScale::UTC => "UTC",
        }
    }

    #[getter]
    fn format(&self) -> &str {
        // &self.inner.format
        match self.inner.format {
            spacerocks::time::TimeFormat::JD => "JD",
        }
    }

    fn calendar(&self) -> String {
        self.inner.calendar()
    }

    // define __add__ and __sub__ here
    fn __add__(&self, dt: f64) -> PyTime {
        PyTime { inner: &self.inner + dt }
    }

    fn __sub__(&self, dt: f64) -> PyTime {
        PyTime { inner: &self.inner - dt }
    }

    // fn __sub__(&self, other: &PyTime) -> f64 {
    //     self.inner.epoch - other.inner.epoch
    // }

    
}


