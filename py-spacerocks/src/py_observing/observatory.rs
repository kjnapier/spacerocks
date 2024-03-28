use pyo3::prelude::*;
use pyo3::types::PyType;
use pyo3::exceptions::PyValueError;

use spacerocks::Observatory;

use crate::py_observing::observer::PyObserver;
use crate::py_time::time::PyTime;
use crate::py_spacerock::origin::PyOrigin;

#[pyclass]
#[pyo3(name = "Observatory")]
pub struct PyObservatory {
    pub inner: Observatory,
}

#[pymethods]
impl PyObservatory {


    #[classmethod]
    fn from_obscode(_cls: &PyType, obscode: &str) -> PyResult<Self> {
        match Observatory::from_obscode(obscode) {
            Ok(o) => Ok(PyObservatory { inner: o }),
            Err(e) => Err(PyValueError::new_err(e))
        }
    }

    #[classmethod]
    pub fn from_parallax(_cls: &PyType, lon: f64, rho_sin: f64, rho_cos: f64) -> Self {
        PyObservatory { inner: Observatory::from_parallax(lon, rho_sin, rho_cos) }
    }

   //  #[pyo3(signature = (epoch, origin="SSB"))]
    fn at(&self, epoch: &PyTime, origin: &PyOrigin) -> PyObserver {
        let ep = &epoch.inner;
        let or = &origin.inner;
        let observer = self.inner.at(&ep, &or);

        PyObserver { inner: observer }
    }

    #[getter]
    fn lat(&self) -> f64 {
        self.inner.lat()
    }

    #[getter]
    fn lon(&self) -> f64 {
        self.inner.lon
    }

    #[getter]
    fn rho(&self) -> f64 {
        self.inner.rho()
    }


}