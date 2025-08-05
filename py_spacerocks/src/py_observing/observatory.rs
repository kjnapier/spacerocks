use pyo3::prelude::*;
use pyo3::types::PyType;
use pyo3::exceptions::PyValueError;

use spacerocks::Observatory;
use crate::py_spice::spicekernel::PySpiceKernel;

use crate::py_observing::observer::PyObserver;
use crate::py_time::time::PyTime;
// use crate::py_spacerock::origin::PyOrigin;

#[pyclass]
#[pyo3(name = "Observatory")]
pub struct PyObservatory {
    pub inner: Observatory,
}

#[pymethods]
impl PyObservatory {


    #[classmethod]
    fn from_obscode(_cls: Py<PyType>, obscode: &str) -> PyResult<PyObservatory> {
        match Observatory::from_obscode(obscode) {
            Ok(o) => Ok(PyObservatory { inner: o }),
            Err(e) => Err(PyValueError::new_err(e.to_string()))
        }
    }

    #[classmethod]
    fn from_name(_cls: Py<PyType>, name: &str) -> PyResult<PyObservatory> {
        match Observatory::from_name(name) {
            Ok(o) => Ok(PyObservatory { inner: o }),
            Err(e) => Err(PyValueError::new_err(e.to_string()))
        }
    }

    #[pyo3(signature = (epoch, kernel, reference_plane="J2000", origin="SSB"))]
    fn at(&self, epoch: &PyTime, kernel: &PySpiceKernel, reference_plane: &str, origin: &str) -> PyResult<PyObserver> {
        let ep = &epoch.inner;
        let ker = &kernel.inner;
        match self.inner.at(ep, reference_plane, origin, ker) {
            Ok(o) => Ok(PyObserver { inner: o }),
            Err(e) => Err(PyValueError::new_err(e.to_string()))
        }
    }

    #[getter]
    fn lat(&self) -> Option<f64> {
        self.inner.lat()
    }

    #[getter]
    fn lon(&self) -> Option<f64> {
        self.inner.lon()
    }

    #[getter]
    fn rho(&self) -> Option<f64> {
        self.inner.rho()
    }


}