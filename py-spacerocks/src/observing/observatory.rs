use pyo3::prelude::*;
use pyo3::types::PyType;

use spacerocks::Observatory;

use crate::observing::observer::PyObserver;
use crate::time::time::PyTime;

#[pyclass]
#[pyo3(name = "Observatory")]
pub struct PyObservatory {
    pub inner: Observatory,
}

#[pymethods]
impl PyObservatory {

    // #[classmethod]
    // fn from_obscode(name: &str, obscode: &str) -> Self {
    //     PyObservatory { inner: Observatory::new(name, code) }
    // }

    #[classmethod]
    pub fn from_coordinates(_cls: &PyType, lat: f64, lon: f64, elevation: f64) -> Self {
        PyObservatory { inner: Observatory::from_coordinates(lat, lon, elevation) }
    }

    // #[classmethod]
    // pub fn from_parallax(_cls: &PyType, lon: f64, rho_sin: f64, rho_cos: f64, elevation: f64) -> Self {
    //     PyObservatory { inner: Observatory::from_parallax(lat, lon, elevation) }
    // }

    fn at(&self, epoch: &PyTime) -> PyObserver {
        let ep = &epoch.inner;
        let observer = self.inner.at(&ep);

        PyObserver { inner: observer }
    }

    #[getter]
    fn lat(&self) -> f64 {
        self.inner.lat
    }

    #[getter]
    fn lon(&self) -> f64 {
        self.inner.lon
    }

    #[getter]
    fn elevation(&self) -> f64 {
        self.inner.elevation
    }


}