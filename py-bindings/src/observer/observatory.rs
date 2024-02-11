use pyo3::prelude::*;
use pyo3::types::PyType;

use spacerocks::Observatory;

use crate::spacerock::spacerock::PySpaceRock;
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

    fn at(&self, epoch: &PyTime) -> PySpaceRock {
        let ep = &epoch.inner;
        let rock = self.inner.at(&ep);
        PySpaceRock { inner: rock }
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