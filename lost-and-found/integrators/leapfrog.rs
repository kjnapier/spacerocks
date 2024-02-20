use pyo3::prelude::*;

use spacerocks::SpaceRock;
use spacerocks::nbody::leapfrog::Leapfrog;

#[pyclass]
#[pyo3(name = "Leapfrog")]
pub struct PyLeapfrog {
    pub inner: Leapfrog,
}

#[pymethods]
impl PyLeapfrog {

    #[new]
    pub fn new(timestep: f64) -> Self {
        PyLeapfrog { inner: Leapfrog::new(timestep) }
    }

    #[getter]
    pub fn timestep(&self) -> f64 {
        self.inner.timestep
    }

}