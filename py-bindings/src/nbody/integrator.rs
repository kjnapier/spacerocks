use pyo3::prelude::*;
use pyo3::types::PyType;

use spacerocks::nbody::integrator::Integrator;
use spacerocks::nbody::freezer::Freezer;
use spacerocks::nbody::leapfrog::Leapfrog;
use spacerocks::SpaceRock;

#[pyclass]
#[pyo3(name = "Integrator")]
#[derive(Clone)]
pub struct PyIntegrator {
    pub inner: Integrator,
}

#[pymethods]
impl PyIntegrator {

    #[classmethod]
    pub fn leapfrog(cls: &PyType, timestep: f64) -> PyResult<Self> {
        Ok(PyIntegrator { inner: Integrator::Leapfrog(Leapfrog::new(timestep)) })
    }

    #[classmethod]
    pub fn freezer(cls: &PyType, timestep: f64) -> PyResult<Self> {
        Ok(PyIntegrator { inner: Integrator::Freezer(Freezer::new(timestep)) })
    }

    #[getter]
    pub fn timestep(&self) -> f64 {
        self.inner.timestep()
    }

    pub fn set_timestep(&mut self, timestep: f64) {
        self.inner.set_timestep(timestep);
    }

}