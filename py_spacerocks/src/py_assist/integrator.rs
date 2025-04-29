use pyo3::prelude::*;
use pyo3::types::PyType;

use spacerocks::assist::integrators::{Integrator, IAS15};

#[pyclass]
#[pyo3(name = "Integrator")]
// #[derive(Clone)]
pub struct PyIntegrator {
    pub inner: Box<dyn Integrator + Send + Sync>,
}

#[pymethods]
impl PyIntegrator {

    #[classmethod]
    pub fn ias15(_cls: Py<PyType>, timestep: f64) -> PyResult<Self> {
        Ok(PyIntegrator { inner: Box::new(IAS15::new(timestep)) })
    }

    #[getter]
    pub fn timestep(&self) -> f64 {
        self.inner.timestep()
    }

    pub fn set_timestep(&mut self, timestep: f64) {
        self.inner.set_timestep(timestep);
    }

}