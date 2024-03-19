use pyo3::prelude::*;
use pyo3::types::PyType;


// get Arc
use std::sync::Arc;
use std::rc::Rc; 

use spacerocks::nbody::integrators::{Integrator, Leapfrog, IAS15};

#[pyclass]
#[pyo3(name = "Integrator")]
// #[derive(Clone)]
pub struct PyIntegrator {
    pub inner: Box<dyn Integrator + Send + Sync>,
}

#[pymethods]
impl PyIntegrator {

    #[classmethod]
    pub fn leapfrog(_cls: &PyType, timestep: f64) -> PyResult<Self> {
        Ok(PyIntegrator { inner: Box::new(Leapfrog::new(timestep)) })
    }

    #[classmethod]
    pub fn ias15(_cls: &PyType, timestep: f64) -> PyResult<Self> {
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