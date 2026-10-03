use pyo3::prelude::*;
use pyo3::types::PyType;

use spacerocks::nbody::integrators::{Integrator, Leapfrog, IAS15, Trace, WisdomHolman};

#[pyclass]
#[pyo3(name = "Integrator")]
// #[derive(Clone)]
pub struct PyIntegrator {
    pub inner: Box<dyn Integrator + Send + Sync>,
}

#[pymethods]
impl PyIntegrator {

    #[classmethod]
    pub fn leapfrog(_cls: Py<PyType>, timestep: f64) -> PyResult<Self> {
        Ok(PyIntegrator { inner: Box::new(Leapfrog::new(timestep)) })
    }

    #[classmethod]
    pub fn wisdom_holman(_cls: Py<PyType>, timestep: f64) -> PyResult<Self> {
        Ok(PyIntegrator { inner: Box::new(WisdomHolman::new(timestep)) })
    }

    #[classmethod]
    pub fn trace(_cls: Py<PyType>, timestep: f64) -> PyResult<Self> {
        Ok(PyIntegrator { inner: Box::new(Trace::new(timestep)) })
    }

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