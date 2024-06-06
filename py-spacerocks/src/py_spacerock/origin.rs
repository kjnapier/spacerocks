use pyo3::prelude::*;
use pyo3::types::PyType;

use spacerocks::Origin;

#[pyclass]
#[pyo3(name = "Origin")]
#[derive(Clone, Debug, PartialEq)]
pub struct PyOrigin {
    pub inner: Origin,
}

#[pymethods]
impl PyOrigin {

    #[new]
    pub fn new(name: &str, mu: f64) -> Self {
        PyOrigin { inner: Origin::new_custom(mu, name) }
    }

    #[classmethod]
    pub fn sun(_cls: &PyType) -> Self {
        PyOrigin { inner: Origin::SUN }
    }

    #[classmethod]
    pub fn ssb(_cls: &PyType) -> Self {
        PyOrigin { inner: Origin::SSB }
    }

    pub fn __repr__(&self) -> String {
        format!("Origin: {} with mu = {}", self.inner.name(), self.inner.mu())
    }

    #[getter]
    pub fn mu(&self) -> f64 {
        self.inner.mu()
    }

    #[getter]
    pub fn name(&self) -> String {
        self.inner.name().to_string()
    }

}
