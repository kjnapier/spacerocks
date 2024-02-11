use pyo3::prelude::*;
use pyo3::types::PyType;

use spice;
use spacerocks::spice::SpiceKernel;

#[pyclass]
#[pyo3(name = "SpiceKernel")]
pub struct PySpiceKernel {
    pub inner: SpiceKernel,
}

#[pymethods]
impl PySpiceKernel {

    #[new]
    fn new() -> Self {
        PySpiceKernel { inner: SpiceKernel::new() }
    }

    fn load(&mut self, path: &str) {
        self.inner.load(path);
    }

    fn unload(&mut self) {
        self.inner.unload();
    }

    fn __repr__(&self) -> String {
        format!("SpiceKernel: {:?}", self.inner.loaded_files)
    }

    #[getter]
    fn loaded_files(&self) -> Vec<String> {
        self.inner.loaded_files.clone()
    }
    
}

