use pyo3::prelude::*;
use pyo3::types::PyType;

use nalgebra::Vector3;
use spacerocks::observing::Observer;

use numpy::{PyArray1, IntoPyArray};

#[pyclass]
#[pyo3(name = "Observer")]
pub struct PyObserver {
    pub inner: Observer,
}

#[pymethods]
impl PyObserver {

    // fn __repr__(&self) -> String {
    //     format!("Observer: {} at position: {:?}", self.inner.name, self.inner.position)
    // }

    #[getter]
    fn position(&self, py: Python) -> Py<PyArray1<f64>> {
        // self.inner.position.clone().to_vec().into_pyarray(py).to_owned()
        let pos = vec![self.inner.position.x, self.inner.position.y, self.inner.position.z];
        pos.into_pyarray(py).to_owned()
    }

    #[getter]
    fn velocity(&self, py: Python) -> Py<PyArray1<f64>> {
        // self.inner.velocity.clone().to_vec().into_pyarray(py).to_owned()
        let vel = vec![self.inner.velocity.x, self.inner.velocity.y, self.inner.velocity.z];
        vel.into_pyarray(py).to_owned()
    }

    #[getter]
    fn lat(&self) -> Option<f64> {
        self.inner.lat
    }

    #[getter]
    fn lon(&self) -> Option<f64> {
        self.inner.lon
    }

    #[getter]
    fn elevation(&self) -> Option<f64> {
        self.inner.elevation
    }
    
}