use pyo3::prelude::*;

use spacerocks::observing::Observer;

use crate::py_time::time::PyTime;

use numpy::{PyArray1, IntoPyArray};

#[pyclass]
#[pyo3(name = "Observer")]
#[derive(Clone)]
pub struct PyObserver {
    pub inner: Observer,
}

#[pymethods]
impl PyObserver {

    // #[classmethod]
    // pub fn from_spacerock(_cls: &PyType, rock: &PySpaceRock) -> Self {
    //     PyObserver { inner: Observer::from_spacerock(&rock.inner) }
    // }

    #[getter]
    fn position(&self, py: Python) -> Py<PyArray1<f64>> {
        let pos = vec![self.inner.position.x, self.inner.position.y, self.inner.position.z];
        pos.into_pyarray(py).to_owned().into()
    }

    #[getter]
    fn velocity(&self, py: Python) -> Py<PyArray1<f64>> {
        let vel = match self.inner.velocity {
            Some(vel) => vec![vel.x, vel.y, vel.z],
            None => vec![0.0, 0.0, 0.0],
        };
        vel.into_pyarray(py).to_owned().into()
    }

    #[getter]
    fn origin(&self) -> String {
        self.inner.origin.to_string()
    }

    // fn change_frame(&mut self, new_frame: &str) {
    //     let frame = CoordinateFrame::from_str(new_frame).unwrap();
    //     self.inner.change_frame(&frame);
    // }

    #[getter]
    fn lat(&self) -> Option<f64> {
        self.inner.observatory.lat()
    }

    #[getter]
    fn lon(&self) -> Option<f64> {
        self.inner.observatory.lon()
    }

    #[getter]
    fn rho(&self) -> Option<f64> {
        self.inner.observatory.rho()
    }

    #[getter]
    fn reference_plane(&self) -> String {
        self.inner.reference_plane.to_string()
    }

    #[getter]
    fn epoch(&self) -> PyTime {
        PyTime { inner: self.inner.epoch.clone() }
    }

    // display the observer
    fn __str__(&self) -> String {
        format!("Observer at epoch: {}", self.inner.epoch)
    }

    // display the observer
    fn __repr__(&self) -> String {
        format!("Observer at epoch: {}", self.inner.epoch)
    }
    
}