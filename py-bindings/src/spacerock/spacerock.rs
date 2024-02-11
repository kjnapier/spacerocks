use pyo3::prelude::*;
use pyo3::types::PyType;

use nalgebra::Vector3;

use spacerocks::spacerock::SpaceRock;
use spacerocks::time::Time;

use crate::time::time::PyTime;

use spice;

#[pyclass]
#[pyo3(name = "SpaceRock")]
pub struct PySpaceRock {
    pub inner: SpaceRock,
}

#[pymethods]
impl PySpaceRock {

    #[classmethod]
    fn from_horizons(cls: &PyType, name: &str) -> PyResult<Self> {
        let rock = SpaceRock::from_horizons(name);
        if rock.is_err() {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("Failed to create SpaceRock from Horizons for name: {}", name)));
        }
        Ok(PySpaceRock { inner: rock.unwrap() })
    }

    #[classmethod]
    fn from_spice(cls: &PyType, name: &str, epoch: &PyTime) -> PyResult<Self> {

        let ep = &epoch.inner;
        let rock = SpaceRock::from_spice(name, &ep);
        // if rock.is_err() {
        //     return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("Failed to create SpaceRock from Spice for name: {}", name)));
        // }

        Ok(PySpaceRock { inner: rock })
    }

    #[classmethod]
    fn random(cls: &PyType) -> Self {
        PySpaceRock { inner: SpaceRock::random() }
    }

    fn __repr__(&self) -> String {
        format!("SpaceRock: {}", self.inner.name)
    }

    fn analytic_propagate(&mut self, t: &PyTime) {
        // let ep = Time::new(t, "jd", "utc");
        self.inner.analytic_propagate(&t.inner);
    }

    fn change_frame(&mut self, frame: &str) {
        self.inner.change_frame(frame);
    }

    fn calculate_orbit(&mut self) {
        self.inner.calculate_orbit();
    }

    #[getter]
    fn epoch(&self) -> PyTime {
        PyTime { inner: self.inner.epoch.clone() }
    }

    #[getter]
    fn frame(&self) -> String {
        self.inner.frame.clone()
    }

    #[getter]
    fn r(&self) -> f64 {
        self.inner.r()
    }

    #[getter]
    fn e(&self) -> PyResult<f64> {
        match self.inner.orbit.as_ref() {
            Some(orbit) => Ok(orbit.e),
            None => Err(PyErr::new::<pyo3::exceptions::PyAttributeError, _>("Orbit not calculated")),
        }
    }

    #[getter]
    fn a(&self) -> PyResult<f64> {
        match self.inner.orbit.as_ref() {
            Some(orbit) => Ok(orbit.a),
            None => Err(PyErr::new::<pyo3::exceptions::PyAttributeError, _>("Orbit not calculated")),
        }
    }

    #[getter]
    fn inc(&self) -> PyResult<f64> {
        match self.inner.orbit.as_ref() {
            Some(orbit) => Ok(orbit.inc),
            None => Err(PyErr::new::<pyo3::exceptions::PyAttributeError, _>("Orbit not calculated")),
        }
    }

    #[getter]
    fn node(&self) -> PyResult<f64> {
        match self.inner.orbit.as_ref() {
            Some(orbit) => Ok(orbit.node),
            None => Err(PyErr::new::<pyo3::exceptions::PyAttributeError, _>("Orbit not calculated")),
        }
    }


    #[getter]
    fn name(&self) -> String {
        self.inner.name.clone()
    }

    #[getter]
    fn position(&self) -> (f64, f64, f64) {
        (self.inner.position.x, self.inner.position.y, self.inner.position.z)
    }

    #[getter]
    fn velocity(&self) -> (f64, f64, f64) {
        (self.inner.velocity.x, self.inner.velocity.y, self.inner.velocity.z)
    }

    #[getter]
    fn x(&self) -> f64 {
        self.inner.position.x
    }

    #[getter]
    fn y(&self) -> f64 {
        self.inner.position.y
    }

    #[getter]
    fn z(&self) -> f64 {
        self.inner.position.z
    }

    #[getter]
    fn vx(&self) -> f64 {
        self.inner.velocity.x
    }

    #[getter]
    fn vy(&self) -> f64 {
        self.inner.velocity.y
    }

    #[getter]
    fn vz(&self) -> f64 {
        self.inner.velocity.z
    }


}

