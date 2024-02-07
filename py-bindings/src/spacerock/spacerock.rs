use pyo3::prelude::*;
use pyo3::types::PyType;

use nalgebra::Vector3;

use spacerocks::spacerock::SpaceRock;

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

        // PySpaceRock { inner: SpaceRock::from_horizons(name) })
    }

    fn __repr__(&self) -> String {
        format!("SpaceRock: {}", self.inner.name)
    }

    // #[classmethod]
    // fn from_state(cls: &PyType, name: &str, position: (f64, f64, f64), velocity: (f64, f64, f64)) -> Self {
    //     PySpaceRock { inner: SpaceRock { name: name.to_string(), position: Vector3::new(position.0, position.1, position.2), velocity: Vector3::new(velocity.0, velocity.1, velocity.2) } }
    // }

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

