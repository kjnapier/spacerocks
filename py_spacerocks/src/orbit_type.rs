use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;
use pyo3::types::PyType;
use spacerocks::OrbitType;

#[pyclass]
#[derive(Clone)]
#[pyo3(name = "OrbitType")]
pub struct PyOrbitType {
    pub inner: OrbitType,
}

#[pymethods]
#[allow(non_snake_case)]
impl PyOrbitType {
    #[classattr]
    fn HYPERBOLIC() -> Self {
        Self { inner: OrbitType::Hyperbolic }
    }

    #[classattr]
    fn PARABOLIC() -> Self {
        Self { inner: OrbitType::Parabolic }
    }

    #[classattr]
    fn ELLIPTICAL() -> Self {
        Self { inner: OrbitType::Elliptical }
    }

    #[classattr]
    fn CIRCULAR() -> Self {
        Self { inner: OrbitType::Circular }
    }

    #[classattr]
    fn RADIAL() -> Self {
        Self { inner: OrbitType::Radial }
    }

    #[classmethod]
    fn from_eccentricity(_cls: &Bound<'_, PyType>, e: f64, threshold: f64) -> PyResult<PyOrbitType> {
        match OrbitType::from_eccentricity(e, threshold) {
            Ok(orbit_type) => Ok(PyOrbitType { inner: orbit_type }),
            Err(e) => Err(PyValueError::new_err(format!("Failed to create OrbitType: {:?}", e)))
        }
    }

    fn __eq__(&self, other: &Self) -> bool {
        self.inner == other.inner
    }


    fn __str__(&self) -> PyResult<String> {
        Ok(format!("{:?}", self.inner))
    }

    fn __repr__(&self) -> PyResult<String> {
        Ok(format!("OrbitType::{:?}", self.inner))
    }
}