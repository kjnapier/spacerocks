use pyo3::prelude::*;
use pyo3::types::PyType;

use spacerocks::CoordinateFrame;

#[pyclass]
#[pyo3(name = "CoordinateFrame")]
#[derive(Clone, Debug, PartialEq)]
pub struct PyCoordinateFrame {
    pub inner: CoordinateFrame,
}

#[pymethods]
impl PyCoordinateFrame {
    #[new]
    fn new() -> Self {
        PyCoordinateFrame {
            inner: CoordinateFrame::default(),
        }
    }

    #[classmethod]
    fn j2000(_cls: &PyType) -> Self {
        PyCoordinateFrame {
            inner: CoordinateFrame::J2000,
        }
    }

    #[classmethod]
    fn eclipj2000(_cls: &PyType) -> Self {
        PyCoordinateFrame {
            inner: CoordinateFrame::ECLIPJ2000,
        }
    }

    #[classmethod]
    fn invariable(_cls: &PyType) -> Self {
        PyCoordinateFrame {
            inner: CoordinateFrame::INVARIABLE,
        }
    }

    #[classmethod]
    fn galactic(_cls: &PyType) -> Self {
        PyCoordinateFrame {
            inner: CoordinateFrame::GALACTIC,
        }
    }

    #[classmethod]
    fn fk4(_cls: &PyType) -> Self {
        PyCoordinateFrame {
            inner: CoordinateFrame::FK4,
        }
    }

}