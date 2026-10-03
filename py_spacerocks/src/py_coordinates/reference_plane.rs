use pyo3::prelude::*;
use pyo3::types::PyType;

use spacerocks::ReferencePlane;

#[pyclass(from_py_object)]
#[pyo3(name = "ReferencePlane")]
#[derive(Clone, Debug, PartialEq)]
pub struct PyReferencePlane {
    pub inner: ReferencePlane,
}

#[pymethods]
impl PyReferencePlane {
    #[new]
    fn new() -> Self {
        PyReferencePlane {
            inner: ReferencePlane::default(),
        }
    }

    #[classmethod]
    fn j2000(_cls: Py<PyType>) -> Self {
        PyReferencePlane {
            inner: ReferencePlane::J2000,
        }
    }

    #[classmethod]
    fn eclipj2000(_cls: Py<PyType>) -> Self {
        PyReferencePlane {
            inner: ReferencePlane::ECLIPJ2000,
        }
    }

    #[classmethod]
    fn invariable(_cls: Py<PyType>) -> Self {
        PyReferencePlane {
            inner: ReferencePlane::INVARIABLE,
        }
    }

    #[classmethod]
    fn galactic(_cls: Py<PyType>) -> Self {
        PyReferencePlane {
            inner: ReferencePlane::GALACTIC,
        }
    }

    #[classmethod]
    fn fk4(_cls: Py<PyType>) -> Self {
        PyReferencePlane {
            inner: ReferencePlane::FK4,
        }
    }

}