use pyo3::prelude::*;

use crate::py_observing::observation::PyObservation;
use crate::PySpaceRock;

use spacerocks::orbfit::gauss;

#[pyfunction]
#[pyo3(name = "gauss")]
pub fn gauss_py(_py: Python, o1: &PyObservation, o2: &PyObservation, o3: &PyObservation, min_distance: f64) -> PyResult<Option<Vec<PySpaceRock>>> {
    let rocks = gauss(&o1.inner, &o2.inner, &o3.inner, min_distance);
    match rocks {
        Some(rocks) => {
            let py_rocks: Vec<PySpaceRock> = rocks.iter().map(|rock| PySpaceRock { inner: rock.clone() }).collect();
            Ok(Some(py_rocks))
        },
        None => {
            // Err(PyErr::new::<pyo3::exceptions::PyValueError, _>("No solutions found"))
            Ok(None)
        }
    }
}

// #[pyfunction]
// #[pyo3(name = "gauss2")]
// pub fn gauss2_py(py: Python<'_>, input: Vec<PyRef<PyObservation>>) -> PyResult<()> {
//     for obs in input {
//         println!("{:?}", obs.inner);
//     }
//     Ok(())
// }
