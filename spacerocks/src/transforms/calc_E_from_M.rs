use pyo3::prelude::*;

use spacerocks::transforms::calc_E_from_M;

#[pyfunction]
#[pyo3(name = "calc_E_from_M")]
pub fn calc_E_from_M_py(e: f64, M: f64) -> PyResult<f64> {
    Ok(calc_E_from_M(e, M))
}