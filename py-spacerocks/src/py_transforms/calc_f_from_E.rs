use pyo3::prelude::*;

use spacerocks::transforms::calc_f_from_E;

#[pyfunction]
#[pyo3(name = "calc_f_from_E")]
pub fn calc_f_from_E_py(e: f64, f: f64) -> PyResult<f64> {
    Ok(calc_f_from_E(e, f))
}