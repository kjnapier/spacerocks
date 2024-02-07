use pyo3::prelude::*;

use spacerocks::transforms::calc_E_from_f;

#[pyfunction]
#[pyo3(name = "calc_E_from_f")]
pub fn calc_E_from_f_py(e: f64, f: f64) -> PyResult<f64> {
    Ok(calc_E_from_f(e, f))
}