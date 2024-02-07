use pyo3::prelude::*;

// pub mod calc_M_from_E;
// pub use self::calc_M_from_E::calc_M_from_E;

use spacerocks::transforms::calc_M_from_E;

#[pyfunction]
#[pyo3(name = "calc_M_from_E")]
pub fn calc_M_from_E_py(e: f64, E: f64) -> PyResult<f64> {
    Ok(calc_M_from_E(e, E))
}