// src/transforms/mod.rs
use pyo3::prelude::*;
use pyo3::wrap_pyfunction;

pub mod calculate;  // Assuming the Rust function is in calculate.rs
pub use calculate::calc_E_from_M;  // The actual Rust function

#[pyfunction]
#[pyo3(name = "calc_E_from_M")]
pub fn calc_E_from_M_py(e: f64, M: f64) -> PyResult<f64> {
    Ok(calc_E_from_M(e, M))
}

// #[pymodule]
// fn transforms(py: Python, m: &PyModule) -> PyResult<()> {
//     m.add_function(wrap_pyfunction!(calc_E_from_M_py, m)?)?;
//     Ok(())
// }
