use pyo3::prelude::*;
use pyo3::wrap_pymodule;
// use crate::transforms::__PYO3_PYMODULE_DEF_TRANSFORMS;

pub mod transforms;
// use crate::transforms::transforms;

// #[pyfunction]
// fn func() -> PyResult<()> {
//     println!("Hello, world!");
//     Ok(())
// }

use crate::transforms::calc_E_from_M_py;

#[pymodule]
pub fn spacerocks(py: Python, m: &PyModule) -> PyResult<()> {
    // m.add_wrapped(wrap_pymodule!(transforms))?;

    let child_module = PyModule::new(py, "transforms")?;
    child_module.add_function(wrap_pyfunction!(calc_E_from_M_py, child_module)?)?;
    m.add_submodule(child_module)?;
    // m.add_submodule(transforms)?;
    Ok(())
}


// #[pymodule]
// pub fn spacerocks(py: Python, m: &PyModule) -> PyResult<()> {
//     m.add_wrapped(wrap_pymodule!(transforms))?;
//     Ok(())
// }