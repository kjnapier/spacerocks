use pyo3::prelude::*;

pub mod calc_E_from_M;
pub use self::calc_E_from_M::calc_E_from_M_py;

pub mod calc_E_from_f;
pub use self::calc_E_from_f::calc_E_from_f_py;

pub mod calc_M_from_E;
pub use self::calc_M_from_E::calc_M_from_E_py;

pub mod calc_f_from_E;
pub use self::calc_f_from_E::calc_f_from_E_py;

pub fn make_transforms_submodule(py: Python, m: &PyModule) -> PyResult<()> {
    // Add the `transforms` submodule
    let submodule = PyModule::new(py, "transforms")?;

    submodule.add_function(wrap_pyfunction!(calc_E_from_M_py, submodule)?)?;
    submodule.add_function(wrap_pyfunction!(calc_E_from_f_py, submodule)?)?;
    submodule.add_function(wrap_pyfunction!(calc_M_from_E_py, submodule)?)?;
    submodule.add_function(wrap_pyfunction!(calc_f_from_E_py, submodule)?)?;

    m.add_submodule(submodule)?;
    py.import("sys")?
        .getattr("modules")?
        .set_item("spacerocks.transforms", submodule)?;
    submodule.setattr("__name__", "spacerocks.transforms")?;
    Ok(())
}