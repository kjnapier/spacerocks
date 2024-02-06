use pyo3::prelude::*;

pub mod calc_E_from_M;
pub use self::calc_E_from_M::calc_E_from_M;

#[pyfunction]
#[pyo3(name = "calc_E_from_M")]
pub fn calc_E_from_M_py(e: f64, M: f64) -> PyResult<f64> {
    Ok(calc_E_from_M(e, M))
}

pub mod calc_E_from_f;
pub use self::calc_E_from_f::calc_E_from_f;

#[pyfunction]
#[pyo3(name = "calc_E_from_f")]
pub fn calc_E_from_f_py(e: f64, f: f64) -> PyResult<f64> {
    Ok(calc_E_from_f(e, f))
}

pub mod calc_M_from_E;
pub use self::calc_M_from_E::calc_M_from_E;

#[pyfunction]
#[pyo3(name = "calc_M_from_E")]
pub fn calc_M_from_E_py(e: f64, E: f64) -> PyResult<f64> {
    Ok(calc_M_from_E(e, E))
}

pub fn make_transforms_submodule(py: Python, m: &PyModule) -> PyResult<()> {
    // Add the `transforms` submodule
    let submodule = PyModule::new(py, "transforms")?;

    submodule.add_function(wrap_pyfunction!(calc_E_from_M_py, submodule)?)?;
    submodule.add_function(wrap_pyfunction!(calc_E_from_f_py, submodule)?)?;
    submodule.add_function(wrap_pyfunction!(calc_M_from_E_py, submodule)?)?;

    m.add_submodule(submodule)?;
    py.import("sys")?
        .getattr("modules")?
        .set_item("spacerocks.transforms", submodule)?;
    submodule.setattr("__name__", "spacerocks.transforms")?;
    Ok(())
}