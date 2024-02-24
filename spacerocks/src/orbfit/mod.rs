use pyo3::prelude::*;

pub mod gauss;
pub use gauss::gauss_fit;

pub fn make_orbfit_submodule(py: Python, m: &PyModule) -> PyResult<()> {
    let submodule = PyModule::new(py, "orbfit")?;

    submodule.add_function(wrap_pyfunction!(gauss_fit, submodule)?)?;

    m.add_submodule(submodule)?;
    py.import("sys")?
        .getattr("modules")?
        .set_item("spacerocks.orbfit", submodule)?;
    submodule.setattr("__name__", "spacerocks.orbfit")?;
    Ok(())
}
