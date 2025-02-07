use pyo3::prelude::*;

pub mod gauss;
pub mod lmfit;

pub fn make_orbfit_submodule(py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    // Create a submodule named "orbfit"
    let submodule = PyModule::new(m.py(), "orbfit")?;
    
    // Add the submodule’s contents
    submodule.add_function(wrap_pyfunction!(gauss::gauss_py, submodule.clone())?)?;
    // submodule.add_function(wrap_pyfunction!(gauss::gauss2_py, submodule.clone())?)?;

    submodule.add_function(wrap_pyfunction!(lmfit::fit_orbit_lm_py, submodule.clone())?)?;
    submodule.add_function(wrap_pyfunction!(lmfit::randj_py, submodule.clone())?)?;



    // Add the submodule to the parent module (requires PyO3 ≥ 0.14)
    m.add_submodule(&submodule)?;

    // For a fully importable submodule in Python, register it in sys.modules
    py.import("sys")?
        .getattr("modules")?
        .set_item("spacerocks.orbfit", submodule.clone())?;

    // Update the submodule’s __name__ to reflect the fully qualified name
    submodule.setattr("__name__", "spacerocks.orbfit")?;

    Ok(())
}