use pyo3::prelude::*;

pub mod gauss;
pub mod fit;
pub mod weights;
pub mod observers;

pub fn make_orbfit_submodule(py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    let submodule = PyModule::new(m.py(), "orbfit")?;

    submodule.add_function(wrap_pyfunction!(gauss::gauss_py, submodule.clone())?)?;
    submodule.add_function(wrap_pyfunction!(fit::fit, submodule.clone())?)?;
    submodule.add_function(wrap_pyfunction!(fit::fit_many, submodule.clone())?)?;
    submodule.add_function(wrap_pyfunction!(fit::residuals, submodule.clone())?)?;
    submodule.add_function(wrap_pyfunction!(fit::bk_iod, submodule.clone())?)?;
    submodule.add_function(wrap_pyfunction!(fit::herget_iod, submodule.clone())?)?;
    submodule.add_function(wrap_pyfunction!(fit::sequential_update, submodule.clone())?)?;
    submodule.add_function(wrap_pyfunction!(fit::predict, submodule.clone())?)?;
    submodule.add_function(wrap_pyfunction!(fit::comet_orbits, submodule.clone())?)?;
    submodule.add_function(wrap_pyfunction!(weights::veres_sigma, submodule.clone())?)?;
    submodule.add_function(wrap_pyfunction!(weights::debias, submodule.clone())?)?;
    submodule.add_function(wrap_pyfunction!(weights::occultation_radec, submodule.clone())?)?;
    submodule.add_function(wrap_pyfunction!(observers::observers, submodule.clone())?)?;
    submodule.add_class::<fit::PyOrbitFit>()?;

    m.add_submodule(&submodule)?;
    py.import("sys")?
        .getattr("modules")?
        .set_item("spacerocks.orbfit", submodule.clone())?;
    submodule.setattr("__name__", "spacerocks.orbfit")?;

    Ok(())
}
