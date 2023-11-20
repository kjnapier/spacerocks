use pyo3::prelude::*;


#[pyfunction]
#[pyo3(name = "calc_E_from_M")]
pub fn calc_E_from_M_py(e: f64, M: f64) -> PyResult<f64> {
    Ok(0.0)
}

#[pymodule]
pub fn spacerocks(py: Python, m: &PyModule) -> PyResult<()> {

    // Add the `transforms` submodule
    let submodule = PyModule::new(py, "transforms")?;
    submodule.add_function(wrap_pyfunction!(calc_E_from_M_py, submodule)?)?;
    m.add_submodule(submodule)?;
    py.import("sys")?
        .getattr("modules")?
        .set_item("spacerocks.transforms", submodule)?;
    submodule.setattr("__name__", "spacerocks.transforms")?;

    // add the `orbfit` submodule
    let submodule = PyModule::new(py, "orbfit")?;
    m.add_submodule(submodule)?;
    py.import("sys")?
        .getattr("modules")?
        .set_item("spacerocks.orbfit", submodule)?;
    submodule.setattr("__name__", "spacerocks.orbfit")?;


    // add the `time` submodule
    let submodule = PyModule::new(py, "time")?;
    m.add_submodule(submodule)?;
    py.import("sys")?
        .getattr("modules")?
        .set_item("spacerocks.time", submodule)?;
    submodule.setattr("__name__", "spacerocks.time")?;

    Ok(())
}