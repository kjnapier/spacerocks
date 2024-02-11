use pyo3::prelude::*;

pub mod time;

pub fn make_time_submodule(py: Python, m: &PyModule) -> PyResult<()> {
    // Add the `spice` submodule
    let submodule = PyModule::new(py, "time")?;

    submodule.add_class::<time::PyTime>()?;

    m.add_submodule(submodule)?;
    py.import("sys")?
        .getattr("modules")?
        .set_item("spacerocks.time", submodule)?;
    submodule.setattr("__name__", "spacerocks.time")?;
    Ok(())
}