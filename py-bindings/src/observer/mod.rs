use pyo3::prelude::*;

pub mod observatory;

pub fn make_observer_submodule(py: Python, m: &PyModule) -> PyResult<()> {
    // Add the `spice` submodule
    let submodule = PyModule::new(py, "observer")?;

    submodule.add_class::<observatory::PyObservatory>()?;

    m.add_submodule(submodule)?;
    py.import("sys")?
        .getattr("modules")?
        .set_item("spacerocks.observer", submodule)?;
    submodule.setattr("__name__", "spacerocks.observer")?;
    Ok(())
}