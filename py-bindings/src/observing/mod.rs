use pyo3::prelude::*;

pub mod observatory;
pub mod detection;
pub mod detectioncatalog;
pub mod observer;

pub fn make_observing_submodule(py: Python, m: &PyModule) -> PyResult<()> {
    let submodule = PyModule::new(py, "observing")?;

    submodule.add_class::<observatory::PyObservatory>()?;
    submodule.add_class::<detectioncatalog::DetectionCatalog>()?;
    submodule.add_class::<observer::PyObserver>()?;
    submodule.add_class::<detection::PyDetection>()?;

    m.add_submodule(submodule)?;
    py.import("sys")?
        .getattr("modules")?
        .set_item("spacerocks.observing", submodule)?;
    submodule.setattr("__name__", "spacerocks.observing")?;
    Ok(())
}