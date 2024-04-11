use pyo3::prelude::*;

pub mod spacerock;
pub mod rockcollection;
pub mod origin;
pub mod coordinate_frame;

pub fn make_spacerock_submodule(py: Python, m: &PyModule) -> PyResult<()> {
    // Add the `spacerock` submodule
    let submodule = PyModule::new(py, "spacerock")?;

    submodule.add_class::<spacerock::PySpaceRock>()?;
    submodule.add_class::<rockcollection::RockCollection>()?;
    submodule.add_class::<origin::PyOrigin>()?;
    submodule.add_class::<coordinate_frame::PyCoordinateFrame>()?;


    m.add_submodule(submodule)?;
    py.import("sys")?
        .getattr("modules")?
        .set_item("spacerocks.spacerock", submodule)?;
    submodule.setattr("__name__", "spacerocks.spacerock")?;
    Ok(())
}