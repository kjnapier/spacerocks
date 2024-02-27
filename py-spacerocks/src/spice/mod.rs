use pyo3::prelude::*;

pub mod spicekernel;

pub fn make_spice_submodule(py: Python, m: &PyModule) -> PyResult<()> {
    // Add the `spice` submodule
    let submodule = PyModule::new(py, "spice")?;

    submodule.add_class::<spicekernel::PySpiceKernel>()?;

    m.add_submodule(submodule)?;
    py.import("sys")?
        .getattr("modules")?
        .set_item("spacerocks.spice", submodule)?;
    submodule.setattr("__name__", "spacerocks.spice")?;
    Ok(())
}