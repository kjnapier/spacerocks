use pyo3::prelude::*;

pub mod simulation;
pub mod integrator;
use crate::nbody::integrator::PyIntegrator;


pub fn make_nbody_submodule(py: Python, m: &PyModule) -> PyResult<()> {
    let submodule = PyModule::new(py, "nbody")?;

    submodule.add_class::<simulation::PySimulation>()?;
    submodule.add_class::<PyIntegrator>()?;

    m.add_submodule(submodule)?;
    py.import("sys")?
        .getattr("modules")?
        .set_item("spacerocks.nbody", submodule)?;
    submodule.setattr("__name__", "spacerocks.nbody")?;
    Ok(())
}