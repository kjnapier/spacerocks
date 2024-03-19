use pyo3::prelude::*;

pub mod simulation;
pub mod integrator;
pub mod force;

use crate::py_nbody::integrator::PyIntegrator;
use crate::py_nbody::force::PyForce;


pub fn make_nbody_submodule(py: Python, m: &PyModule) -> PyResult<()> {
    let submodule = PyModule::new(py, "nbody")?;

    submodule.add_class::<simulation::PySimulation>()?;
    submodule.add_class::<PyIntegrator>()?;
    submodule.add_class::<PyForce>()?;

    m.add_submodule(submodule)?;
    py.import("sys")?
        .getattr("modules")?
        .set_item("spacerocks.nbody", submodule)?;
    submodule.setattr("__name__", "spacerocks.nbody")?;
    Ok(())
}