use pyo3::prelude::*;

pub mod simulation;
pub mod integrator;
pub mod force;

use crate::py_nbody::integrator::PyIntegrator;
use crate::py_nbody::force::PyForce;


pub fn make_assist_submodule(py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    let submodule = PyModule::new(py, "assist")?;

    submodule.add_class::<simulation::PySpiceSimulation>()?;
    submodule.add_class::<PyIntegrator>()?;
    submodule.add_class::<PyForce>()?;

    m.add_submodule(&submodule)?;
    py.import("sys")?
        .getattr("modules")?
        .set_item("spacerocks.assist", submodule.clone())?;
    submodule.setattr("__name__", "spacerocks.assist")?;
    Ok(())
}