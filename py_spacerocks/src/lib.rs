use pyo3::prelude::*;
// use crate::mpc::MPCHandler;

mod py_transforms;
use py_transforms::make_transforms_submodule;

mod py_coordinates;
use py_coordinates::make_coordinates_submodule;

mod py_spice;  
use py_spice::{make_spice_submodule, spicekernel::PySpiceKernel};

mod py_time;
use py_time::make_time_submodule;

mod spacerock;
use spacerock::PySpaceRock;

mod rockcollection;
use rockcollection::RockCollection;

mod py_observing;
use py_observing::make_observing_submodule;

mod py_nbody;
use py_nbody::make_nbody_submodule;

mod py_orbfit;
use py_orbfit::make_orbfit_submodule;
// use py_orbfit::lmfit::PyFitResult;

pub mod py_assist;
use py_assist::make_assist_submodule;


mod mpc;
// use mpc::MPC;

mod orbit_type;  // declare the module
use orbit_type::PyOrbitType;  // import the type

#[pymodule]
#[pyo3(name = "spacerocks")]
pub fn py_spacerocks(py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {

    // Add the `transforms` submodule
    make_transforms_submodule(py, m)?;

    // Add the `spice` submodule
    make_spice_submodule(py, m)?;

    // Add the `time` submodule
    make_time_submodule(py, m)?;

    // Add the `coordinates` submodule
    make_coordinates_submodule(py, m)?;

    // Add the `observing` submodule
    make_observing_submodule(py, m)?;

    // Add the `nbody` submodule
    make_nbody_submodule(py, m)?;

    // Add the `orbfit` submodule
    make_orbfit_submodule(py, m)?;

    // Add the `assist` submodule
    make_assist_submodule(py, m)?;


    m.add_class::<PySpaceRock>()?;
    m.add_class::<RockCollection>()?;
    // m.add_class::<MPCHandler>()?;
    m.add_class::<PySpiceKernel>()?;
    m.add_class::<PyOrbitType>()?; 

    Ok(())
}
