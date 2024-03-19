use pyo3::prelude::*;

mod py_transforms;
use py_transforms::make_transforms_submodule;

mod py_spacerock;
use py_spacerock::make_spacerock_submodule;

mod py_spice;
use py_spice::make_spice_submodule;

mod py_time;
use py_time::make_time_submodule;

mod py_observing;
use py_observing::make_observing_submodule;

mod py_nbody;
use py_nbody::make_nbody_submodule;

// mod orbfit;
// use orbfit::make_orbfit_submodule;

#[pymodule]
pub fn spacerocks(py: Python, m: &PyModule) -> PyResult<()> {

    // Add the `transforms` submodule
    make_transforms_submodule(py, m)?;

    // Add the `spacerock` submodule
    make_spacerock_submodule(py, m)?;

    // Add the `spice` submodule
    make_spice_submodule(py, m)?;

    // // Add the `time` submodule
    make_time_submodule(py, m)?;

    // Add the `observing` submodule
    make_observing_submodule(py, m)?;

    // Add the `nbody` submodule
    make_nbody_submodule(py, m)?;

    // // Add the `orbfit` submodule
    // make_orbfit_submodule(py, m)?;

    Ok(())
}
