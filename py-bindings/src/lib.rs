use pyo3::prelude::*;

mod transforms;
use transforms::make_transforms_submodule;

// mod spacerock;
// use spacerock::make_spacerock_submodule;

// mod time;
// use time::make_time_submodule;



#[pymodule]
pub fn spacerocks(py: Python, m: &PyModule) -> PyResult<()> {

    // Add the `transforms` submodule
    make_transforms_submodule(py, m)?;

    // // Add the `spacerock` submodule
    // make_spacerock_submodule(py, m)?;

    // // Add the `time` submodule
    // make_time_submodule(py, m)?;

    Ok(())
}
