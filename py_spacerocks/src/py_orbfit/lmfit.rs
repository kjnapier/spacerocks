use pyo3::prelude::*;

use crate::py_observing::observation::PyObservation;
use crate::PySpaceRock;
use crate::py_nbody::simulation::PySimulation;

use spacerocks::orbfit::fitter::{fit_orbit_lm, residuals_and_derivatives, orbit_chisq};

// initial_guess: &[f64; 7], sim: Simulation

#[pyfunction]
#[pyo3(name = "fit_orbit_lm")]
pub fn fit_orbit_lm_py(py: Python<'_>, detections: Vec<PyRef<PyObservation>>, initial_guess: PyRef<PySpaceRock>, sim: PyRef<PySimulation>) -> PyResult<()> {
    let detections = detections.iter().map(|obs| &obs.inner).collect();

    let initial_guess = [
        initial_guess.inner.position.x,
        initial_guess.inner.position.y,
        initial_guess.inner.position.z,
        initial_guess.inner.velocity.x,
        initial_guess.inner.velocity.y,
        initial_guess.inner.velocity.z,
        initial_guess.inner.epoch.tdb().jd()
    ];

    let fit_result = fit_orbit_lm(&detections, &initial_guess, sim.inner.clone()).unwrap();
    // println!("{:?}", fit_result);
    Ok(())
}

#[pyfunction]
#[pyo3(name = "orbit_chisq")]
pub fn orbit_chisq_py(py: Python<'_>, detections: Vec<PyRef<PyObservation>>, initial_guess: PyRef<PySpaceRock>, sim: PyRef<PySimulation>) -> PyResult<f64> {
    let detections = detections.iter().map(|obs| &obs.inner).collect();

    let initial_guess = [
        initial_guess.inner.position.x,
        initial_guess.inner.position.y,
        initial_guess.inner.position.z,
        initial_guess.inner.velocity.x,
        initial_guess.inner.velocity.y,
        initial_guess.inner.velocity.z,
        initial_guess.inner.epoch.tdb().jd()
    ];

    let csq = orbit_chisq(&detections, &initial_guess, sim.inner.clone());
    
    Ok(csq)
}

#[pyfunction]
#[pyo3(name = "randj")]
pub fn randj_py(py: Python<'_>, detections: Vec<PyRef<PyObservation>>, initial_guess: PyRef<PySpaceRock>, sim: PyRef<PySimulation>) -> PyResult<()> {
    let detections = detections.iter().map(|obs| &obs.inner).collect();

    let initial_guess = [
        initial_guess.inner.position.x,
        initial_guess.inner.position.y,
        initial_guess.inner.position.z,
        initial_guess.inner.velocity.x,
        initial_guess.inner.velocity.y,
        initial_guess.inner.velocity.z,
        initial_guess.inner.epoch.tdb().jd()
    ];

    let fit_result = residuals_and_derivatives(&detections, &initial_guess, sim.inner.clone());
    println!("{:?}", fit_result);
    Ok(())
}