use pyo3::prelude::*;

use crate::py_observing::observation::PyObservation;
use crate::PySpaceRock;
use crate::py_nbody::simulation::PySimulation;

use spacerocks::orbfit::fitter::{fit_orbit_lm, residuals_and_derivatives, orbit_chisq, FitResult};
use pyo3::exceptions::PyValueError;
use pyo3::exceptions::PyRuntimeError;

// initial_guess: &[f64; 7], sim: Simulation

#[pyclass]
#[derive(Clone)]
pub struct PyFitResult {
    pub inner: FitResult,
}

#[pymethods]
// #[pyo3(name = "FitResult")]
impl PyFitResult {
    #[getter]
    fn chisq(&self) -> f64 {
        self.inner.chisq
    }

    #[getter]
    fn dof(&self) -> f64 {
        self.inner.dof
    }

    #[getter]
    fn niter(&self) -> usize {
        self.inner.niter
    }

    #[getter]
    fn residuals(&self) -> Vec<f64> {
        self.inner.residuals.clone()
    }

    #[getter]
    fn ra_residuals(&self) -> Vec<f64> {
        self.inner.ra_residuals.clone()
    }

    #[getter]
    fn dec_residuals(&self) -> Vec<f64> {
        self.inner.dec_residuals.clone()
    }

    #[getter]
    fn rock(&self) -> PyResult<Py<PySpaceRock>> {
        Python::with_gil(|py| {
            let rock = PySpaceRock { inner: self.inner.rock.clone() };
            Py::new(py, rock)  // This ensures correct type conversion
        })
    }

    #[getter]
    fn covariance(&self) -> Vec<f64> {
        let rows = self.inner.covariance.nrows();
        let cols = self.inner.covariance.ncols();
        let mut result = Vec::with_capacity(rows * cols);
        for i in 0..rows {
            for j in 0..cols {
                result.push(self.inner.covariance[(i, j)]);
            }
        }
        result
    }

    #[getter]
    fn keplerian_covariance(&self) -> Vec<f64> {
        let rows = self.inner.keplerian_covariance.nrows();
        let cols = self.inner.keplerian_covariance.ncols();
        let mut result = Vec::with_capacity(rows * cols);
        for i in 0..rows {
            for j in 0..cols {
                result.push(self.inner.keplerian_covariance[(i, j)]);
            }
        }
        result
    }
}

#[pyfunction]
#[pyo3(name = "fit_orbit_lm")]
pub fn fit_orbit_lm_py(py: Python<'_>, detections: Vec<PyRef<PyObservation>>, initial_guess: PyRef<PySpaceRock>, sim: PyRef<PySimulation>) -> PyResult<PyFitResult> {
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

    match fit_orbit_lm(&detections, &initial_guess, sim.inner.clone()) {
        Ok(Some(result)) => Ok(PyFitResult { inner: result }),
        Ok(None) => Err(PyErr::new::<PyValueError, _>("No valid fit found")),
        Err(e) => Err(PyErr::new::<PyRuntimeError, _>(e.to_string()))
    }
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