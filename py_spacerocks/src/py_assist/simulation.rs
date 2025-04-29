use pyo3::prelude::*;
use pyo3::types::PyType;

use spacerocks::assist::SpiceSimulation;
use spacerocks::Time;   

use crate::PySpaceRock;
// use crate::py_spacerock::rockcollection::RockCollection;
use crate::py_time::time::PyTime;
use crate::py_assist::force::PyForce;
use crate::py_coordinates::origin::PyOrigin;
use crate::PySpiceKernel;

use std::sync::Arc;

#[pyclass]
#[pyo3(name = "SpiceSimulation")]
pub struct PySpiceSimulation {
    pub inner: SpiceSimulation,
}

#[pymethods]
impl PySpiceSimulation {

    // #[new]
    // pub fn new() -> PyResult<Self> {
    //     match Simulation::new(&Time::now(), "J2000", "SSB") {
    //         Ok(sim) => Ok(PySimulation { inner: sim }),
    //         Err(e) => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))
    //     }
    // }

    /// Create a new simulation with the giant planets.
    ///
    /// # Arguments
    ///
    /// * `epoch` - The epoch of the simulation.
    /// * `reference_plane` - The reference plane of the simulation.
    /// * `origin` - The origin of the simulation.
    ///
    /// # Returns
    ///
    /// * `Result<Simulation, &'static str>` - The Simulation object.
    #[classmethod]
    pub fn giants(_cls: Py<PyType>, epoch: &PyTime, kernel: &PySpiceKernel) -> PyResult<Self> {
        match SpiceSimulation::giants(&epoch.inner, &kernel.inner) {
            Ok(sim) => Ok(PySpiceSimulation { inner: sim }),
            Err(e) => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))
        }
    }

    /// Create a new simulation with the JPL horizons perturbers.
    ///
    /// # Arguments
    ///
    /// * `epoch` - The epoch of the simulation.
    /// * `reference_plane` - The reference plane of the simulation.
    /// * `origin` - The origin of the simulation.
    ///
    /// # Returns
    ///
    /// * `Result<Simulation, &'static str>` - The Simulation object.
    #[classmethod]
    pub fn horizons(_cls: Py<PyType>, epoch: &PyTime, kernel: &PySpiceKernel) -> PyResult<Self> {
        // need to Arc the kernel to pass it to the simulation
        match SpiceSimulation::horizons(&epoch.inner, &kernel.inner) {
            
            Ok(sim) => Ok(PySpiceSimulation { inner: sim }),
            Err(e) => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))
        }
    }

    /// Add a SpaceRock to the simulation.
    ///
    /// # Arguments
    ///
    /// * `rock` - The SpaceRock to add.
    ///
    /// # Returns
    ///
    /// * `Result<(), &'static str>` - The result of the operation.
    pub fn add(&mut self, rock: &PySpaceRock) -> PyResult<()> {
        match self.inner.add(rock.inner.clone()) {
            Ok(_) => Ok(()),
            Err(e) => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))
        }
    }

    /// Integrate the simulation to a specific epoch.
    ///
    /// # Arguments
    ///
    /// * `epoch` - The epoch to integrate to.
    pub fn integrate(&mut self, epoch: &PyTime, kernel: &PySpiceKernel) {
        self.inner.integrate(&epoch.inner.clone(), &kernel.inner);
    }
    

    /// Step the simulation by one timestep.
    pub fn step(&mut self, kernel: &PySpiceKernel) {
        self.inner.step(&kernel.inner);
    }

    /// Add a force to the simulation.
    ///
    /// # Arguments
    ///
    /// * `force` - The force to add.
    pub fn add_force(&mut self, force: PyRef<PyForce>) {
        self.inner.add_force(force.inner.clone());
    }

    // /// Get a single particle from the simulation by name.
    // pub fn get_particle(&self, name: &str) -> PyResult<PySpaceRock> {
    //     let rock = self.inner.get_particle(name);
    //     match rock {
    //         Ok(r) => Ok(PySpaceRock { inner: r.clone() }),
    //         Err(e) => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))
    //     }
    // }

    #[getter]
    pub fn timestep(&self) -> f64 {
        self.inner.integrator.timestep()
    }

    #[getter]
    pub fn particles(&self) -> Vec<PySpaceRock> {
        self.inner.state.particles.iter().map(|rock| PySpaceRock { inner: rock.clone() }).collect()
    }

}