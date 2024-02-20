use pyo3::prelude::*;
use spacerocks::nbody::Simulation;

use spacerocks::spacerock::SpaceRock;
use crate::spacerock::spacerock::PySpaceRock;
use crate::spacerock::rockcollection::RockCollection;
use crate::time::time::PyTime;
use crate::nbody::integrator::PyIntegrator;

use std::collections::HashMap;

#[pyclass]
#[pyo3(name = "Simulation")]
pub struct PySimulation {
    pub inner: Simulation,
}

#[pymethods]
impl PySimulation {

    #[new]
    pub fn new() -> Self {
        PySimulation { inner: Simulation::new() }
    }

    pub fn add(&mut self, rock: &PySpaceRock) -> PyResult<()> {
        match self.inner.add(rock.inner.clone()) {
            Ok(_) => Ok(()),
            Err(e) => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))
        }
    }

    pub fn integrate(&mut self, epoch: &PyTime) {
        self.inner.integrate(&epoch.inner.clone());
    }

    pub fn step(&mut self) {
        self.inner.step();
    }

    pub fn threaded_step(&mut self) {
        self.inner.threaded_step();
    }

    pub fn move_to_center_of_mass(&mut self) {
        self.inner.move_to_center_of_mass();
    }

    pub fn change_origin(&mut self, origin: &str) -> PyResult<()> {
        match self.inner.change_origin(origin) {
            Ok(_) => Ok(()),
            Err(e) => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))
        }
    }

    pub fn set_epoch(&mut self, epoch: &PyTime) {
        self.inner.epoch = epoch.inner.clone();
    }

    pub fn set_frame(&mut self, frame: &str) {
        self.inner.frame = Some(frame.to_string());
    }

    pub fn set_origin(&mut self, origin: &str) {
        self.inner.origin = Some(origin.to_string());
    }

    pub fn set_integrator(&mut self, integrator: PyIntegrator) {
        self.inner.integrator = integrator.inner;
    }

    pub fn energy(&self) -> f64 {
        self.inner.energy()
    }

    pub fn get_particle(&self, name: &str) -> PyResult<PySpaceRock> {
        let rock = self.inner.get_particle(name);
        match rock {
            Ok(r) => Ok(PySpaceRock { inner: r.clone() }),
            Err(e) => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))
        }
       
    }


    #[getter]
    pub fn particles(&self) -> RockCollection {
        RockCollection { rocks: self.inner.particles.clone(), name_hash_map: self.inner.particle_index_map.clone() }
    }

    #[getter]
    pub fn perturbers(&self) -> RockCollection {
        RockCollection { rocks: self.inner.perturbers.clone(), name_hash_map: self.inner.perturber_index_map.clone() }
    }

    #[getter]
    pub fn epoch(&self) -> PyTime {
        PyTime { inner: self.inner.epoch.clone() }
    }

    #[getter]
    pub fn frame(&self) -> Option<String> {
        self.inner.frame.clone()
    }

    #[getter]
    pub fn mass(&self) -> f64 {
        self.inner.mass
    }

    #[getter]
    pub fn timestep(&self) -> f64 {
        self.inner.integrator.timestep()
    }

    #[getter]
    pub fn origin(&self) -> Option<String> {
        self.inner.origin.clone()
    }

    #[getter]
    pub fn integrator(&self) -> PyIntegrator {
        PyIntegrator { inner: self.inner.integrator.clone() }
    }

}