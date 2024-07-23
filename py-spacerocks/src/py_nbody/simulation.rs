use pyo3::prelude::*;
use pyo3::types::PyType;

use spacerocks::nbody::Simulation;
use spacerocks::spacerock::CoordinateFrame;

use crate::py_spacerock::spacerock::PySpaceRock;
use crate::py_spacerock::rockcollection::RockCollection;
use crate::py_time::time::PyTime;
use crate::py_nbody::integrator::PyIntegrator;
use crate::py_nbody::force::PyForce;
use crate::py_spacerock::origin::PyOrigin;

use std::sync::Arc;

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

    #[classmethod]
    pub fn giants(_cls: &PyType, epoch: &PyTime, frame: &str, origin: &PyOrigin) -> PyResult<Self> {
        let ep = &epoch.inner;
        let or = &origin.inner;

        let frame = CoordinateFrame::from_str(frame).unwrap();
        let sim = Simulation::giants(ep, &frame, or);
        Ok(PySimulation { inner: sim.unwrap() })
    }

    #[classmethod]
    pub fn planets(_cls: &PyType, epoch: &PyTime, frame: &str, origin: &PyOrigin) -> PyResult<Self> {
        let ep = &epoch.inner;
        let or = &origin.inner;

        let frame = CoordinateFrame::from_str(frame).unwrap();
        let sim = Simulation::planets(ep, &frame, or);
        Ok(PySimulation { inner: sim.unwrap() })
    }

    #[classmethod]
    pub fn horizons(_cls: &PyType, epoch: &PyTime, frame: &str, origin: &PyOrigin) -> PyResult<Self> {
        let ep = &epoch.inner;
        let or = &origin.inner;

        let frame = CoordinateFrame::from_str(frame).unwrap();
        let sim = Simulation::horizons(ep, &frame, or);
        Ok(PySimulation { inner: sim.unwrap() })
    }

    pub fn add(&mut self, rock: &PySpaceRock) -> PyResult<()> {
        match self.inner.add(rock.inner.clone()) {
            Ok(_) => Ok(()),
            Err(e) => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))
        }
    }

    pub fn remove(&mut self, name: &str) -> PyResult<()> {
        match self.inner.remove(name) {
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
        let frame = CoordinateFrame::from_str(frame).unwrap();
        self.inner.frame = Some(frame);
    }

    pub fn set_origin(&mut self, origin: &PyOrigin) {
        self.inner.origin = Some(origin.inner.clone());
    }

    pub fn set_integrator(&mut self, integrator: PyRef<PyIntegrator>) {
        self.inner.integrator = integrator.inner.clone();
    }

    pub fn add_force(&mut self, force: PyRef<PyForce>) {
        self.inner.add_force(force.inner.clone());
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

    // #[getter]
    // pub fn perturbers(&self) -> RockCollection {
    //     RockCollection { rocks: self.inner.perturbers.clone(), name_hash_map: self.inner.perturber_index_map.clone() }
    // }

    #[getter]
    pub fn epoch(&self) -> PyTime {
        PyTime { inner: self.inner.epoch.clone() }
    }

    #[getter]
    pub fn frame(&self) -> Option<String> {
        let f = self.inner.frame.clone();
        match f {
            Some(frame) => Some(frame.to_string()),
            None => None
        }
    }


    #[getter]
    pub fn timestep(&self) -> f64 {
        self.inner.integrator.timestep()
    }

    #[getter]
    pub fn origin(&self) -> Option<PyOrigin> {
        let o = self.inner.origin.clone();
        match o {
            Some(origin) => Some(PyOrigin { inner: origin.clone() }),
            None => None
        }
    }

    // #[getter]
    // pub fn integrator(&self) -> PyIntegrator {
    //     PyIntegrator { inner: self.inner.integrator.clone() }
    // }

}