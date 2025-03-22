use pyo3::prelude::*;
use pyo3::types::PyType;

use spacerocks::nbody::Simulation;
use spacerocks::coordinates::{ReferencePlane, Origin};
use spacerocks::Time;   

use crate::PySpaceRock;
// use crate::py_spacerock::rockcollection::RockCollection;
use crate::py_time::time::PyTime;
use crate::py_nbody::integrator::PyIntegrator;
use crate::py_nbody::force::PyForce;
use crate::py_coordinates::origin::PyOrigin;
use crate::PySpiceKernel;

#[pyclass]
#[pyo3(name = "Simulation")]
pub struct PySimulation {
    pub inner: Simulation,
}

#[pymethods]
impl PySimulation {

    #[new]
    pub fn new() -> PyResult<Self> {
        match Simulation::new(&Time::now(), "J2000", "SSB") {
            Ok(sim) => Ok(PySimulation { inner: sim }),
            Err(e) => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))
        }
    }

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
    pub fn giants(_cls: Py<PyType>, epoch: &PyTime, reference_plane: &str, origin: &str, kernel: &PySpiceKernel) -> PyResult<Self> {
        match Simulation::giants(&epoch.inner, reference_plane, origin, &kernel.inner) {
            Ok(sim) => Ok(PySimulation { inner: sim }),
            Err(e) => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))
        }
    }

    /// Create a new simulation with the planets.
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
    pub fn planets(_cls: Py<PyType>, epoch: &PyTime, reference_plane: &str, origin: &str, kernel: &PySpiceKernel) -> PyResult<Self> {
        match Simulation::planets(&epoch.inner, reference_plane, origin, &kernel.inner) {
            Ok(sim) => Ok(PySimulation { inner: sim }),
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
    pub fn horizons(_cls: Py<PyType>, epoch: &PyTime, reference_plane: &str, origin: &str, kernel: &PySpiceKernel) -> PyResult<Self> {
        match Simulation::horizons(&epoch.inner, reference_plane, origin, &kernel.inner) {
            Ok(sim) => Ok(PySimulation { inner: sim }),
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

    /// Remove a SpaceRock from the simulation.
    ///
    /// # Arguments
    ///
    /// * `name` - The name of the SpaceRock to remove.
    ///
    /// # Returns
    ///
    /// * `Result<(), &'static str>` - The result of the operation.
    pub fn remove(&mut self, name: &str) -> PyResult<()> {
        match self.inner.remove(name) {
            Ok(_) => Ok(()),
            Err(e) => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))
        }
    }

    /// Integrate the simulation to a specific epoch.
    ///
    /// # Arguments
    ///
    /// * `epoch` - The epoch to integrate to.
    pub fn integrate(&mut self, epoch: &PyTime) {
        self.inner.integrate(&epoch.inner.clone());
    }
    
    /// Step the simulation by one timestep.
    pub fn step(&mut self) {
        self.inner.step();
    }

    /// Move the simulation to the center of mass.
    pub fn move_to_center_of_mass(&mut self) -> PyResult<()> {
        match self.inner.move_to_center_of_mass() {
            Ok(_) => Ok(()),
            Err(e) => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))
        }
    }

    /// Change the origin of the simulation. The new origin must be present in the simulation.
    ///
    /// # Arguments
    ///
    /// * `origin` - The new origin of the simulation.
    ///
    /// # Returns
    ///
    /// * `Result<(), &'static str>` - The result of the operation.
    pub fn change_origin(&mut self, origin: &str) -> PyResult<()> {
        match self.inner.change_origin(origin) {
            Ok(_) => Ok(()),
            Err(e) => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))
        }
    }

    /// Change the reference plane of the simulation.
    ///
    /// # Arguments
    ///
    /// * `reference_plane` - The new reference plane of the simulation.
    ///
    /// # Returns
    ///
    /// * `Result<(), &'static str>` - The result of the operation.
    pub fn set_epoch(&mut self, epoch: &PyTime) {
        self.inner.epoch = epoch.inner.clone();
    }

    /// Set the reference plane of the simulation.
    ///
    /// # Arguments
    ///
    /// * `reference_plane` - The new reference plane of the simulation.
    pub fn set_reference_plane(&mut self, reference_plane: &str) {
        let reference_plane = ReferencePlane::from_str(reference_plane).unwrap();
        self.inner.reference_plane = reference_plane;
    }

    /// Set the origin of the simulation.
    ///
    /// # Arguments
    ///
    /// * `origin` - The new origin of the simulation.
    pub fn set_origin(&mut self, origin: &str) -> PyResult<()> {
        match Origin::from_str(origin) {
            Ok(o) => {
                self.inner.origin = o;
                Ok(())
            },
            Err(e) => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))
        }
    }

    /// Set the integrator.
    ///
    /// # Arguments
    ///
    /// * `integrator` - The new integrator.
    pub fn set_integrator(&mut self, integrator: PyRef<PyIntegrator>) {
        self.inner.integrator = integrator.inner.clone();
    }

    /// Add a force to the simulation.
    ///
    /// # Arguments
    ///
    /// * `force` - The force to add.
    pub fn add_force(&mut self, force: PyRef<PyForce>) {
        self.inner.add_force(force.inner.clone());
    }

    /// Calculate the total energy of the simulation.
    pub fn energy(&self) -> f64 {
        self.inner.energy()
    }

    /// Get a single particle from the simulation by name.
    pub fn get_particle(&self, name: &str) -> PyResult<PySpaceRock> {
        let rock = self.inner.get_particle(name);
        match rock {
            Ok(r) => Ok(PySpaceRock { inner: r.clone() }),
            Err(e) => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))
        }
       
    }


    // #[getter]
    // pub fn particles(&self) -> RockCollection {
    //     RockCollection { rocks: self.inner.particles.clone(), name_hash_map: self.inner.particle_index_map.clone() }
    // }


    #[getter]
    pub fn epoch(&self) -> PyTime {
        PyTime { inner: self.inner.epoch.clone() }
    }

    #[getter]
    pub fn reference_plane(&self) -> String {
        self.inner.reference_plane.to_string()
    }


    #[getter]
    pub fn timestep(&self) -> f64 {
        self.inner.integrator.timestep()
    }

    #[getter]
    pub fn origin(&self) -> PyOrigin {
        let o = self.inner.origin.clone();
        PyOrigin { inner: o }
    }


    // repr
    pub fn __repr__(&self) -> String {
        // format!("Simulation:\n");
        // format!("Epoch: {}\n", self.inner.epoch);
        // format!("Reference Plane: {}\n", self.inner.reference_plane);
        // format!("Origin: {}\n", self.inner.origin);
        // format!("Timestep: {}\n", self.inner.integrator.timestep());
        // format!("Forces: {:?}\n", self.inner.forces);
        // format!("Integrator: {:?}\n", self.inner.integrator);
        // format!("Particles: {}\n", self.inner.particles.len());

        let mut s = String::new();
        s.push_str("Simulation:\n");
        s.push_str(&format!("    Epoch: {}\n", self.inner.epoch));
        s.push_str(&format!("    Reference Plane: {}\n", self.inner.reference_plane));
        s.push_str(&format!("    Origin: {}\n", self.inner.origin));
        s.push_str(&format!("    Timestep: {}\n", self.inner.integrator.timestep()));
        s
    }

    #[getter]
    pub fn integrator(&self) -> PyIntegrator {
        PyIntegrator { inner: self.inner.integrator.clone() }
    }

}