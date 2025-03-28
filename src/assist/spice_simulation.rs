//! Set up and manage a collection of gravitationally interacting particles.
use std::collections::HashMap;
use std::sync::Arc;

use crate::SpaceRock;
use crate::constants::GRAVITATIONAL_CONSTANT;
use crate::time::Time;
use crate::{ReferencePlane, Origin};
use crate::errors::SimulationError;

use crate::assist::forces::{Force, NewtonianGravity};
use crate::assist::integrators::{Integrator, IAS15};

use crate::spice::{SpiceBody, SpiceKernel};

use nalgebra::Vector3;

#[derive(Debug, Clone)]
pub struct VariationalParticle {
    pub position: Vector3<f64>,
    pub velocity: Vector3<f64>,
    pub acceleration: Vector3<f64>,
}

#[derive(Debug, Clone)]
pub struct SimulationParticle {
    pub mass: f64,
    pub position: Vector3<f64>,
    pub velocity: Vector3<f64>,
    pub acceleration: Vector3<f64>,
    pub epoch: f64,
    pub variational_particles: Option<Vec<VariationalParticle>>,
    pub interpolation_coefficients: Option<Vec<f64>>,
}

#[derive(Debug, Clone)]
pub struct SimulationState {
    pub epoch: f64,
    pub particles: Vec<SimulationParticle>,
    pub spice_particles: Vec<SimulationParticle>,
    pub particle_index_map: HashMap<String, usize>,
    pub spice_particle_index_map: HashMap<String, usize>,
    pub spice_bodies: Vec<SpiceBody>,
    pub kernel: Arc<SpiceKernel>,
    pub reference_plane: ReferencePlane,
    pub origin: SpiceBody,
}

pub struct SpiceSimulation {
    pub state: SimulationState,
    pub integrator: Box<dyn Integrator + Send + Sync>,
    pub forces: Vec<Box<dyn Force + Send + Sync>>,
}

impl SpiceSimulation {

    pub fn horizons(
        epoch: &Time,
        reference_plane: &str,
        origin: &str,
        kernel: &Arc<SpiceKernel>,
    ) -> Result<SpiceSimulation, Box<dyn std::error::Error>> {
        let reference_plane = ReferencePlane::from_str(reference_plane)?;
        let origin = SpiceBody::from_name(origin)?;

        let mut state = SimulationState {
            epoch: epoch.tdb().jd(),
            particles: Vec::new(),
            spice_particles: Vec::new(),
            particle_index_map: HashMap::new(),
            spice_particle_index_map: HashMap::new(),
            spice_bodies: Vec::new(),
            kernel: kernel.clone(),
            reference_plane,
            origin,
        };

        let mut sim = SpiceSimulation {
            state,
            integrator: Box::new(IAS15::new(0.001)),
            forces: vec![Box::new(NewtonianGravity)],
        };

        let names = ["sun", "mercury barycenter", "venus barycenter", "earth", "moon", "mars barycenter", "jupiter barycenter", 
                     "saturn barycenter", "uranus barycenter", "neptune barycenter", "pluto barycenter", 
                     "2000001", 
                     "2000002", 
                     "2000003", 
                     "2000004", 
                     "2000007",
                     "2000010", 
                     "2000015", 
                     "2000016", 
                     "2000031", 
                    //  "2000048", 
                     "2000052", 
                     "2000065", 
                     "2000087",
                     "2000088", 
                     "2000107",
                    //  "2000451", 
                     "2000511", 
                     "2000704"];

        for name in names.iter() {
            // uppercase the name
            let name = name.to_uppercase();
            let particle = SpiceBody::from_name(&name)?;
            sim.add_spice_body(particle)?;
        }

        Ok(sim)
    }

    /// Instantiate a simulation with the solar system giants.
    ///
    /// # Arguments
    ///
    /// * `epoch` - The epoch of the simulation.
    /// * `reference_plane` - The reference plane of the simulation.
    /// * `origin` - The origin of the simulation.
    ///
    /// # Returns
    ///
    /// * `Result<Simulation, Box<dyn std::error::Error>>` - The simulation with the solar system giants.
    pub fn giants(epoch: &Time, reference_plane: &str, origin: &str, kernel: &Arc<SpiceKernel>) -> Result<SpiceSimulation, Box<dyn std::error::Error>> {

        let reference_plane = ReferencePlane::from_str(reference_plane)?;
        let origin = SpiceBody::from_name(origin)?;

        let mut state = SimulationState {
            epoch: epoch.tdb().jd(),
            particles: Vec::new(),
            spice_particles: Vec::new(),
            particle_index_map: HashMap::new(),
            spice_particle_index_map: HashMap::new(),
            spice_bodies: Vec::new(),
            kernel: kernel.clone(),
            reference_plane,
            origin,
        };

        let mut sim = SpiceSimulation {
            state,
            integrator: Box::new(IAS15::new(0.001)),
            forces: vec![Box::new(NewtonianGravity)],
        };

        // add sun, jupiter barycenter, saturn barycenter, uranus barycenter, neptune barycenter.
        for name in ["sun", "jupiter barycenter", "saturn barycenter", "uranus barycenter", "neptune barycenter"].iter() {
            // uppercase the name
            let name = name.to_uppercase();
            let particle = SpiceBody::from_name(&name)?;
            sim.add_spice_body(particle);
        }
        Ok(sim)
    }

    pub fn add(&mut self, body: SpaceRock) -> Result<(), Box<dyn std::error::Error>> {
        // Check if the body is already in the simulation
        // if self.state.particle_index_map.contains_key(&body.name) {
        //     return Err(SimulationError::BodyAlreadyExists(body.name.clone()));
        // }

        // Create a new simulation particle
        let simulation_particle = SimulationParticle {
            mass: body.mass(),
            position: body.position,
            velocity: body.velocity,
            acceleration: Vector3::zeros(),
            epoch: self.state.epoch,
            variational_particles: None,
            interpolation_coefficients: None,
        };

        // Add the simulation particle to the simulation
        self.state.particles.push(simulation_particle);
        self.state.particle_index_map.insert(body.name.clone(), self.state.particles.len() - 1);

        Ok(())
    }

    pub fn add_spice_body(&mut self, body: SpiceBody) -> Result<(), Box<dyn std::error::Error>> {
        // Check if the body is already in the simulation
        // if self.state.spice_particle_index_map.contains_key(&body.name) {
        //     return Err(SimulationError::BodyAlreadyExists(body.name.clone()));
        // }

        // Get the state of the body at the current epoch
        let (x, y, z, vx, vy, vz) = self.state.kernel.compute_state(&body, &self.state.origin, self.state.epoch)?;
        let position = Vector3::new(x, y, z);
        let velocity = Vector3::new(vx, vy, vz);

        // Create a new spice particle
        let spice_particle = SimulationParticle {
            mass: body.mass,
            position: position,
            velocity: velocity,
            acceleration: Vector3::zeros(),
            epoch: self.state.epoch,
            variational_particles: None,
            interpolation_coefficients: None,
        };

        // Add the spice particle to the simulation
        self.state.spice_bodies.push(body.clone());
        self.state.spice_particles.push(spice_particle);
        self.state.spice_particle_index_map.insert(body.name.clone(), self.state.spice_particles.len() - 1);

        Ok(())
    }

    /// Add a force to the simulation.
    ///
    /// # Arguments
    ///
    /// * `force` - The force to add to the simulation.
    pub fn add_force(&mut self, force: Box<dyn Force + Send + Sync>) {
        self.forces.push(force);
    }


    pub fn step(&mut self) {
        self.integrator.step(&mut self.state, &self.forces);
    }
}