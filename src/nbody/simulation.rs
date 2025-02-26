//! Set up and manage a collection of gravitationally interacting particles.
use std::collections::HashMap;

use crate::SpaceRock;
use crate::constants::GRAVITATIONAL_CONSTANT;
use crate::time::Time;
use crate::{ReferencePlane, Origin};
use crate::errors::SimulationError;


use crate::nbody::forces::{Force, NewtonianGravity};
use crate::nbody::integrators::{Integrator, IAS15};


use nalgebra::Vector3;

/// A simulation maintains:
/// - A collection of particles and their states
/// - The current simulation epoch
/// - Reference frame and origin specifications
/// - Integration method and forces
#[derive(Clone)]
pub struct Simulation {
    pub particles: Vec<SpaceRock>,
    pub epoch: Time,
    pub particle_index_map: HashMap<String, usize>,

    pub reference_plane: ReferencePlane,
    pub origin: Origin,

    pub integrator: Box<dyn Integrator + Send + Sync>,
    pub forces: Vec<Box<dyn Force + Send + Sync>>,
}

impl Default for Simulation {
    fn default() -> Self {
        Self::new(&Time::now(), "J2000", "SSB").unwrap()
    }
}

impl Simulation {
    /// Creates a new simulation at the specified epoch and reference frame.
    pub fn new(epoch: &Time, reference_plane: &str, origin: &str) -> Result<Simulation, Box<dyn std::error::Error>> {
        let mut t = Time::now();
        t.to_tdb();

        let reference_plane = ReferencePlane::from_str(reference_plane)?;
        let origin = Origin::from_str(origin)?;

        Ok(Simulation {
            particles: Vec::new(), 
            epoch: epoch.clone(),
            forces: vec![Box::new(NewtonianGravity)],
            reference_plane: reference_plane,  
            origin: origin,
            integrator: Box::new(IAS15::new(1.0)),
            particle_index_map: HashMap::new()
        })
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
    pub fn giants(epoch: &Time, reference_plane: &str, origin: &str) -> Result<Simulation, Box<dyn std::error::Error>> {

        let mut sim = Simulation::new(epoch, reference_plane, origin)?;
        sim.epoch = epoch.clone();
        sim.epoch.to_tdb();
        sim.integrator = Box::new(IAS15::new(1.0));

        // add sun, jupiter barycenter, saturn barycenter, uranus barycenter, neptune barycenter.
        for name in ["sun", "jupiter barycenter", "saturn barycenter", "uranus barycenter", "neptune barycenter"].iter() {
            let particle = SpaceRock::from_spice(name, epoch, reference_plane, origin)?;
            sim.add(particle)?;
        }
        Ok(sim)
    }

    /// Instantiate a simulation with the solar system planets.
    /// Includes the sun, mercury barycenter, venus barycenter, earth barycenter, mars barycenter, jupiter barycenter, saturn barycenter, uranus barycenter, neptune barycenter.
    ///
    /// # Arguments
    ///
    /// * `epoch` - The epoch of the simulation.
    /// * `reference_plane` - The reference plane of the simulation.
    /// * `origin` - The origin of the simulation.
    ///
    /// # Returns
    ///
    /// * `Result<Simulation, Box<dyn std::error::Error>>` - The simulation with the solar system planets.
    pub fn planets(epoch: &Time, reference_plane: &str, origin: &str) -> Result<Simulation, Box<dyn std::error::Error>> {
        let mut sim = Simulation::new(epoch, reference_plane, origin)?;
        sim.epoch = epoch.clone();
        sim.epoch.to_tdb();
        sim.integrator = Box::new(IAS15::new(1.0));

        let names = ["sun", "mercury barycenter", "venus barycenter", "earth barycenter", "mars barycenter", "jupiter barycenter", 
                     "saturn barycenter", "uranus barycenter", "neptune barycenter"];
        for name in names.iter() {
            let particle = SpaceRock::from_spice(name, epoch, reference_plane, origin)?;
            sim.add(particle)?;
        }
        Ok(sim)
    }

    /// Instantiate a simulation with the solar system planets and moons.
    /// Includes the sun, mercury barycenter, venus barycenter, earth, moon, mars barycenter, jupiter barycenter, saturn barycenter, uranus barycenter, neptune barycenter, pluto barycenter.
    ///
    /// # Arguments
    ///
    /// * `epoch` - The epoch of the simulation.
    /// * `reference_plane` - The reference plane of the simulation.
    /// * `origin` - The origin of the simulation.
    ///
    /// # Returns
    ///
    /// * `Result<Simulation, Box<dyn std::error::Error>>` - The simulation with the solar system planets and moons.
    pub fn horizons(epoch: &Time, reference_plane: &str, origin: &str) -> Result<Simulation, Box<dyn std::error::Error>> {
        let mut sim = Simulation::new(epoch, reference_plane, origin)?;
        sim.epoch = epoch.clone();
        sim.epoch.to_tdb();
        sim.integrator = Box::new(IAS15::new(0.001));

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
            let particle = SpaceRock::from_spice(name, epoch, reference_plane, origin)?;
            sim.add(particle)?;
        }
        Ok(sim)
    }

    
    /// Add a particle to the simulation.
    ///
    /// # Arguments
    ///
    /// * `particle` - The particle to add to the simulation.
    pub fn add(&mut self, mut particle: SpaceRock) -> Result<(), Box<dyn std::error::Error>> {

        if self.epoch.tdb().jd() != particle.epoch.tdb().jd() {
            let err = SimulationError::EpochMismatch(particle.epoch.clone(), self.epoch.clone(), particle.name.clone());
            return Err(err.into());
        }

        if particle.origin.clone() != self.origin {
            if !self.particle_index_map.contains_key(&particle.origin.to_string()) {
                let err = SimulationError::OriginMismatch(particle.origin.clone(), self.origin.clone(), particle.name.clone());
                return Err(err.into());
            }
            let origin = &self.particles[self.particle_index_map[&particle.origin.to_string()]];
            particle.change_origin(origin);
            println!("Changing origin of {} from {} to {}", particle.name, particle.origin, origin.name);
        }

        particle.change_reference_plane(self.reference_plane.as_str())?;
        particle.epoch.to_tdb();
        // self.particle_index_map.insert((*particle.name).to_string(), self.particles.len());
        // self.particles.push(particle);

        if particle.mass() == 0.0 {
            self.particle_index_map.insert((*particle.name).to_string(), self.particles.len());
            self.particles.push(particle);
            return Ok(());
        }

        self.particle_index_map.insert((*particle.name).to_string(), self.particles.len());
        self.particles.push(particle);
        
        // make sure to sort the particles by mass
        self.particles.sort_by(|a, b| b.mass().partial_cmp(&a.mass()).unwrap());
        // update the particle index map
        for (idx, particle) in self.particles.iter().enumerate() {
            self.particle_index_map.insert((*particle.name).to_string(), idx);
        }

        Ok(())
    }

    /// Remove a particle from the simulation.
    ///
    /// # Arguments
    ///
    /// * `name` - The name of the particle to remove.
    pub fn remove(&mut self, name: &str) -> Result<(), SimulationError> {
        if self.particle_index_map.contains_key(name) {
            let idx = self.particle_index_map[name];
            self.particles.remove(idx);
            self.particle_index_map.remove(name);
            for value in self.particle_index_map.values_mut() {
                if *value > idx {
                    *value -= 1;
                }
            }
        } else {
            return Err(SimulationError::ParticleNotFound(name.to_string()));
        }
        Ok(())
    }

    /// Move the simulation to the center of mass.
    pub fn move_to_center_of_mass(&mut self) -> Result<(), Box<dyn std::error::Error>> {
        let mut total_mass = 0.0;
        let mut center_of_mass = Vector3::new(0.0, 0.0, 0.0);
        let mut center_of_mass_velocity = Vector3::new(0.0, 0.0, 0.0);

        for particle in &self.particles {
            if particle.mass() == 0.0 {
                continue;
            }
            center_of_mass += particle.mass() * particle.position;
            center_of_mass_velocity += particle.mass() * particle.velocity;
            total_mass += particle.mass();
        }

        center_of_mass /= total_mass;
        center_of_mass_velocity /= total_mass;

        let x = center_of_mass.x;
        let y = center_of_mass.y;
        let z = center_of_mass.z;
        let vx = center_of_mass_velocity.x;
        let vy = center_of_mass_velocity.y;
        let vz = center_of_mass_velocity.z;
       
        let mut origin_rock = SpaceRock::from_xyz("simulation_barycenter", 
                                                  x, y, z, vx, vy, vz, 
                                                  self.epoch.clone(), 
                                                  self.reference_plane.as_str(),
                                                  self.origin.as_str())?;
        origin_rock.set_mass(total_mass);

        for particle in &mut self.particles {
            particle.change_origin(&origin_rock);
        }

        let origin = Origin::new_custom(total_mass * GRAVITATIONAL_CONSTANT, "simulation_barycenter");
        self.origin = origin;
        Ok(())
    }

    /// Change the origin of the simulation.
    ///
    /// # Arguments
    ///     
    /// * `origin` - The name of the particle to set as the origin.
    pub fn change_origin(&mut self, origin: &str) -> Result<(), String> {

        if !self.particle_index_map.contains_key(origin) {
           return Err(format!("Origin {} not found in perturbers", origin));
        }

        let new_origin = Origin::new_custom(self.particles[self.particle_index_map[origin]].mass() * GRAVITATIONAL_CONSTANT, origin);

        self.origin = new_origin;

        let origin_position = self.particles[self.particle_index_map[origin]].position;
        let origin_velocity = self.particles[self.particle_index_map[origin]].velocity;

        for particle in &mut self.particles {
            particle.position -= origin_position;
            particle.velocity -= origin_velocity;
        }

        Ok(())
    }

    /// Step the simulation forward in time by one timestep.
    pub fn step(&mut self) {
        self.integrator.step(&mut self.particles, &mut self.epoch, &self.forces);
    }

    /// Integrate the simulation to a new epoch.
    ///
    /// Advances the simulation time by taking steps until reaching the target epoch.
    /// The integration direction (forward/backward) is determined automatically based
    /// on the difference between current and target epochs. The timestep is adjusted
    /// automatically when approaching the target epoch to hit it exactly.
    ///
    /// # Arguments
    ///
    /// * `epoch` - The new epoch to integrate to.
    pub fn integrate(&mut self, epoch: &Time) {

        
        let dt = epoch.tdb().jd() - self.epoch.tdb().jd();
        if dt.abs() < 1e-16 {
            return;
        }

        if dt < 0.0 && self.integrator.timestep() > 0.0 {
            self.integrator.set_timestep(-self.integrator.timestep());
        }


        loop {
            let dt = epoch.tdb().jd() - self.epoch.tdb().jd();

            // done integrating
            if dt.abs() < 1e-16 {
                break;
            }

            // if we're within a timestep of the epoch, just take a step of that size
            if dt.abs() < self.integrator.timestep().abs() {
                let last_timestep = self.integrator.timestep().clone();
                self.integrator.set_timestep(dt);
                self.step();
                self.integrator.set_timestep(last_timestep);
                continue;

                // let dt = epoch.tdb().jd() - self.epoch.tdb().jd();
                // if dt.abs() < 1e-16 {
                //     self.integrator.set_timestep(last_timestep);
                //     break;
                // }
                
            }

            // if the timestep is negative, make sure the integrator is set to negative
            if dt < 0.0 {
                if self.integrator.timestep() > 0.0 {
                    self.integrator.set_timestep(-self.integrator.timestep());
                }
            } else if self.integrator.timestep() < 0.0 {
                self.integrator.set_timestep(-self.integrator.timestep());
            }
            self.step();
        }
        
        // let dt = &epoch - &self.epoch;
        // let dt = epoch.tdb().jd() - self.epoch.tdb().jd();
        // if dt.abs() < 1e-16 {
        //     return;
        // }
        // // create an exact match for the epoch
        // let old_timestep = self.integrator.timestep();
        // self.integrator.set_timestep(dt);
        // self.step();
        // // reset the timestep
        // self.integrator.set_timestep(old_timestep);
    }

    /// Get a particle from the simulation by name.
    ///
    /// # Arguments
    ///
    /// * `name` - The name of the particle to get.
    ///
    /// # Returns
    ///
    /// * `Result<&SpaceRock, SimulationError>` - The particle with the given name.
    pub fn get_particle(&self, name: &str) -> Result<&SpaceRock, SimulationError> {
        if self.particle_index_map.contains_key(name) { 
            let idx = self.particle_index_map[name];
            let p = &self.particles[idx];
            return Ok(p);
        }
        Err(SimulationError::ParticleNotFound(name.to_string()))
    }

    /// Get the energy of the simulation.
    pub fn energy(&self) -> f64 {
        let mut kinetic_energy = 0.0;
        let mut potential_energy = 0.0;

        for idx in 0..self.particles.len() {
            kinetic_energy += 0.5 * self.particles[idx].mass() * self.particles[idx].velocity.norm_squared();
            for jdx in (idx + 1)..self.particles.len() {
                let r = (self.particles[idx].position - self.particles[jdx].position).norm();
                potential_energy -= GRAVITATIONAL_CONSTANT * self.particles[idx].mass() * self.particles[jdx].mass() / r;
            }
        }
        kinetic_energy + potential_energy
    }

    /// Add a force to the simulation.
    ///
    /// # Arguments
    ///
    /// * `force` - The force to add to the simulation.
    pub fn add_force(&mut self, force: Box<dyn Force + Send + Sync>) {
        self.forces.push(force);
    }

}