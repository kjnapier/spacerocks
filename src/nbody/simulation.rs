use std::collections::HashMap;

use crate::SpaceRock;
use crate::constants::GRAVITATIONAL_CONSTANT;
use crate::time::Time;


use crate::nbody::forces::{Force, NewtonianGravity, Drag, SolarGR};
use crate::nbody::integrators::{Integrator, Leapfrog, IAS15};


use nalgebra::Vector3;

use rayon::prelude::*;

pub struct Simulation {
    pub particles: Vec<SpaceRock>,
    pub epoch: Time,
    pub particle_index_map: HashMap<String, usize>,

    pub frame: Option<String>,
    pub origin: Option<String>,

    pub integrator: Box<dyn Integrator + Send + Sync>,
    pub forces: Vec<Box<dyn Force + Send + Sync>>,
}

impl Simulation {

    pub fn new() -> Simulation {
        Simulation {particles: Vec::new(), 
                    epoch: Time::now(), 
                    forces: vec![],
                    frame: None,
                    origin: None,
                    integrator: Box::new(IAS15::new(1.0)),
                    particle_index_map: HashMap::new()}
    }

    pub fn add(&mut self, mut particle: SpaceRock) -> Result<(), String> {

        // if the origin is None, set it to be the origin of the first rock
        if self.origin.is_none() {
            self.origin = Some(particle.origin.clone());
            println!("Setting origin to {}", particle.origin.clone());
        }

        if self.frame.is_none() {
            self.frame = Some(particle.frame.clone());
            println!("Setting frame to {}", particle.frame);
        }

        if self.epoch != particle.epoch {
            return Err(format!("The epoch of {} particle ({:?}) did not match the simulation epoch", particle.name, particle.epoch));
        }

        if Some(particle.origin.clone()) != self.origin {
            if !self.particle_index_map.contains_key(&particle.origin.clone()) {
                return Err(format!("The origin of the particle ({}) did not match the simulation origin, and was not found in perturbers", particle.origin.clone()));
            }

            let origin = &self.particles[self.particle_index_map[&particle.origin.clone()]];
            particle.change_origin(&origin);

            // log the change of origin to console
            println!("Changing origin of {} from {} to {}", particle.name, particle.origin, origin.name);
        }

        particle.change_frame(self.frame.clone().unwrap().as_str());
        particle.orbit = None;

        self.particle_index_map.insert(particle.name.clone(), self.particles.len());
        self.particles.push(particle);

        Ok(())
    }

    pub fn remove(&mut self, name: &str) -> Result<(), String> {
        if self.particle_index_map.contains_key(name) {
            let idx = self.particle_index_map[name];
            self.particles.remove(idx);
            self.particle_index_map.remove(name);
        } else {
            return Err(format!("No particle found with name {}", name));
        }
        Ok(())
    }
   

    pub fn zero_accelerations(&mut self) {
        for particle in &mut self.particles {
            particle.acceleration = Vector3::new(0.0, 0.0, 0.0);
        }
    }

    pub fn step(&mut self) {
        // Zero the accelerations of each particle
        self.zero_accelerations();
        self.integrator.step(&mut self.particles, &mut self.epoch, &self.forces);
    }

    pub fn move_to_center_of_mass(&mut self) {
        let mut total_mass = 0.0;
        let mut center_of_mass = Vector3::new(0.0, 0.0, 0.0);
        let mut center_of_mass_velocity = Vector3::new(0.0, 0.0, 0.0);

        for particle in &self.particles {
            if particle.mass == 0.0 {
                continue;
            }
            center_of_mass += particle.mass * particle.position;
            center_of_mass_velocity += particle.mass * particle.velocity;
            total_mass += particle.mass;
        }

        center_of_mass /= total_mass;
        center_of_mass_velocity /= total_mass;

        let x = center_of_mass.x;
        let y = center_of_mass.y;
        let z = center_of_mass.z;
        let vx = center_of_mass_velocity.x;
        let vy = center_of_mass_velocity.y;
        let vz = center_of_mass_velocity.z;
       
        let mut origin = SpaceRock::from_xyz("ssb", 
                                             x, y, z, vx, vy, vz, 
                                             self.epoch.clone(), 
                                             self.frame.clone().unwrap().as_str(), 
                                             self.origin.clone().unwrap().as_str());
        origin.mass = total_mass;

        for particle in &mut self.particles {
            particle.change_origin(&origin);
        }

        self.origin = Some("ssb".to_string());
    }

    pub fn change_origin(&mut self, origin: &str) -> Result<(), String> {

        if !self.particle_index_map.contains_key(origin) {
           return Err(format!("Origin {} not found in perturbers", origin));
        }

        self.origin = Some(origin.to_string());

        let origin_position = self.particles[self.particle_index_map[origin]].position;
        let origin_velocity = self.particles[self.particle_index_map[origin]].velocity;

        for particle in &mut self.particles {
            particle.position -= origin_position;
            particle.velocity -= origin_velocity;
        }

        Ok(())
    }

    pub fn integrate(&mut self, epoch: &Time) {

        let mut epoch = epoch.clone();
        epoch.change_timescale(self.epoch.timescale.as_str());

        let dt = &epoch - &self.epoch;

        if dt == 0.0 {
            return;
        }

        if dt < 0.0 {
            self.integrator.set_timestep(-self.integrator.timestep());
        }

        while (&self.epoch - &epoch).abs() >= self.integrator.timestep().abs() {
            self.step();
        }       

        // create an exact match for the epoch
        let old_timestep = self.integrator.timestep();
        self.integrator.set_timestep(&epoch - &self.epoch);
        self.step();

        // reset the timestep
        self.integrator.set_timestep(old_timestep);
    }

    pub fn get_particle(&self, name: &str) -> Result<&SpaceRock, Box<dyn std::error::Error>> {
        if self.particle_index_map.contains_key(name) { 
            let idx = self.particle_index_map[name];
            let p = &self.particles[idx];
            return Ok(p);
        }
        return Err(format!("{} not found in simulation", name).into());
    }

    pub fn energy(&self) -> f64 {
        let mut kinetic_energy = 0.0;
        let mut potential_energy = 0.0;

        for idx in 0..self.particles.len() {
            kinetic_energy += 0.5 * self.particles[idx].mass * self.particles[idx].velocity.norm_squared();
            for jdx in (idx + 1)..self.particles.len() {
                let r = (self.particles[idx].position - self.particles[jdx].position).norm();
                potential_energy -= GRAVITATIONAL_CONSTANT * self.particles[idx].mass * self.particles[jdx].mass / r;
            }
        }
        return kinetic_energy + potential_energy;
    }

    pub fn add_force(&mut self, force: Box<dyn Force + Send + Sync>) {
        self.forces.push(force);
    }

}