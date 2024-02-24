use std::collections::HashMap;

use crate::SpaceRock;
use crate::constants::GRAVITATIONAL_CONSTANT;
use crate::time::Time;
use crate::nbody::integrator::Integrator;
use crate::nbody::leapfrog::Leapfrog;

use std::cell::RefCell;

use nalgebra::Vector3;

#[derive(Debug, Clone, PartialEq)]
pub enum Force {
    NewtonianGravity,
}

#[derive(Debug)]
pub struct Simulation {
    pub perturbers: Vec<SpaceRock>,
    pub particles: Vec<SpaceRock>,
    pub epoch: Time,
    // pub integrator: Integrator,
    // make integrator be a trait object
    // pub integrator: Box<dyn Integrator>,
    pub particle_index_map: HashMap<String, usize>,
    pub perturber_index_map: HashMap<String, usize>,
    pub forces: Vec<Force>,

    pub frame: Option<String>,
    pub origin: Option<String>,
}

impl Simulation {

    pub fn new() -> Simulation {
        Simulation {perturbers: Vec::new(), 
                    particles: Vec::new(), 
                    epoch: Time::now(), 
                    forces: vec![Force::NewtonianGravity],
                    frame: None,
                    origin: None,
                    integrator: Box::new(Leapfrog::new(20.0)),
                    particle_index_map: HashMap::new(),
                    perturber_index_map: HashMap::new()}
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
            if !self.perturber_index_map.contains_key(&particle.origin.clone()) {
                return Err(format!("The origin of the particle ({}) did not match the simulation origin, and was not found in perturbers", particle.origin.clone()));
            }

            let origin = &self.perturbers[self.perturber_index_map[&particle.origin.clone()]];
            particle.change_origin(&origin);

            // log the change of origin to console
            println!("Changing origin of {} from {} to {}", particle.name, particle.origin, origin.name);
        }

        particle.change_frame(self.frame.clone().unwrap().as_str());
        particle.orbit = None;

        if particle.mass > 0.0 {
            self.perturber_index_map.insert(particle.name.clone(), self.perturbers.len());
            self.perturbers.push(particle);
        } else {
            self.particle_index_map.insert(particle.name.clone(), self.particles.len());
            self.particles.push(particle);
        }

        Ok(())
    }

    pub fn remove(&mut self, name: &str) -> Result<(), String> {

        if self.perturber_index_map.contains_key(name) {
            let idx = self.perturber_index_map[name];
            self.perturbers.remove(idx);
            self.perturber_index_map.remove(name);
        } else if self.particle_index_map.contains_key(name) {
            let idx = self.particle_index_map[name];
            self.particles.remove(idx);
            self.particle_index_map.remove(name);
        } else {
            return Err(format!("No particle found with name {}", name));
        }

        Ok(())

    }

    pub fn apply_forces(&mut self) -> Result<(), String> {
        self.zero_accelerations();
        for force in &self.forces.clone() {
            match force {
                Force::NewtonianGravity => self.apply_newtonian_gravity(),
                _ => Err(format!("Force {:?} not implemented", force))?,
            }
        }
        Ok(())
    }
   
    pub fn apply_newtonian_gravity(&mut self) {
        let n_perturbers = self.perturbers.len();
        for idx in 0..self.perturbers.len() {
            for jdx in (idx + 1)..self.perturbers.len() {
                let r_vec = self.perturbers[idx].position - self.perturbers[jdx].position;
                let r = r_vec.norm();

                let xi = -GRAVITATIONAL_CONSTANT * r_vec / (r * r * r);
                let idx_acceleration = xi * self.perturbers[jdx].mass;
                let jdx_acceleration = -xi * self.perturbers[idx].mass;
                self.perturbers[idx].acceleration += idx_acceleration;
                self.perturbers[jdx].acceleration += jdx_acceleration;
            }

            for particle in self.particles.iter_mut() {
                let r_vec = particle.position - self.perturbers[idx].position;
                let r = r_vec.norm();
                let xi = -GRAVITATIONAL_CONSTANT * r_vec / (r * r * r);
                let a = xi * self.perturbers[idx].mass;
                particle.acceleration += a;
            }
        }
    }

    pub fn zero_accelerations(&mut self) {
        for perturber in &mut self.perturbers {
            perturber.acceleration = Vector3::new(0.0, 0.0, 0.0);
        }
        for particle in &mut self.particles {
            particle.acceleration = Vector3::new(0.0, 0.0, 0.0);
        }
    }

    pub fn step(&mut self, integrator: &dyn Integrator) {
        // use Option to take the integrator out of the box
        
        integrator.step(self);      
    }

    pub fn move_to_center_of_mass(&mut self) {
        let mut total_mass = 0.0;
        let mut center_of_mass = Vector3::new(0.0, 0.0, 0.0);
        let mut center_of_mass_velocity = Vector3::new(0.0, 0.0, 0.0);

        for perturber in &self.perturbers {
            center_of_mass += perturber.mass * perturber.position;
            center_of_mass_velocity += perturber.mass * perturber.velocity;
            total_mass += perturber.mass;
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

        for perturber in &mut self.perturbers {
            perturber.change_origin(&origin);
        }

        for particle in &mut self.particles {
            particle.change_origin(&origin);
        }

        self.origin = Some("ssb".to_string());
    }

    pub fn change_origin(&mut self, origin: &str) -> Result<(), String> {

        if !self.perturber_index_map.contains_key(origin) {
           return Err(format!("Origin {} not found in perturbers", origin));
        }

        self.origin = Some(origin.to_string());
        let origin_position = self.perturbers[self.perturber_index_map[origin]].position;
        let origin_velocity = self.perturbers[self.perturber_index_map[origin]].velocity;

        for perturber in &mut self.perturbers {
            perturber.position -= origin_position;
            perturber.velocity -= origin_velocity;
        }

        for particle in &mut self.particles {
            particle.position -= origin_position;
            particle.velocity -= origin_velocity;
        }

        Ok(())

    }

    pub fn integrate(&mut self, epoch: &Time, integrator: &dyn Integrator) {

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

        if self.perturber_index_map.contains_key(name) {
            let idx = self.perturber_index_map[name];
            let p = &self.perturbers[idx];
            return Ok(p);
        }

        return Err(format!("{} not found in simulation", name).into());
    }

    pub fn energy(&self) -> f64 {
        let mut kinetic_energy = 0.0;
        let mut potential_energy = 0.0;

        for idx in 0..self.perturbers.len() {
            kinetic_energy += 0.5 * self.perturbers[idx].mass * self.perturbers[idx].velocity.norm_squared();
            for jdx in (idx + 1)..self.perturbers.len() {
                let r = (self.perturbers[idx].position - self.perturbers[jdx].position).norm();
                potential_energy -= GRAVITATIONAL_CONSTANT * self.perturbers[idx].mass * self.perturbers[jdx].mass / r;
            }
        }
        return kinetic_energy + potential_energy;
    }

}