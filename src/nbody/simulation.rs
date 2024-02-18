use std::collections::HashMap;

use crate::spacerock::SpaceRock;
use crate::constants::gravitational_constant;
use crate::time::Time;
use crate::nbody::integrators::Integrator;
use crate::nbody::leapfrog::Leapfrog;
use crate::nbody::freezer::Freezer;

use nalgebra::Vector3;

use spice;

#[derive(Debug, Clone, PartialEq)]
pub struct Simulation {
    pub perturbers: Vec<SpaceRock>,
    pub particles: Vec<SpaceRock>,
    pub epoch: Time,
    pub integrator: Integrator,
    pub particle_index_map: HashMap<String, usize>,
    pub perturber_index_map: HashMap<String, usize>,

    pub frame: Option<String>,
    pub origin: Option<String>,
    pub mass: f64,
}

impl Simulation {

    pub fn new() -> Simulation {
        Simulation {perturbers: Vec::new(), 
                    particles: Vec::new(), 
                    epoch: Time::now(), 
                    frame: None,
                    origin: None,
                    mass: 0.0,
                    // integrator: Integrator::Leapfrog(Leapfrog::new(20.0)),
                    integrator: Integrator::Freezer(Freezer::new(20.0)),
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
            particle.change_origin(origin.name.as_str(), &origin.position, &origin.velocity, &(origin.mass * gravitational_constant));
            // log the change of origin to console
            println!("Changing origin of {} from {} to {}", particle.name, particle.origin, origin.name);
        }

        particle.change_frame(self.frame.clone().unwrap().as_str());
        particle.orbit = None;

        if particle.mass > 0.0 {
            self.mass += particle.mass;
            self.perturber_index_map.insert(particle.name.clone(), self.perturbers.len());
            self.perturbers.push(particle);
        } else {
            self.particle_index_map.insert(particle.name.clone(), self.particles.len());
            self.particles.push(particle);
        }

        Ok(())
    }

    pub fn remove(&mut self, name: &str) {
        if self.perturber_index_map.contains_key(name) {
            let idx = self.perturber_index_map[name];
            self.perturbers.remove(idx);
            self.perturber_index_map.remove(name);
        } else if self.particle_index_map.contains_key(name) {
            let idx = self.particle_index_map[name];
            self.mass -= self.particles[idx].mass;
            self.particles.remove(idx);
            self.particle_index_map.remove(name);
        }
    }

    pub fn step(&mut self) {
        self.integrator.step(&mut self.perturbers, &mut self.particles);
        self.epoch += self.integrator.timestep();
    }

    pub fn threaded_step(&mut self) {
        self.integrator.threaded_step(&mut self.perturbers, &mut self.particles);
        self.epoch += self.integrator.timestep();
    }

    pub fn move_to_center_of_mass(&mut self) {
        let mut total_mass = 0.0;
        let mut center_of_mass = Vector3::new(0.0, 0.0, 0.0);
        let mut center_of_mass_velocity = Vector3::new(0.0, 0.0, 0.0);

        for perturber in &self.perturbers {
            center_of_mass += perturber.mass * perturber.position;
            center_of_mass_velocity += perturber.mass * perturber.velocity;
        }

        center_of_mass /= self.mass;
        center_of_mass_velocity /= self.mass;

        let origin = "ssb";
        let mu = self.mass * gravitational_constant;

        for perturber in &mut self.perturbers {
            perturber.change_origin(origin, &center_of_mass, &center_of_mass_velocity, &mu);
        }

        for particle in &mut self.particles {
            particle.change_origin(origin, &center_of_mass, &center_of_mass_velocity, &mu);
        }
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

    pub fn integrate(&mut self, epoch: &Time) {

        let mut epoch = epoch.clone();
        epoch.change_timescale(self.epoch.timescale.as_str());

        let old_timestep = self.integrator.timestep();
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
        self.integrator.set_timestep(&epoch - &self.epoch);
        self.step();

        // reset the timestep
        self.integrator.set_timestep(old_timestep);

    }

    //pub fn get_particle(&mut self, name: &str) -> Option<&SpaceRock> {
    // pub fn get_particle(&self, name: &str) -> Option<SpaceRock> {
    //     if self.test_particle_index_map.contains_key(name) { 
    //         let idx = self.test_particle_index_map[name];
    //         let p = &self.test_particles[idx];
    //         let mut rock = self.test_spacerocks[idx].clone();
    //         rock.position = p.position;
    //         rock.velocity = p.velocity;
    //         rock.epoch.epoch = self.time;
    //         rock.calculate_kepler();
    //         return Some(rock);
    //     }
    //     if self.perturber_index_map.contains_key(name) {
    //         let idx = self.perturber_index_map[name];
    //         let p = &self.perturbers[idx];
    //         let mut rock = self.perturber_spacerocks[idx].clone();
    //         rock.position = p.position;
    //         rock.velocity = p.velocity;
    //         rock.epoch.epoch = self.time;
    //         rock.calculate_kepler();
    //         return Some(rock);
    //     }
    //     return None;
    // }

    pub fn energy(&self) -> f64 {
        let mut kinetic_energy = 0.0;
        let mut potential_energy = 0.0;

        for idx in 0..self.perturbers.len() {
            kinetic_energy += 0.5 * self.perturbers[idx].mass * self.perturbers[idx].velocity.norm_squared();
            for jdx in (idx + 1)..self.perturbers.len() {
                let r = (self.perturbers[idx].position - self.perturbers[jdx].position).norm();
                potential_energy -= gravitational_constant * self.perturbers[idx].mass * self.perturbers[jdx].mass / r;
            }
        }
        return kinetic_energy + potential_energy;
    }

}