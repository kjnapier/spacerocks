use std::collections::HashMap;

use crate::spacerock::SpaceRock;
use crate::nbody::leapfrog::leapfrog_step;
use crate::nbody::simparticle::SimParticle;
use crate::constants::G;
use crate::time::Time;

use nalgebra::Vector3;
use nalgebra;

use spice;

#[derive(Debug, Clone, PartialEq)]
pub struct Simulation {
    pub perturbers: Vec<SimParticle>,
    pub test_particles: Vec<SimParticle>,
    pub time: f64,
    pub timestep: f64,
    pub integrator: String,
    pub test_particle_index_map: HashMap<String, usize>,
    pub perturber_index_map: HashMap<String, usize>,
    pub test_spacerocks: Vec<SpaceRock>,
    pub perturber_spacerocks: Vec<SpaceRock>,
}

impl Simulation {

    // Instantiation Methods
    pub fn new() -> Simulation {
        Simulation {perturbers: Vec::new(), 
                    test_particles: Vec::new(), 
                    perturber_spacerocks: Vec::new(),
                    test_spacerocks: Vec::new(),

                    time: 0.0, 
                    timestep: 1.0, 
                    integrator: "leapfrog".to_string(),
                    test_particle_index_map: HashMap::new(),
                    perturber_index_map: HashMap::new()}
    }    

    pub fn giants(t: &Time) -> Simulation {
        let mut sim = Simulation::new();
        let sun = SpaceRock::from_spice("Sun", t);
        let jupiter = SpaceRock::from_spice("Jupiter Barycenter", t);
        let saturn = SpaceRock::from_spice("Saturn Barycenter", t);
        let uranus = SpaceRock::from_spice("Uranus Barycenter", t);
        let neptune = SpaceRock::from_spice("Neptune Barycenter", t);

        sim.time = t.epoch;
        sim.timestep = 0.001;

        let perturbers = vec![sun, jupiter, saturn, uranus, neptune];
        for mut perturber in perturbers {
            sim.add(perturber);
        }

        sim.move_to_center_of_mass();
        sim
    }

    pub fn add(&mut self, particle: SpaceRock) {

        let mut p = SimParticle::from_spacerock(&particle);
        if p.mass > 0.0 {
            self.perturber_index_map.insert(p.name.clone(), self.perturbers.len());
            self.perturbers.push(p);
            self.perturber_spacerocks.push(particle);
        } else {
            self.test_particle_index_map.insert(p.name.clone(), self.test_particles.len());
            self.test_particles.push(p);
            self.test_spacerocks.push(particle);
        }
    }

    pub fn remove(&mut self, name: &str) {
        if self.perturber_index_map.contains_key(name) {
            let idx = self.perturber_index_map[name];
            self.perturbers.remove(idx);
            self.perturber_index_map.remove(name);
            self.perturber_spacerocks.remove(idx);
        } else if self.test_particle_index_map.contains_key(name) {
            let idx = self.test_particle_index_map[name];
            self.test_particles.remove(idx);
            self.test_particle_index_map.remove(name);
            self.test_spacerocks.remove(idx);
        }
    }

    pub fn move_to_center_of_mass(&mut self) {
        let mut total_mass = 0.0;
        let mut center_of_mass = Vector3::new(0.0, 0.0, 0.0);
        let mut center_of_mass_velocity = Vector3::new(0.0, 0.0, 0.0);
        for perturber in &self.perturbers {
            total_mass += perturber.mass;
            center_of_mass += perturber.mass * perturber.position;
            center_of_mass_velocity += perturber.mass * perturber.velocity;
        }

        center_of_mass /= total_mass;
        center_of_mass_velocity /= total_mass;

        for perturber in &mut self.perturbers {
            perturber.position -= center_of_mass;
            perturber.velocity -= center_of_mass_velocity;
        }

        for test_particle in &mut self.test_particles {
            test_particle.position -= center_of_mass;
            test_particle.velocity -= center_of_mass_velocity;
        }
    }

    pub fn integrate(&mut self, epoch: f64) {
        let mut timesteps = Vec::new();
        let dt = epoch - self.time;
        let mut t = 0.0;

        if dt == 0.0 {
            return;
        }
        
        if dt < 0.0 {
            while t > dt {
                timesteps.push(-self.timestep);
                t -= self.timestep;
            }
        } else {
            while t < dt {
                timesteps.push(self.timestep);
                t += self.timestep;
            }
        }
        if t != dt {
            timesteps.push(dt - t);
        }
        for timestep in timesteps {
            self.step(timestep);
        }
    }

    // pub fn get_particle(&self, name: &str) -> Option<&SimParticle> {
    //     if self.test_particle_index_map.contains_key(name) {
    //         return Some(&self.test_particles[self.test_particle_index_map[name]]);
    //     }
    //     if self.perturber_index_map.contains_key(name) {
    //         return Some(&self.perturbers[self.perturber_index_map[name]]);
    //     }
    //     return None;
    // }

    //pub fn get_particle(&mut self, name: &str) -> Option<&SpaceRock> {
    pub fn get_particle(&self, name: &str) -> Option<SpaceRock> {
        if self.test_particle_index_map.contains_key(name) { 
            let idx = self.test_particle_index_map[name];
            let p = &self.test_particles[idx];
            let mut rock = self.test_spacerocks[idx].clone();
            rock.position = p.position;
            rock.velocity = p.velocity;
            rock.epoch.epoch = self.time;
            rock.calculate_kepler();
            return Some(rock);
        }
        if self.perturber_index_map.contains_key(name) {
            let idx = self.perturber_index_map[name];
            let p = &self.perturbers[idx];
            let mut rock = self.perturber_spacerocks[idx].clone();
            rock.position = p.position;
            rock.velocity = p.velocity;
            rock.epoch.epoch = self.time;
            rock.calculate_kepler();
            return Some(rock);
        }
        return None;
    }

    pub fn step(&mut self, timestep: f64) {
        leapfrog_step(&mut self.perturbers, &mut self.test_particles, timestep);
        self.time += timestep;
    }

    pub fn energy(&self) -> f64 {
        let mut kinetic_energy = 0.0;
        let mut potential_energy = 0.0;

        for idx in 0..self.perturbers.len() {
            kinetic_energy += 0.5 * self.perturbers[idx].mass * self.perturbers[idx].velocity.norm_squared();
            for jdx in (idx + 1)..self.perturbers.len() {
                let r = (self.perturbers[idx].position - self.perturbers[jdx].position).norm();
                potential_energy -= G * self.perturbers[idx].mass * self.perturbers[jdx].mass / r;
            }
        }
        return kinetic_energy + potential_energy;
    }

}