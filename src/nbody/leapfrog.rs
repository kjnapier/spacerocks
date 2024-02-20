
use crate::SpaceRock;
use crate::nbody::acceleration::{calculate_and_set_test_particle_accelerations, calculate_and_set_perturber_accelerations};

use rayon::prelude::*;

#[derive(PartialEq, Debug, Clone)]
pub struct Leapfrog {
    pub timestep: f64,
}

impl Leapfrog {
    pub fn new(timestep: f64) -> Leapfrog {
        Leapfrog { timestep }
    }

    pub fn step(&self, perturbers: &mut Vec<SpaceRock>, particles: &mut Vec<SpaceRock>) {
        // drift
        for particle in &mut *particles {
            particle.position += particle.velocity * 0.5 * self.timestep;
            particle.epoch += 0.5 * self.timestep;
        }

        // kick + drift
        calculate_and_set_test_particle_accelerations(perturbers, particles);       
        for particle in &mut *particles {
            particle.velocity += self.timestep * particle.acceleration;
            particle.position += particle.velocity * 0.5 * self.timestep;
            particle.epoch += 0.5 * self.timestep;
        }
        
        // drift
        for mut perturber in &mut *perturbers {
            perturber.position += perturber.velocity * 0.5 * self.timestep;
            perturber.epoch += 0.5 * self.timestep;
        }

        // kick + drift
        calculate_and_set_perturber_accelerations(perturbers);
        for mut perturber in &mut *perturbers {
            perturber.velocity += self.timestep * perturber.acceleration;
            perturber.position += perturber.velocity * 0.5 * self.timestep;
            perturber.epoch += 0.5 * self.timestep;
        }
       
    }

    pub fn threaded_step(&self, perturbers: &mut Vec<SpaceRock>, particles: &mut Vec<SpaceRock>) {
        // drift
        particles.par_iter_mut().for_each(|particle| {
            particle.position += particle.velocity * 0.5 * self.timestep;
            particle.epoch += 0.5 * self.timestep;
        });

        // kick + drift
        calculate_and_set_test_particle_accelerations(perturbers, particles);
        particles.par_iter_mut().for_each(|particle| {
            particle.velocity += self.timestep * particle.acceleration;
            particle.position += particle.velocity * 0.5 * self.timestep;
            particle.epoch += 0.5 * self.timestep;
        });

        // drift
        perturbers.par_iter_mut().for_each(|perturber| {
            perturber.position += perturber.velocity * 0.5 * self.timestep;
            perturber.epoch += 0.5 * self.timestep;
        });

        // kick + drift
        calculate_and_set_perturber_accelerations(perturbers);
        perturbers.par_iter_mut().for_each(|perturber| {
            perturber.velocity += self.timestep * perturber.acceleration;
            perturber.position += perturber.velocity * 0.5 * self.timestep;
            perturber.epoch += 0.5 * self.timestep;
        });
    }

}