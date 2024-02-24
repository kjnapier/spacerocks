
use crate::SpaceRock;
use crate::Simulation;
use crate::nbody::integrators::Integrator;
use crate::nbody::forces::Force;

use rayon::prelude::*;

#[derive(PartialEq, Debug, Clone, Copy)]
pub struct Leapfrog {
    pub timestep: f64,
}

impl Leapfrog {
    pub fn new(timestep: f64) -> Leapfrog {
        Leapfrog { timestep }
    }
}

impl Integrator for Leapfrog {

    fn step(&mut self, particles: &mut Vec<SpaceRock>, perturbers: &mut Vec<SpaceRock>, forces: &Vec<Box<dyn Force + Send + Sync>>) {
        // drift
        for particle in &mut *particles {
            particle.position += particle.velocity * 0.5 * self.timestep;
            particle.epoch += 0.5 * self.timestep;
        }
        for perturber in &mut *perturbers {
            perturber.position += perturber.velocity * 0.5 * self.timestep;
            perturber.epoch += 0.5 * self.timestep;
        }

        for force in forces {
            force.apply(particles, perturbers);
        }

        for particle in &mut *particles {
            particle.velocity += self.timestep * particle.acceleration;
            particle.position += particle.velocity * 0.5 * self.timestep;
            particle.epoch += 0.5 * self.timestep;
        }

        for perturber in &mut *perturbers {
            perturber.velocity += self.timestep * perturber.acceleration;
            perturber.position += perturber.velocity * 0.5 * self.timestep;
            perturber.epoch += 0.5 * self.timestep;
        }

    }

    fn timestep(&self) -> f64 {
        self.timestep
    }

    fn set_timestep(&mut self, timestep: f64) {
        self.timestep = timestep;
    }
}
