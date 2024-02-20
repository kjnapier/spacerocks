
use crate::SpaceRock;
use crate::nbody::acceleration::{calculate_and_set_test_particle_accelerations, calculate_and_set_perturber_accelerations};

use rayon::prelude::*;

#[derive(PartialEq, Debug, Clone)]
pub struct Freezer {
    pub timestep: f64,
}

impl Freezer {
    pub fn new(timestep: f64) -> Freezer {
        Freezer { timestep }
    }

    pub fn step(&self, perturbers: &mut Vec<SpaceRock>, particles: &mut Vec<SpaceRock>) {
        for particle in &mut *particles {
            particle.epoch += self.timestep;
        }

        for mut perturber in &mut *perturbers {
            perturber.epoch += self.timestep;
        }
    }

}