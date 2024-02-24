use crate::nbody::forces::Force;
use crate::spacerock::SpaceRock;
use crate::constants::GRAVITATIONAL_CONSTANT;

use rayon::prelude::*;
use nalgebra::Vector3;


pub struct NewtonianGravity;

impl Force for NewtonianGravity {
    fn apply(&self, particles: &mut Vec<SpaceRock>, perturbers: &mut Vec<SpaceRock>) {
        let n_perturbers = perturbers.len();
        for idx in 0..perturbers.len() {
            for jdx in (idx + 1)..perturbers.len() {
                let r_vec = perturbers[idx].position - perturbers[jdx].position;
                let r = r_vec.norm();

                let xi = -GRAVITATIONAL_CONSTANT * r_vec / (r * r * r);
                let idx_acceleration = xi * perturbers[jdx].mass;
                let jdx_acceleration = -xi * perturbers[idx].mass;
                perturbers[idx].acceleration += idx_acceleration;
                perturbers[jdx].acceleration += jdx_acceleration;
            }

            for particle in particles.iter_mut() {
                let r_vec = particle.position - perturbers[idx].position;
                let r = r_vec.norm();
                let xi = -GRAVITATIONAL_CONSTANT * r_vec / (r * r * r);
                let a = xi * perturbers[idx].mass;
                particle.acceleration += a;
            }
        }
    }
}