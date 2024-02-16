use crate::constants::gravitational_constant;
use crate::spacerock::SpaceRock;
use nalgebra::Vector3;
use rayon::prelude::*;

pub fn calculate_and_set_test_particle_accelerations(perturbers: &Vec<SpaceRock>, particles: &mut Vec<SpaceRock>) {
    for particle in particles.iter_mut() {
        particle.acceleration = Vector3::new(0.0, 0.0, 0.0);
        for perturber in perturbers {
            let r_vec = particle.position - perturber.position;
            let r = r_vec.norm();
            let r_hat = r_vec / r;
            let a = -gravitational_constant * perturber.mass * r_hat / r.powi(2);
            particle.acceleration += a;
        }
    }
}

pub fn calculate_and_set_perturber_accelerations(perturbers: &mut Vec<SpaceRock>) {
    // First, zero out the accelerations
    for perturber in perturbers.iter_mut() {
        perturber.acceleration = Vector3::new(0.0, 0.0, 0.0);
    }

    // Then, calculate the accelerations.
    let N = perturbers.len();
    for idx in 0..N {
        for jdx in (idx + 1)..N {
            let r_vec = perturbers[idx].position - perturbers[jdx].position;
            let r = r_vec.norm();
            let xi = -gravitational_constant * r_vec / (r * r * r);
            let da_i = xi * perturbers[jdx].mass;
            let da_j = -xi * perturbers[idx].mass;

            perturbers[idx].acceleration += da_i;
            perturbers[jdx].acceleration += da_j;
        }
    }
}