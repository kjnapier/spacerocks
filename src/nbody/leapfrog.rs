
use crate::spacerock::SpaceRock;
use crate::nbody::simparticle::SimParticle;
use crate::nbody::acceleration::{calculate_and_set_test_particle_accelerations, calculate_and_set_perturber_accelerations};

pub fn leapfrog_step(perturbers: &mut Vec<SimParticle>, test_particles: &mut Vec<SimParticle>, timestep: f64) {
    
    // drift
    for particle in &mut *test_particles {
        particle.position += particle.velocity * 0.5 * timestep;
    }

    // kick + drift
    calculate_and_set_test_particle_accelerations(perturbers, test_particles);
    for particle in &mut *test_particles {
        particle.velocity += timestep * particle.acceleration;
        particle.position += particle.velocity * 0.5 * timestep;
        particle.epoch += timestep;
    }


    // drift
    for mut perturber in &mut *perturbers {
        perturber.position += perturber.velocity * 0.5 * timestep;
    }

    // kick + drift
    calculate_and_set_perturber_accelerations(perturbers);
    for mut perturber in &mut *perturbers {
        perturber.velocity += timestep * perturber.acceleration;
        perturber.position += perturber.velocity * 0.5 * timestep;
        perturber.epoch += timestep;
    }

}