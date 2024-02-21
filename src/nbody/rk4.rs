use crate::SpaceRock;
use crate::nbody::acceleration::{calculate_and_set_test_particle_accelerations, calculate_and_set_perturber_accelerations};

use nalgebra::Vector3;

#[derive(PartialEq, Debug, Clone, Copy)]
struct K {
    pub velocity: Vector3<f64>,
    pub acceleration: Vector3<f64>,
}

#[derive(PartialEq, Debug, Clone)]
pub struct RK4 {
    pub timestep: f64,

}

impl RK4 {
    pub fn new(timestep: f64) -> RK4 {
        RK4 { timestep }
    }

    pub fn step(&self, perturbers: &mut Vec<SpaceRock>, particles: &mut Vec<SpaceRock>) {
        
        // let initial_particles = particles.clone();
        // let initial_perturbers = perturbers.clone();

        // let k1_particles = initial_particles.clone();
        // let k1_perturbers = initial_perturbers.clone();

        // calculate_and_set_test_particle_accelerations(&k1_perturbers, &mut particles.clone());
        // calculate_and_set_perturber_accelerations(&mut perturbers.clone());

        0.0;
    }
       

}