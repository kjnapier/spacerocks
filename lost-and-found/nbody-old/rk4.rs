use crate::SpaceRock;
// use crate::nbody::acceleration::{calculate_and_set_test_particle_accelerations, calculate_and_set_perturber_accelerations};

use crate::Simulation;
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

    pub fn step(&self, sim: &mut Simulation) {
        0.0;
    }
       

}