use crate::nbody::forces::Force;

use crate::spacerock::SpaceRock;
use crate::constants::{GRAVITATIONAL_CONSTANT, SPEED_OF_LIGHT};

// use rayon::prelude::*;
use nalgebra::Vector3;


#[derive(Debug, Clone, Copy)]
/// Implementation of general relativistic corrections to solar gravity.
pub struct SolarGR;

impl Force for SolarGR {
    
    /// Calculates relativistic acceleration.
    fn calculate_acceleration(&self, entities: &mut Vec<SpaceRock>) -> Vec<Vector3<f64>> {

        let mut acceleration = vec![Vector3::zeros(); entities.len()];

        let sun_index = entities.iter().position(|x| x.name == *"sun").unwrap();
        let sun = entities[sun_index].clone();
        let mu = GRAVITATIONAL_CONSTANT * sun.mass();

        for idx in 0..entities.len() {
            if idx == sun_index {
                acceleration[idx] = Vector3::zeros();
                continue;
            }

            let entity = &mut entities[idx];

            let r_vec = entity.position - sun.position;
            let r = r_vec.norm();

            let v_vec = entity.velocity - sun.velocity;
            let v = v_vec.norm();

            let s0 = mu / (SPEED_OF_LIGHT.powi(2) * r * r * r);
            let s1 = ((4.0 * mu) / r - v * v) * r_vec;
            let s2 = 4.0 * (r_vec.dot(&v_vec)) * v_vec;

            let xi = s0 * (s1 + s2);
            // entity.acceleration += xi;
            acceleration[idx] = xi;
        }
        acceleration
    }
}