use crate::nbody::forces::Force;

use crate::spacerock::SpaceRock;
use crate::constants::{GRAVITATIONAL_CONSTANT, SPEED_OF_LIGHT};

use rayon::prelude::*;
use nalgebra::Vector3;
use std::sync::Arc;


#[derive(Debug, Clone, Copy)]
pub struct SolarGR;

impl Force for SolarGR {
    
    fn apply(&self, entities: &mut Vec<SpaceRock>) {

        let sun_index = entities.iter().position(|x| x.name == Arc::new("sun".to_string())).unwrap();
        let sun = entities[sun_index].clone();
        let mu = GRAVITATIONAL_CONSTANT * sun.mass;

        for idx in 0..entities.len() {
            if idx == sun_index {
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
            entity.acceleration += xi;
        }
    }
}