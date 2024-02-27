use crate::nbody::forces::Force;
use crate::spacerock::SpaceRock;
use crate::constants::GRAVITATIONAL_CONSTANT;

use rayon::prelude::*;
use nalgebra::Vector3;

use std::sync::Mutex;

#[derive(Debug, Clone, Copy)]
pub struct NewtonianGravity;

impl Force for NewtonianGravity {
   
    fn apply(&self, entities: &mut Vec<SpaceRock>) {
        // Naive implementation of Newtonian gravity. O(0.5 * n^2) complexity.
        // Speed it up if you want!

        let n_entities = entities.len();
        for idx in 0..n_entities {
            for jdx in (idx + 1)..n_entities {

                if (entities[idx].mass == 0.0) & (entities[jdx].mass == 0.0) {
                    continue;
                }

                let r_vec = entities[idx].position - entities[jdx].position;
                let r = r_vec.norm();

                let xi = -GRAVITATIONAL_CONSTANT * r_vec / (r * r * r);
                let idx_acceleration = xi * entities[jdx].mass;
                let jdx_acceleration = -xi * entities[idx].mass;
                entities[idx].acceleration += idx_acceleration;
                entities[jdx].acceleration += jdx_acceleration;
            }
        }
    }
}