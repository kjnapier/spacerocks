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

        // let positions: Vec<Vector3<f64>> = entities.iter().map(|entity| entity.position).collect();
        // let masses: Vec<f64> = entities.iter().map(|entity| entity.mass).collect();

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

    //     // calculate the acceleration for each entity in parallel
    //     let acceleration: Vec<Mutex<Vector3<f64>>> = entities.par_iter().map(|entity| {
    //         let acceleration = Mutex::new(Vector3::new(0.0, 0.0, 0.0));
    //         for other_entity in entities.iter() {
    //             if entity.name != other_entity.name {

    //                 if (entity.mass == 0.0) & (other_entity.mass == 0.0) {
    //                     continue;
    //                 }

    //                 let r_vec = entity.position - other_entity.position;
    //                 let r = r_vec.norm();
    //                 let xi = -GRAVITATIONAL_CONSTANT * r_vec / (r * r * r);
    //                 let mut acceleration = acceleration.lock().unwrap();
    //                 *acceleration += xi * other_entity.mass;
    //             }
    //         }
    //         acceleration
    //     }).collect();     

    //     // apply the acceleration to each entity
    //     entities.par_iter_mut().enumerate().for_each(|(idx, entity)| {
    //         let mut acceleration = acceleration[idx].lock().unwrap();
    //         entity.acceleration += *acceleration;
    //     });
    // }
}