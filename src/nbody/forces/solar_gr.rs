use crate::nbody::forces::Force;

use crate::spacerock::SpaceRock;
use crate::constants::{GRAVITATIONAL_CONSTANT, SPEED_OF_LIGHT};

use rayon::prelude::*;
use nalgebra::Vector3;


#[derive(Debug, Clone, Copy)]
pub struct SolarGR;

impl Force for SolarGR {
    // fn apply(&self, particles: &mut Vec<SpaceRock>, perturbers: &mut Vec<SpaceRock>) {

    //     let sun_index = perturbers.iter().position(|x| x.name == "sun").unwrap();
    //     let sun = perturbers[sun_index].clone();
    //     let mu = GRAVITATIONAL_CONSTANT * sun.mass;

    //     for idx in 0..perturbers.len() {
    //         if idx == sun_index {
    //             continue;
    //         }

    //         let perturber = &mut perturbers[idx];

    //         let r_vec = perturber.position - sun.position;
    //         let r = r_vec.norm();

    //         let v_vec = perturber.velocity - sun.velocity;
    //         let v = v_vec.norm();

    //         let s0 = mu / (SPEED_OF_LIGHT.powi(2) * r * r * r);
    //         let s1 = ((4.0 * mu) / r - v * v) * r_vec;
    //         let s2 = 4.0 * (r_vec.dot(&v_vec)) * v_vec;

    //         let xi = s0 * (s1 + s2);
    //         let idx_acceleration = xi;
    //         perturber.acceleration += idx_acceleration;
    //     }

    //     for particle in particles.iter_mut() {
    //         let r_vec = particle.position - sun.position;
    //         let r = r_vec.norm();

    //         let v_vec = particle.velocity - sun.velocity;
    //         let v = v_vec.norm();

    //         let s0 = mu / (SPEED_OF_LIGHT.powi(2) * r * r * r);
    //         let s1 = ((4.0 * mu) / r - v * v) * r_vec;
    //         let s2 = 4.0 * (r_vec.dot(&v_vec)) * v_vec;

    //         let a = s0 * (s1 + s2);
    //         particle.acceleration += a;
    //     }

    // }

    fn apply(&self, entities: &mut Vec<SpaceRock>) {
        let sun_index = entities.iter().position(|x| x.name == "sun").unwrap();
        let sun = entities[sun_index].clone();
        let mu = GRAVITATIONAL_CONSTANT * sun.mass;

        entities.par_iter_mut().enumerate().for_each(|(idx, entity)| {
            if idx == sun_index {
                return;
            }

            let r_vec = entity.position - sun.position;
            let r = r_vec.norm();

            let v_vec = entity.velocity - sun.velocity;
            let v = v_vec.norm();

            let s0 = mu / (SPEED_OF_LIGHT.powi(2) * r * r * r);
            let s1 = ((4.0 * mu) / r - v * v) * r_vec;
            let s2 = 4.0 * (r_vec.dot(&v_vec)) * v_vec;

            let xi = s0 * (s1 + s2);
            entity.acceleration += xi;
        });
    }
}