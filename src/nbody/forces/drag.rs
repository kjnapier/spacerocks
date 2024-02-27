use crate::nbody::forces::Force;

use crate::spacerock::SpaceRock;
use crate::constants::GRAVITATIONAL_CONSTANT;

use rayon::prelude::*;
use nalgebra::Vector3;

/*
An ad-hoc drag force that is meant to simulate the effects of a fluid on the motion of the particles.

Only used for testing purposes.
*/

#[derive(Debug, Clone, Copy)]
pub struct Drag;

impl Force for Drag {
    // fn apply(&self, particles: &mut Vec<SpaceRock>, perturbers: &mut Vec<SpaceRock>) {

    //     let mut system_mass = 0.0;
    //     for perturber in perturbers.iter() {
    //         system_mass += perturber.mass;
    //     }

    //     let mu = GRAVITATIONAL_CONSTANT * system_mass;

    //     for particle in particles.iter_mut() {
    //         let v = particle.velocity.norm();
    //         let v_hat = particle.velocity / v;

    //         // the fluid is moving in a circular orbit at radius r
    //         let r = particle.position.norm();
    //         let fluid_speed = (mu / r).sqrt();
    //         let fluid_velocity_hat = Vector3::new(-particle.position.y, particle.position.x, 0.0) / r;
    //         let fluid_velocity = fluid_speed * fluid_velocity_hat;

    //         let v_rel_vec = particle.velocity - fluid_velocity;
    //         let v_rel = v_rel_vec.norm();
    //         let v_rel_hat = v_rel_vec / v_rel;


    //         particle.acceleration -= 0.01 * v_rel * v_rel * v_rel_hat;
    //     }
    // }

    fn apply(&self, entities: &mut Vec<SpaceRock>) {
        entities.par_iter_mut().for_each(|entity| {
            let v = entity.velocity.norm();
            let v_hat = entity.velocity / v;

            // the fluid is moving in a circular orbit at radius r
            let r = entity.position.norm();
            let fluid_speed = (GRAVITATIONAL_CONSTANT * 1.0 / r).sqrt();
            let fluid_velocity_hat = Vector3::new(-entity.position.y, entity.position.x, 0.0) / r;
            let fluid_velocity = fluid_speed * fluid_velocity_hat;

            let v_rel_vec = entity.velocity - fluid_velocity;
            let v_rel = v_rel_vec.norm();
            let v_rel_hat = v_rel_vec / v_rel;

            entity.acceleration -= 0.01 * v_rel * v_rel * v_rel_hat;
        });
    }

}