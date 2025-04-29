use crate::spacerock::SpaceRock;
use crate::constants::GRAVITATIONAL_CONSTANT;

use crate::assist::forces::Force;
use crate::assist::{SimulationParticle, SimulationState};

use nalgebra::Vector3;


/// Implementation of classical Newtonian gravitational force.
///
/// Calculates gravitational interactions between all pairs of bodies using Newton's
/// law of universal gravitation. 
#[derive(Debug, Clone, Copy)]
pub struct NewtonianGravity;

impl Force for NewtonianGravity {
    /// Calculates gravitational accelerations for a system of bodies using Newton's law
    /// of universal gravitation.
    fn calculate_acceleration(&self, state: &mut SimulationState) -> Vec<Vector3<f64>> {
        // Naive implementation of Newtonian gravity. O(0.5 * n^2) complexity.
        // Speed it up if you want!

        let spice_particles = &state.spice_particles;
        let particles = &state.particles_1;

        // we only need to calculate the accelerations for the particles. assume for now that they're all massless.

        let mut acceleration = vec![Vector3::zeros(); particles.len()];

        let n_particles = particles.len();
        let n_spice_particles = spice_particles.len();

        for idx in 0..n_particles {

            let particle = &particles[idx];

            for jdx in 0..n_spice_particles {

                let p2 = &spice_particles[jdx];

                let r_vec = particle.position - p2.position;
                let r = r_vec.norm();

                let xi = -GRAVITATIONAL_CONSTANT * r_vec / (r * r * r);
                let idx_acceleration = xi * p2.mass;
                acceleration[idx] += idx_acceleration;
            }
        }
        acceleration
    }

    fn apply_acceleration(&self, state: &mut SimulationState) {
        let spice_particles = &state.spice_particles;
        let particles = &mut state.particles_1;

        let n_particles = particles.len();
        let n_spice_particles = spice_particles.len();

        for idx in 0..n_particles {

            let particle = &mut particles[idx];

            for jdx in 0..n_spice_particles {

                let p2 = &spice_particles[jdx];

                let r_vec = particle.position - p2.position;
                let r = r_vec.norm();

                let xi = -GRAVITATIONAL_CONSTANT * r_vec / (r * r * r);
                let idx_acceleration = xi * p2.mass;
                particle.acceleration += idx_acceleration;
            }
        }
        
    }

}






// let n_entities = entities.len();
        // for idx in 0..n_entities {

        //     let idx_massless = entities[idx].mass() == 0.0;

        //     for jdx in (idx + 1)..n_entities {

        //         // if (entities[idx].mass() == 0.0) & (entities[jdx].mass() == 0.0) {
        //         //     continue;
        //         // }

        //         if idx_massless & (entities[jdx].mass() == 0.0) {
        //             continue;
        //         }

        //         let r_vec = entities[idx].position - entities[jdx].position;
        //         let r = r_vec.norm();

        //         let xi = -GRAVITATIONAL_CONSTANT * r_vec / (r * r * r);
        //         let idx_acceleration = xi * entities[jdx].mass();
        //         let jdx_acceleration = -xi * entities[idx].mass();
        //         acceleration[idx] += idx_acceleration;
        //         acceleration[jdx] += jdx_acceleration;
        //     }
        // }



//     let mut massless_indices = Vec::new();
    //     let mut massive_indices = Vec::new();
    //     for (idx, entity) in entities.iter().enumerate() {
    //         if entity.mass() == 0.0 {
    //             massless_indices.push(idx);
    //         } else {
    //             massive_indices.push(idx);
    //         }
    //     }

    //     let n_massive = massive_indices.len();

    //     for ii in 0..n_massive {
    //         let idx = massive_indices[ii];

    //         // loop over all massive entities
    //         for jj in (ii + 1)..n_massive {
    //             let jdx = massive_indices[jj];
    //             let r_vec = entities[idx].position - entities[jdx].position;
    //             let r = r_vec.norm();
    //             let xi = -GRAVITATIONAL_CONSTANT * r_vec / (r * r * r);
    //             let idx_acceleration = xi * entities[jdx].mass();
    //             let jdx_acceleration = -xi * entities[idx].mass();
    //             acceleration[idx] += idx_acceleration;
    //             acceleration[jdx] += jdx_acceleration;
    //         }

    //         // loop over all massless entities
    //         for kdx in &massless_indices {
    //             let r_vec = entities[idx].position - entities[*kdx].position;
    //             let r = r_vec.norm();
    //             let xi = -GRAVITATIONAL_CONSTANT * r_vec / (r * r * r);
    //             let kdx_acceleration = -xi * entities[idx].mass();
    //             acceleration[*kdx] += kdx_acceleration;
    //         }
    //     }