use crate::nbody::forces::Force;
use crate::spacerock::SpaceRock;
use crate::constants::GRAVITATIONAL_CONSTANT;

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
    fn calculate_acceleration(&self, entities: &mut Vec<SpaceRock>) -> Vec<Vector3<f64>> {
        let mut acceleration = vec![Vector3::zeros(); entities.len()];
        self.add_acceleration(entities, &mut acceleration);
        acceleration
    }

    fn add_acceleration(&self, entities: &mut Vec<SpaceRock>, acc: &mut [Vector3<f64>]) {
        // Naive implementation of Newtonian gravity. O(0.5 * n^2) complexity. Particles are
        // sorted by mass, so the loop stops at the first massless one: test particles don't
        // pull on each other.
        let n_entities = entities.len();
        for idx in 0..n_entities {
            let m_idx = entities[idx].mass();
            if m_idx == 0.0 {
                break;
            }
            let x_idx = entities[idx].position;
            let mut acc_idx = Vector3::zeros();
            for (jdx, other) in entities.iter().enumerate().skip(idx + 1) {
                let r_vec = x_idx - other.position;
                let r2 = r_vec.norm_squared();
                let xi = (GRAVITATIONAL_CONSTANT / (r2 * r2.sqrt())) * r_vec;
                acc_idx -= other.mass() * xi;
                acc[jdx] += m_idx * xi;
            }
            acc[idx] += acc_idx;
        }
    }

    fn is_newtonian_gravity(&self) -> bool {
        true
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