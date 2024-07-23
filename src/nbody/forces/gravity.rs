use crate::nbody::forces::Force;
use crate::spacerock::SpaceRock;
use crate::constants::GRAVITATIONAL_CONSTANT;

use rayon::prelude::*;
use nalgebra::Vector3;

use std::sync::Mutex;

// use std::sync::atomic::{AtomicU64, Ordering};

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


// impl NewtonianGravity {
//     const PARALLEL_THRESHOLD: usize = 1000;  // Adjust based on benchmarking

//     pub fn apply(&self, entities: &mut Vec<SpaceRock>) {
//         if entities.len() >= Self::PARALLEL_THRESHOLD {
//             self.apply_parallel(entities);
//         } else {
//             self.apply_serial(entities);
//         }
//     }

//     fn apply_serial(&self, entities: &mut Vec<SpaceRock>) {
//         let n_entities = entities.len();
//         for idx in 0..n_entities {
//             for jdx in (idx + 1)..n_entities {

//                 if (entities[idx].mass == 0.0) & (entities[jdx].mass == 0.0) {
//                     continue;
//                 }

//                 let r_vec = entities[idx].position - entities[jdx].position;
//                 let r = r_vec.norm();

//                 let xi = -GRAVITATIONAL_CONSTANT * r_vec / (r * r * r);
//                 let idx_acceleration = xi * entities[jdx].mass;
//                 let jdx_acceleration = -xi * entities[idx].mass;
//                 entities[idx].acceleration += idx_acceleration;
//                 entities[jdx].acceleration += jdx_acceleration;
//             }
//         }
//     }

//     fn apply_parallel(&self, entities: &mut Vec<SpaceRock>) {
//         let n_entities = entities.len();
//         let inv_masses: Vec<f64> = entities.iter().map(|e| 1.0 / e.mass).collect();
        
//         // Use AtomicU64 for thread-safe updates
//         let accelerations: Vec<[AtomicU64; 3]> = (0..n_entities)
//             .map(|_| [AtomicU64::new(0), AtomicU64::new(0), AtomicU64::new(0)])
//             .collect();

//         accelerations.par_iter().enumerate().for_each(|(i, acc_i)| {
//             let pos_i = &entities[i].position;
//             let mass_i = entities[i].mass;

//             for j in (i + 1)..n_entities {
//                 let pos_j = &entities[j].position;
//                 let mass_j = entities[j].mass;

//                 let dx = pos_i.x - pos_j.x;
//                 let dy = pos_i.y - pos_j.y;
//                 let dz = pos_i.z - pos_j.z;

//                 let dist_sq = dx * dx + dy * dy + dz * dz;
                
//                 if dist_sq > 0.0 {
//                     let inv_dist = 1.0 / dist_sq.sqrt();
//                     let inv_dist_cube = inv_dist * inv_dist * inv_dist;

//                     let force = GRAVITATIONAL_CONSTANT * mass_i * mass_j * inv_dist_cube;

//                     let fx = force * dx;
//                     let fy = force * dy;
//                     let fz = force * dz;

//                     atomic_add(&acc_i[0], -fx * inv_masses[i]);
//                     atomic_add(&acc_i[1], -fy * inv_masses[i]);
//                     atomic_add(&acc_i[2], -fz * inv_masses[i]);

//                     atomic_add(&accelerations[j][0], fx * inv_masses[j]);
//                     atomic_add(&accelerations[j][1], fy * inv_masses[j]);
//                     atomic_add(&accelerations[j][2], fz * inv_masses[j]);
//                 }
//             }
//         });

//         // Update entity accelerations
//         entities.par_iter_mut().zip(accelerations.par_iter()).for_each(|(entity, acc)| {
//             entity.acceleration.x += f64::from_bits(acc[0].load(Ordering::Relaxed));
//             entity.acceleration.y += f64::from_bits(acc[1].load(Ordering::Relaxed));
//             entity.acceleration.z += f64::from_bits(acc[2].load(Ordering::Relaxed));
//         });
//     }
// }

// impl Force for NewtonianGravity {
//     fn apply(&self, entities: &mut Vec<SpaceRock>) {
//         NewtonianGravity::apply(self, entities);
//     }
// }

// // Atomic addition for f64 using AtomicU64
// fn atomic_add(atomic: &AtomicU64, val: f64) {
//     let mut old_bits = atomic.load(Ordering::Relaxed);
//     loop {
//         let old_val = f64::from_bits(old_bits);
//         let new_val = old_val + val;
//         let new_bits = new_val.to_bits();
//         match atomic.compare_exchange_weak(
//             old_bits,
//             new_bits,
//             Ordering::Relaxed,
//             Ordering::Relaxed,
//         ) {
//             Ok(_) => break,
//             Err(x) => old_bits = x,
//         }
//     }
// }