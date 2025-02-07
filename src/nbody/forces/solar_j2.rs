use crate::nbody::forces::Force;

use crate::spacerock::SpaceRock;
use crate::constants::{GRAVITATIONAL_CONSTANT, SPEED_OF_LIGHT};

// use rayon::prelude::*;
use nalgebra::Vector3;


/// Implements perturbations due to the Sun's oblateness (J2).
/// Accounts for the Sun's slight equatorial bulge which creates
/// a non-spherical component to its gravitational field.
#[derive(Debug, Clone, Copy)]
pub struct SolarJ2;

const sun_j2: f64 = 2.17e-7;
const sun_radius: f64 = 696_342.0 / 149_597_870.7;

impl Force for SolarJ2 {
    /// Calculates acceleration from solar J2 effect. The force depends on 
    /// latitude with respect to the solar equator and falls off as r⁻⁵.
    ///
    /// # Arguments
    /// * `entities` - Vector of bodies (must include "sun")
    ///
    /// # Returns
    /// * Vector of J2 perturbation accelerations
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
            let z2_r2 = r_vec[2] * r_vec[2] / (r * r);

            let factor = 2.0 * sun_j2 * mu * sun_radius * sun_radius / (2.0 * r.powi(5));
            let ax = factor * r_vec.x * (5.0 * z2_r2 - 1.0);
            let ay = factor * r_vec.y * (5.0 * z2_r2 - 1.0);
            let az = factor * r_vec.z * (5.0 * z2_r2 - 3.0);
        
            // entity.acceleration += xi;
            acceleration[idx] = Vector3::new(ax, ay, az);
        }
        acceleration
    }
}