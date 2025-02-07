use crate::SpaceRock;
use crate::time::Time;
use crate::nbody::integrators::Integrator;
use crate::nbody::forces::Force;

use nalgebra::Vector3;

// use rayon::prelude::*;

/// A leapfrog-style integrator for N-body simulations.
/// 
/// The leapfrog integrator advances positions
/// and velocities using a three-stage pattern:
/// 1. Advances positions by half a timestep ("drift")
/// 2. Updates velocities using accelerations over a full timestep ("kick") 
/// 3. Completes the timestep with another half-step position update ("drift")
///
/// This approach provides good stability for orbital dynamics while being
/// simple to implement and computationally efficient.
#[derive(PartialEq, Debug, Clone, Copy)]
pub struct Leapfrog {
    /// Current timestep in simulation time units
    pub timestep: f64,
}

impl Leapfrog {
    /// Creates a new Leapfrog integrator with the specified timestep.
    ///
    /// # Arguments
    ///
    /// * `timestep` - Fixed timestep to use for integration
    pub fn new(timestep: f64) -> Leapfrog {
        Leapfrog { timestep }
    }
}

impl Integrator for Leapfrog {
    /// Advances the system one timestep using the drift-kick-drift sequence.
    ///
    /// The sequence:
    /// 1. Half-step position update
    /// 2. Full-step velocity update using calculated accelerations
    /// 3. Final half-step position update
    fn step(&mut self, particles: &mut Vec<SpaceRock>, epoch: &mut Time, forces: &Vec<Box<dyn Force + Send + Sync>>) {
        // drift
        for particle in &mut *particles {
            particle.position += particle.velocity * 0.5 * self.timestep;
            particle.epoch += 0.5 * self.timestep;
        }
      
        let mut accelerations = vec![Vector3::new(0.0, 0.0, 0.0); particles.len()];
        for force in forces {
            let acc = force.calculate_acceleration(particles);
            for (idx, a) in acc.iter().enumerate() {
                accelerations[idx] += a;
            }
        }

        // for particle in &mut *particles {
        //     particle.velocity += self.timestep * particle.acceleration;
        //     particle.position += particle.velocity * 0.5 * self.timestep;
        //     particle.epoch += 0.5 * self.timestep;
        // }

        *epoch += self.timestep;

        for (particle, acceleration) in particles.iter_mut().zip(accelerations.iter()) {
            particle.velocity += self.timestep * acceleration;
            particle.position += particle.velocity * 0.5 * self.timestep;
            particle.epoch = epoch.clone();
        }

        

    }

    fn timestep(&self) -> f64 {
        self.timestep
    }

    fn set_timestep(&mut self, timestep: f64) {
        self.timestep = timestep;
    }
}
