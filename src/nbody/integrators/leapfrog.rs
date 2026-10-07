use crate::nbody::integrators::Integrator;
use crate::nbody::forces::{total_acceleration, Force};
use crate::state::{from_pv, pv, State};

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
    fn step(&mut self, states: &mut [State], masses: &[f64], t: &mut f64, forces: &[Box<dyn Force + Send + Sync>]) {
        // drift
        for s in states.iter_mut() {
            let (x, v) = pv(s);
            *s = from_pv(&(x + v * 0.5 * self.timestep), &v);
        }

        let mut accelerations = Vec::new();
        total_acceleration(forces, states, masses, &mut accelerations);

        *t += self.timestep;

        for (s, acceleration) in states.iter_mut().zip(accelerations.iter()) {
            let (x, v) = pv(s);
            let v = v + self.timestep * acceleration;
            *s = from_pv(&(x + v * 0.5 * self.timestep), &v);
        }
    }

    fn timestep(&self) -> f64 {
        self.timestep
    }

    fn set_timestep(&mut self, timestep: f64) {
        self.timestep = timestep;
    }

    fn fixed_timestep(&self) -> bool {
        true
    }
}
