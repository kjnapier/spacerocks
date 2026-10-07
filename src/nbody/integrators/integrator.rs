use crate::nbody::forces::Force;
use crate::state::State;


/// A numerical integrator for advancing a system of particles forward in time.
///
/// Integrators work on plain state: a slice of states `[x, y, z, vx, vy, vz]` (AU, AU/day),
/// updated in place, the bodies' masses (solar masses, massive bodies first, 0 for test
/// particles), and the time `t` in days (the TDB Julian date in a
/// [`Simulation`](crate::nbody::Simulation)). Names, epochs, origins and reference planes are
/// the simulation's business. Implementors must be thread-safe (Send + Sync) and clonable.
pub trait Integrator: Send + Sync + IntegratorClone {
    /// Advances the system one timestep forward
    ///
    /// # Arguments
    ///
    /// * `states` - States of the bodies, updated in place
    /// * `masses` - Masses of the bodies
    /// * `t` - Current time, advanced in place
    /// * `forces` - Forces acting on the system
    fn step(&mut self, states: &mut [State], masses: &[f64], t: &mut f64, forces: &[Box<dyn Force + Send + Sync>]);

    /// Advances the system `n` timesteps. The default calls [`Integrator::step`] `n` times;
    /// integrators can override it to share work between steps, with the same result.
    ///
    /// # Arguments
    ///
    /// * `states` - States of the bodies, updated in place
    /// * `masses` - Masses of the bodies
    /// * `t` - Current time, advanced in place
    /// * `forces` - Forces acting on the system
    /// * `n` - Number of steps
    fn steps(&mut self, states: &mut [State], masses: &[f64], t: &mut f64, forces: &[Box<dyn Force + Send + Sync>], n: usize) {
        for _ in 0..n {
            self.step(states, masses, t, forces);
        }
    }

    /// Returns the current timestep of the integrator
    fn timestep(&self) -> f64;

    /// Sets the timestep of the integrator
    ///
    /// # Arguments
    ///
    /// * `timestep` - New timestep value to use
    fn set_timestep(&mut self, timestep: f64);

    /// Whether every step is exactly [`Integrator::timestep`] long (the integrator never
    /// changes it by itself). [`Simulation::integrate`](crate::nbody::Simulation::integrate)
    /// then hands all the full steps to [`Integrator::steps`] at once.
    fn fixed_timestep(&self) -> bool {
        false
    }

    /// Whether [`Integrator::interpolate`] can return states between steps
    fn has_dense_output(&self) -> bool {
        false
    }

    /// States of every body at time `t` (on the same clock as [`Integrator::step`]'s), if `t`
    /// falls within the last completed step. Integrators without dense output return `None`.
    fn interpolate(&self, _t: f64) -> Option<Vec<State>> {
        None
    }
}


pub trait IntegratorClone {
    fn clone_box(&self) -> Box<dyn Integrator + Send + Sync>;
}

impl<T> IntegratorClone for T
where
    T: 'static + Integrator + Clone,
{
    fn clone_box(&self) -> Box<dyn Integrator + Send + Sync> {
        Box::new(self.clone())
    }
}

// We can now implement Clone manually by forwarding to clone_box.
impl Clone for Box<dyn Integrator + Send + Sync>{
    fn clone(&self) -> Box<dyn Integrator + Send + Sync> {
        self.clone_box()
    }
}
