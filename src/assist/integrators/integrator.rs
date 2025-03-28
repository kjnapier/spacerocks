use crate::SpaceRock;
use crate::time::Time;
use crate::assist::forces::Force;
use crate::assist::SimulationState;


/// A numerical integrator for advancing a system of particles forward in time.
/// 
/// This trait defines the core functionality required for any numerical integrator
/// in the system. Implementors must be thread-safe (Send + Sync) and clonable.
/// 
/// The integrator is responsible for:
/// - Advancing particle states (positions and velocities) through time
/// - Managing timestep size
/// - Coordinating force calculations
pub trait Integrator: Send + Sync + IntegratorClone {
    /// Advances the system one timestep forward
    ///
    /// # Arguments
    ///
    /// * `particles` - Vector of particles to be integrated
    /// * `epoch` - Current simulation time, updated in-place to the new time
    /// * `forces` - Vector of forces acting on the system
    fn step(&mut self, state: &mut SimulationState, forces: &Vec<Box<dyn Force + Send + Sync>>);

    /// Returns the current timestep of the integrator
    fn timestep(&self) -> f64;

    /// Sets the timestep of the integrator
    ///
    /// # Arguments
    ///
    /// * `timestep` - New timestep value to use
    fn set_timestep(&mut self, timestep: f64);
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
