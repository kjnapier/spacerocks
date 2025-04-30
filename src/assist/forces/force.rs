use crate::spacerock::SpaceRock;
use nalgebra::Vector3;
use std::fmt::Debug;

use crate::assist::{SimulationParticle, SimulationState};

/// A force that can act on spacerocks in an N-body simulation.
/// 
/// This trait represents any force that can affect the motion of spacerocks,
/// such as gravity, radiation pressure, or non-gravitational forces.
/// Implementors must be thread-safe (Send + Sync) and clonable.
pub trait Force: Send + Sync + ForceClone {
    /// Calculate the acceleration of a set of spacerocks due to some force.
    ///
    /// # Arguments
    ///
    /// * `entities` - A mutable reference to a vector of spacerocks.
    ///
    /// # Returns
    ///
    /// * A vector of accelerations for each spacerock.
    fn calculate_acceleration(&self, state: &mut SimulationState) -> Vec<Vector3<f64>>;
    fn apply_acceleration(&self, state: &mut SimulationState);
    fn apply_acceleration_and_stm(&self, state: &mut SimulationState);
    
}


pub trait ForceClone {
    fn clone_box(&self) -> Box<dyn Force + Send + Sync>;
}

impl<T> ForceClone for T
where
    T: 'static + Force + Clone,
{
    fn clone_box(&self) -> Box<dyn Force + Send + Sync> {
        Box::new(self.clone())
    }
}

// We can now implement Clone manually by forwarding to clone_box.
impl Clone for Box<dyn Force + Send + Sync>{
    fn clone(&self) -> Box<dyn Force + Send + Sync> {
        self.clone_box()
    }
}

impl Debug for Box<dyn Force + Send + Sync> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Box<dyn Force>")
            .field("type", &std::any::type_name::<&dyn Force>())
            .finish()
    }
}
