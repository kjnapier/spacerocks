use crate::spacerock::SpaceRock;
use nalgebra::Vector3;

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
    fn calculate_acceleration(&self, entities: &mut Vec<SpaceRock>) -> Vec<Vector3<f64>>;

    /// Add the acceleration of each spacerock due to this force to `acc`, which has one entry
    /// per spacerock. Integrators call this; the default goes through
    /// [`Force::calculate_acceleration`], and forces can override it to skip the allocation.
    fn add_acceleration(&self, entities: &mut Vec<SpaceRock>, acc: &mut [Vector3<f64>]) {
        for (a, da) in acc.iter_mut().zip(self.calculate_acceleration(entities)) {
            *a += da;
        }
    }

    /// Whether this force is plain pairwise Newtonian gravity
    /// ([`NewtonianGravity`](super::NewtonianGravity)). The Wisdom–Holman and TRACE kicks then
    /// compute it from their own heliocentric coordinates instead of calling the force.
    fn is_newtonian_gravity(&self) -> bool {
        false
    }
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
