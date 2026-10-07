use crate::state::State;
use nalgebra::Vector3;

/// A force acting on the bodies of an N-body simulation.
///
/// Forces are functions of plain state: they see every body's state `[x, y, z, vx, vy, vz]`
/// and mass (solar masses, 0 for test particles), with the massive bodies first and the most
/// massive (the central body) at the start, and nothing else. Implementors must be
/// thread-safe (Send + Sync) and clonable.
pub trait Force: Send + Sync + ForceClone {
    /// Add the acceleration of each body due to this force to `acc`, which has one entry per
    /// body.
    fn add_acceleration(&self, states: &[State], masses: &[f64], acc: &mut [Vector3<f64>]);

    /// Whether this force is plain pairwise Newtonian gravity
    /// ([`NewtonianGravity`](super::NewtonianGravity)). The Wisdom–Holman and TRACE kicks then
    /// compute it from their own heliocentric coordinates instead of calling the force.
    fn is_newtonian_gravity(&self) -> bool {
        false
    }
}

/// The total acceleration of every body due to `forces`, written into `acc` (resized to one
/// entry per body).
pub fn total_acceleration(forces: &[Box<dyn Force + Send + Sync>], states: &[State], masses: &[f64], acc: &mut Vec<Vector3<f64>>) {
    acc.clear();
    acc.resize(states.len(), Vector3::zeros());
    for force in forces {
        force.add_acceleration(states, masses, acc);
    }
}

/// Index of the most massive body (the first, on ties), if any has mass: the central body of
/// the Wisdom–Holman and TRACE splittings and the "Sun" of [`SolarGR`](super::SolarGR) and
/// [`SolarJ2`](super::SolarJ2).
pub fn central_body(masses: &[f64]) -> Option<usize> {
    let best = masses.iter().enumerate().fold(None, |best: Option<(usize, f64)>, (i, &m)| match best {
        Some((_, mb)) if mb >= m => best,
        _ => Some((i, m)),
    });
    best.filter(|&(_, m)| m > 0.0).map(|(i, _)| i)
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
