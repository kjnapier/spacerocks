use crate::{Time, Origin};

#[derive(Debug)]
/// Errors occurring during a simulation.
pub enum SimulationError {
    /// The origin of the particle did not match the simulation origin.
    OriginMismatch(Origin, Origin, String),
    /// The epoch of the particle did not match the simulation epoch.
    EpochMismatch(Time, Time, String),
    /// The particle was not found in the simulation.
    ParticleNotFound(String),
}

impl std::fmt::Display for SimulationError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SimulationError::OriginMismatch(op, osim, p) => write!(f, "The origin of the particle {} ({}) did not match the simulation origin ({}), and was not found in perturbers.", p, op, osim),
            SimulationError::EpochMismatch(tp, tsim, p) => write!(f, "The epoch of particle {} ({:?}) did not match the simulation epoch ({:?}).", p, tp, tsim),
            SimulationError::ParticleNotFound(p) => write!(f, "The particle {} was not found in the simulation.", p),
        }
    }
}

impl std::error::Error for SimulationError {}