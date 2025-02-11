//! Orbital path types for celestial bodies.
use crate::errors::OrbitError;


/// Represents different types of orbital paths that a celestial body can follow.
/// 
/// The orbit type is determined by the eccentricity (e) of the orbit:
/// * `Hyperbolic`: e > 1
/// * `Parabolic`: e ≈ 1
/// * `Elliptical`: 0 < e < 1
/// * `Circular`: e ≈ 0
/// * `Radial`: Special case for radial trajectories
#[derive(Debug, PartialEq, Clone)]
pub enum OrbitType {
    Hyperbolic,
    Parabolic,
    Elliptical,
    Circular,
    Radial,
}

impl OrbitType {
    /// Determines the orbit type based on the eccentricity value
    ///
    /// # Arguments
    /// * `e` - Eccentricity of the orbit
    /// * `threshold` - Numerical threshold for considering an orbit circular (near 0) or parabolic (near 1)
    ///
    /// # Returns
    /// * `Result<OrbitType, OrbitError>` - The determined orbit type or an error for invalid eccentricity
    ///
    /// # Errors
    /// Returns `OrbitError::NegativeEccentricity` if the eccentricity is negative
    pub fn from_eccentricity(e: f64, threshold: f64) -> Result<OrbitType, OrbitError> {
        match e {
            e if e < 0.0 => Err(OrbitError::NegativeEccentricity(e)),
            e if e < threshold => Ok(OrbitType::Circular),
            e if e < 1.0 => Ok(OrbitType::Elliptical),
            e if (e - 1.0).abs() < threshold => Ok(OrbitType::Parabolic),
            _ => Ok(OrbitType::Hyperbolic),
        }
    }
}
