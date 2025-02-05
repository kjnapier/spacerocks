#[derive(Debug, PartialEq)]
/// Errors that can occur when calculating orbits.
pub enum OrbitError {
    /// The eccentricity of the orbit is negative.
    NegativeEccentricity(f64),

    /// The iteration to solve Kepler's equation did not converge.
    ConvergenceFailure(f64, f64),  // (eccentricity, mean_anomaly)
}

impl std::fmt::Display for OrbitError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            OrbitError::NegativeEccentricity(e) => write!(f, "Eccentricity cannot be negative: {}", e),
            OrbitError::ConvergenceFailure(e, m) => write!(f, "Failed to converge for eccentricity {} and mean anomaly {}", e, m),
        }
    }
}

impl std::error::Error for OrbitError {}