use crate::orbit_type::OrbitType;

/// Calculate the eccentric anomaly from the true anomaly.
///
/// # Arguments
///
/// * `eccentricity` - Eccentricity of the orbit.
/// * `true_anomaly` - True anomaly of the orbit in radians.
///
/// # Returns
///
/// * `Result<f64, OrbitError>` - The conic anomaly in radians.
///
/// # Example
///
/// ```
/// use spacerocks::transforms;
/// let e = 0.5;
/// let true_anomaly = 0.5;
/// let result = transforms::calc_conic_anomaly_from_true_anomaly(e, true_anomaly);
/// ```
pub fn calc_conic_anomaly_from_true_anomaly(e: f64, true_anomaly: f64) -> Result<f64, Box<dyn std::error::Error>> {

    let orbit_type = OrbitType::from_eccentricity(e, 1e-10)?; // returns OrbitError error if e < 0.0
    match orbit_type {
        OrbitType::Circular => Ok(true_anomaly),
        // Half-angle form with atan2: well defined at f = pi, and E stays on the same
        // revolution as f (f in [0, 2pi) gives E in [0, 2pi)).
        OrbitType::Elliptical => Ok(2.0 * ((1.0 - e).sqrt() * (true_anomaly / 2.0).sin()).atan2((1.0 + e).sqrt() * (true_anomaly / 2.0).cos())),
        OrbitType::Parabolic => Ok((true_anomaly / 2.0).tan()),
        OrbitType::Hyperbolic => Ok(2.0 * (((e - 1.0) / (e + 1.0)).sqrt() * (true_anomaly / 2.0).tan()).atanh()),
        OrbitType::Radial => unreachable!(),
    }
}