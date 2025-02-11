use std::f64::consts::PI;

/// Represents an orbit using Keplerian orbital elements.
///
/// This structure uses the standard set of Keplerian orbital elements to describe
/// a celestial body's orbit. All angular quantities are in radians, and distances
/// are in astronomical units (AU).
///
/// # Orbital Elements
///
/// * `e`: Eccentricity (dimensionless)
/// * `q`: Perihelion distance (AU)
/// * `inc`: Inclination (rad)
/// * `node`: Longitude of ascending node (rad)
/// * `arg`: Argument of perihelion (rad)
/// * `true_anomaly`: True anomaly (rad)
///
/// # Examples
///
/// ```
/// use spacerocks::structs::KeplerOrbit;
///
/// let orbit = KeplerOrbit::new(
///     0.5,    // eccentricity
///     1.2,    // perihelion distance (AU)
///     0.3,    // inclination (rad)
///     1.5,    // longitude of ascending node (rad)
///     2.1,    // argument of perihelion (rad)
///     0.0,    // true anomaly (rad)
/// );
///
/// println!("Semi-major axis: {} AU", orbit.a());
/// println!("Longitude of perihelion: {} rad", orbit.varpi());
/// ```
pub struct KeplerOrbit {
    pub e: f64,
    pub q: f64,
    pub inc: f64,
    pub node: f64,
    pub arg: f64,
    pub true_anomaly: f64,
}

impl KeplerOrbit {

    /// Creates a new KeplerOrbit from orbital elements.
    ///
    /// # Arguments
    ///
    /// * `e` - Eccentricity (dimensionless)
    /// * `q` - Perihelion distance (AU)
    /// * `inc` - Inclination (rad)
    /// * `node` - Longitude of ascending node (rad)
    /// * `arg` - Argument of perihelion (rad)
    /// * `true_anomaly` - True anomaly (rad)
    pub fn new(e: f64, q: f64, inc: f64, node: f64, arg: f64, true_anomaly: f64) -> KeplerOrbit {
        KeplerOrbit { e, q, inc, node, arg, true_anomaly }
    }

    /// Calculates the semi-major axis.
    ///
    /// # Returns
    ///
    /// Semi-major axis in AU. 
    pub fn a(&self) -> f64 {
        self.q / (1.0 - self.e)
    }

    /// Calculates the longitude of perihelion (ϖ = Ω + ω).
    ///
    /// # Returns
    ///
    /// Longitude of perihelion in radians, normalized to [0, 2π].
    pub fn varpi(&self) -> f64 {
        (self.node + self.arg) % (2.0 * PI)
    }
    
}