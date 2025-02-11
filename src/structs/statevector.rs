use nalgebra::Vector3;

/// Represents position and velocity in Cartesian coordinates.
///
/// This structure holds 3D vectors for both position and velocity components
/// of a celestial body. Positions are in astronomical units (AU) and velocities
/// are in AU/day.
///
/// # Components
///
/// * `position`: 3D position vector [x, y, z] in AU
/// * `velocity`: 3D velocity vector [vx, vy, vz] in AU/day
///
/// # Examples
///
/// ```
/// use spacerocks::structs::StateVector;
/// use nalgebra::Vector3;
///
/// let state = StateVector::new(
///     Vector3::new(1.0, 0.0, 0.0),  // position (AU)
///     Vector3::new(0.0, 1.0, 0.0)   // velocity (AU/day)
/// );
///
/// println!("X position: {} AU", state.x());
/// println!("Y velocity: {} AU/day", state.vy());
/// ```
pub struct StateVector {
    pub position: Vector3<f64>,
    pub velocity: Vector3<f64>,
}

impl StateVector {
    /// Creates a new StateVector from position and velocity vectors.
    ///
    /// # Arguments
    ///
    /// * `position` - 3D position vector [x, y, z] in AU
    /// * `velocity` - 3D velocity vector [vx, vy, vz] in AU/day
    pub fn new(position: Vector3<f64>, velocity: Vector3<f64>) -> StateVector {
        StateVector { position, velocity }
    }

    /// Gets the X position component.
    pub fn x(&self) -> f64 {
        self.position.x
    }

    /// Gets the Y position component.
    pub fn y(&self) -> f64 {
        self.position.y
    }

    /// Gets the Z position component.
    pub fn z(&self) -> f64 {
        self.position.z
    }

    /// Gets the X velocity component.
    pub fn vx(&self) -> f64 {
        self.velocity.x
    }

    /// Gets the Y velocity component.
    pub fn vy(&self) -> f64 {
        self.velocity.y
    }

    /// Gets the Z velocity component.
    pub fn vz(&self) -> f64 {
        self.velocity.z
    }
}
