//! Reference plane for specifying coordinate system.

use nalgebra::Matrix3;
use crate::constants::{ROTATION_J2000, ROTATION_ECLIPJ2000, ROTATION_INVARIABLE, ROTATION_GALACTIC, ROTATION_FK4};

use serde::{Serialize, Deserialize};

/// Defines the orientation of coordinate systems for dynamics calculations.
///
/// A reference plane specifies the fundamental plane and primary direction
/// that orient a coordinate system in space. Each plane provides a rotation
/// matrix that transforms coordinates to the J2000 equatorial system.
///
/// The choice of reference plane depends on the specific problem:
/// - J2000/FK4 
/// - ECLIPJ2000 
/// - GALACTIC 
/// - INVARIABLE 
#[derive(PartialEq, Debug, Clone, Serialize, Deserialize)]
#[derive(Default)]
pub enum ReferencePlane {
    /// Earth's mean equator and equinox at J2000.0 epoch (JD 2451545.0)
    J2000,
    /// Earth's mean ecliptic and equinox at J2000.0 epoch (default)
    #[default]
    ECLIPJ2000,
    /// Invariable plane of the Solar System
    /// 
    /// The plane perpendicular to the solar system's total angular momentum vector
    INVARIABLE,
    /// Galactic reference plane
    ///
    /// A coordinate system aligned with the structure of the Milky Way galaxy, where:
    /// - The fundamental plane (b = 0°) is aligned with the mean plane of the Milky Way's disk
    /// - The primary direction (l = 0°) points toward the galactic center 
    /// - The north galactic pole (b = +90°) points toward the galactic north pole
    GALACTIC,
    /// Earth's mean equator and equinox at B1950.0 epoch
    /// 
    /// Older reference system, mainly used for historical compatibility
    FK4,
}


impl ReferencePlane {

    /// Create a new ReferencePlane from a string.
    ///
    /// # Arguments
    /// * `s` - The string representation of the ReferencePlane (J2000, ECLIPJ2000, INVARIABLE, GALACTIC, or FK4).  
    ///
    /// # Example
    /// ```
    /// use spacerocks::coordinates::ReferencePlane;
    /// let reference_plane = ReferencePlane::from_str("J2000").unwrap();
    /// ```
    pub fn from_str(s: &str) -> Result<Self, String> {
        match s.to_uppercase().as_str() {
            "J2000" => Ok(ReferencePlane::J2000),
            "ECLIPJ2000" => Ok(ReferencePlane::ECLIPJ2000),
            "INVARIABLE" => Ok(ReferencePlane::INVARIABLE),
            "GALACTIC" => Ok(ReferencePlane::GALACTIC),
            "FK4" => Ok(ReferencePlane::FK4),
            _ => Err(format!("Invalid frame: {}", s))
        }
    }

    /// Return the rotation matrix of the ReferencePlane.
    /// These rotation matrices are used to transform the coordinates from the specified reference plane to the J2000 reference plane, 
    /// and can be found in the constants module.
    ///
    /// # Example
    /// ```
    /// let reference_plane = ReferencePlane::from_str("J2000").unwrap();
    /// let rotation_matrix = reference_plane.get_rotation_matrix();
    /// ```
    pub fn get_rotation_matrix(&self) -> Matrix3<f64> {
        match self {
            ReferencePlane::J2000 => ROTATION_J2000,
            ReferencePlane::ECLIPJ2000 => ROTATION_ECLIPJ2000,
            ReferencePlane::INVARIABLE => ROTATION_INVARIABLE,
            ReferencePlane::GALACTIC => ROTATION_GALACTIC,
            ReferencePlane::FK4 => ROTATION_FK4,
        }
    }

    /// Return the string representation of the ReferencePlane.
    ///
    /// # Example
    /// ```
    /// let reference_plane = ReferencePlane::from_str("J2000").unwrap();
    /// assert_eq!(reference_plane.as_str(), "J2000");
    /// ```
    pub fn as_str(&self) -> &str {
        match self {
            ReferencePlane::J2000 => "J2000",
            ReferencePlane::ECLIPJ2000 => "ECLIPJ2000",
            ReferencePlane::INVARIABLE => "INVARIABLE",
            ReferencePlane::GALACTIC => "GALACTIC",
            ReferencePlane::FK4 => "FK4",
        }
    }

}



impl std::fmt::Display for ReferencePlane {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            ReferencePlane::J2000 => write!(f, "J2000"),
            ReferencePlane::ECLIPJ2000 => write!(f, "ECLIPJ2000"),
            ReferencePlane::INVARIABLE => write!(f, "INVARIABLE"),
            ReferencePlane::GALACTIC => write!(f, "GALACTIC"),
            ReferencePlane::FK4 => write!(f, "FK4"),
        }
    }
}