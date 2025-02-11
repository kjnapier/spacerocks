//! Physical and observational properties of a celestial body.
use serde::{Serialize, Deserialize};

/// Represents the physical and observational properties of a celestial body.
///
/// These properties are optional characteristics that can be associated with a SpaceRock:
/// * Mass is in Solar Mass (M☉)
/// * Absolute magnitude (H) follows the IAU photometric system
/// * G-slope is the magnitude-phase relation slope parameter (defaults to 0.15)
/// * Radius is in kilometers (km)
/// * Albedo is the geometric albedo (dimensionless, between 0 and 1)
#[derive(Debug, Clone, PartialEq, Default, Serialize, Deserialize)]
pub struct Properties {
    pub mass: Option<f64>,
    pub absolute_magnitude: Option<f64>,
    pub gslope: Option<f64>,
    pub radius: Option<f64>,
    pub albedo: Option<f64>,
}