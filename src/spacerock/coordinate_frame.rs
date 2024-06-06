use nalgebra::Matrix3;
use crate::constants::{ROTATION_J2000, ROTATION_ECLIPJ2000, ROTATION_INVARIABLE, ROTATION_GALACTIC, ROTATION_FK4};

use serde::{Serialize, Deserialize};

#[derive(PartialEq, Debug, Clone, Serialize, Deserialize)]
pub enum CoordinateFrame {
    J2000,
    ECLIPJ2000,
    INVARIABLE,
    GALACTIC,
    FK4,
}

impl CoordinateFrame {
    pub fn from_str(s: &str) -> Result<Self, String> {
        match s {
            "J2000" => Ok(CoordinateFrame::J2000),
            "ECLIPJ2000" => Ok(CoordinateFrame::ECLIPJ2000),
            "INVARIABLE" => Ok(CoordinateFrame::INVARIABLE),
            "GALACTIC" => Ok(CoordinateFrame::GALACTIC),
            "FK4" => Ok(CoordinateFrame::FK4),
            _ => Err(format!("Invalid frame: {}", s))
        }
    }

    // get the rotation matrix for the frame
    pub fn get_rotation_matrix(&self) -> Matrix3<f64> {
        match self {
            CoordinateFrame::J2000 => ROTATION_J2000,
            CoordinateFrame::ECLIPJ2000 => ROTATION_ECLIPJ2000,
            CoordinateFrame::INVARIABLE => ROTATION_INVARIABLE,
            CoordinateFrame::GALACTIC => ROTATION_GALACTIC,
            CoordinateFrame::FK4 => ROTATION_FK4,
        }
    }

    pub fn as_str(&self) -> &str {
        match self {
            CoordinateFrame::J2000 => "J2000",
            CoordinateFrame::ECLIPJ2000 => "ECLIPJ2000",
            CoordinateFrame::INVARIABLE => "INVARIABLE",
            CoordinateFrame::GALACTIC => "GALACTIC",
            CoordinateFrame::FK4 => "FK4",
        }
    }

}

impl Default for CoordinateFrame {
    fn default() -> Self {
        CoordinateFrame::ECLIPJ2000
    }
}


impl std::fmt::Display for CoordinateFrame {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            CoordinateFrame::J2000 => write!(f, "J2000"),
            CoordinateFrame::ECLIPJ2000 => write!(f, "ECLIPJ2000"),
            CoordinateFrame::INVARIABLE => write!(f, "INVARIABLE"),
            CoordinateFrame::GALACTIC => write!(f, "GALACTIC"),
            CoordinateFrame::FK4 => write!(f, "FK4"),
        }
    }
}