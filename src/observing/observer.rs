use crate::observing::observatory::Observatory;
use crate::SpaceRock;
use crate::{Time};

use nalgebra::Vector3;

/// An observer at a specific location and time, which may be on Earth, 
/// in space, or attached to a moving body
#[derive(Debug, Clone, PartialEq)]
pub struct Observer {
    pub spacerock: SpaceRock,
    pub observatory: Observatory,
}

impl Observer {

    pub fn position(&self) -> Vector3<f64> {
        self.spacerock.position
    }

    pub fn velocity(&self) -> Vector3<f64> {
        self.spacerock.velocity
    }

    pub fn epoch(&self) -> Time {
        self.spacerock.epoch.clone()
    }

    pub fn reference_plane(&self) -> String {
        self.spacerock.reference_plane.to_string()
    }

    pub fn origin(&self) -> String {
        self.spacerock.origin.to_string()
    }

    pub fn lat(&self) -> Option<f64> {
        self.observatory.lat()
    }

    pub fn lon(&self) -> Option<f64> {
        self.observatory.lon()
    }
}
