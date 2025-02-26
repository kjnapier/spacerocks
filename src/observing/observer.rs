use crate::observing::observatory::Observatory;
use crate::SpaceRock;
use crate::{Time};
use crate::coordinates::{ReferencePlane, Origin};
use crate::data::GRAVITATIONAL_CONSTANT;

use nalgebra::Vector3;

/// An observer at a specific location and time, which may be on Earth, 
/// in space, or attached to a moving body.
#[derive(Debug, Clone, PartialEq)]
pub struct Observer {
    pub position: Vector3<f64>,
    pub velocity: Option<Vector3<f64>>,
    pub epoch: Time,
    pub reference_plane: ReferencePlane,
    pub origin: Origin,
    pub observatory: Observatory,
}


impl Observer {

    /// Change the reference plane of the Observer
    ///
    /// # Arguments
    /// * `reference_plane` - The new reference plane
    pub fn change_reference_plane(&mut self, reference_plane: &str) -> Result<(), Box<dyn std::error::Error>> {

        let reference_plane = ReferencePlane::from_str(reference_plane)?;
        if reference_plane == self.reference_plane {
            return Ok(());
        }

        let inv = self.reference_plane.get_rotation_matrix().try_inverse().ok_or("Could not invert rotation matrix")?;
        let rot = reference_plane.get_rotation_matrix() * inv;

        self.position = rot * self.position;

        if let Some(velocity) = self.velocity {
            self.velocity = Some(rot * velocity);
        }
        self.reference_plane = reference_plane;

        Ok(())
    }

    /// Change the origin of the Observer
    ///
    /// # Arguments
    /// * `origin` - The SpaceRock object to change the origin to
    pub fn change_origin(&mut self, origin: &SpaceRock) {

        let origin_position = origin.position;
        self.position -= origin_position;

        if let Some(velocity) = self.velocity {
            self.velocity = Some(velocity - origin.velocity);
        }

        self.origin = Origin::new_custom(origin.mass() * GRAVITATIONAL_CONSTANT, &origin.name);
    }


    /// Returns the latitude of the observatory in radians, if it is a ground-based observatory
    pub fn lat(&self) -> Option<f64> {
        self.observatory.lat()
    }
    
    /// Returns the longitude of the observatory in radians, if it is a ground-based observatory 
    pub fn lon(&self) -> Option<f64> {
        self.observatory.lon()
    }
}
