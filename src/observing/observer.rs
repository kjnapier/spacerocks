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
    pub observatory: Option<Observatory>,
    /// Position of the Sun in the same frame and origin as `position` (AU), used for
    /// phase angles and heliocentric distances. `None` means the origin is taken to be the
    /// Sun's position (exact when `origin` is SUN, ~0.01 AU off for SSB).
    pub sun_position: Option<Vector3<f64>>,
}


impl Observer {

    /// Create a new Observer from Cartesian coordinates.
    ///
    /// # Arguments
    /// * `position` - The position of the observer
    /// * `velocity` - The velocity of the observer
    /// * `epoch` - The epoch of the observer
    /// * `reference_plane` - The reference plane of the observer
    /// * `origin` - The origin of the observer
    /// * `observatory` - The observatory of the observer
    pub fn from_xyz(position: Vector3<f64>, velocity: Option<Vector3<f64>>, epoch: Time, reference_plane: ReferencePlane, origin: Origin, observatory: Option<Observatory>) -> Self {
        Observer {
            position,
            velocity,
            epoch,
            reference_plane,
            origin,
            observatory,
            sun_position: None,
        }
    }

    /// Attach the Sun's position (same frame and origin as the observer).
    pub fn with_sun_position(mut self, sun_position: Vector3<f64>) -> Self {
        self.sun_position = Some(sun_position);
        self
    }

    /// Position of the Sun relative to the observer's origin (zero if unknown).
    pub fn sun(&self) -> Vector3<f64> {
        self.sun_position.unwrap_or_else(Vector3::zeros)
    }

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
        if let Some(sun) = self.sun_position {
            self.sun_position = Some(rot * sun);
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
        if let Some(sun) = self.sun_position {
            self.sun_position = Some(sun - origin_position);
        }

        if let Some(velocity) = self.velocity {
            self.velocity = Some(velocity - origin.velocity);
        }

        self.origin = Origin::new_custom(origin.mass() * GRAVITATIONAL_CONSTANT, &origin.name);
    }


    /// Returns the latitude of the observatory in radians, if it is a ground-based observatory
    pub fn lat(&self) -> Option<f64> {
        self.observatory.as_ref()?.lat()
        // self.observatory.lat()
    }
    
    /// Returns the longitude of the observatory in radians, if it is a ground-based observatory 
    pub fn lon(&self) -> Option<f64> {
        self.observatory.as_ref()?.lon()
        // self.observatory.lon()
    }
}
