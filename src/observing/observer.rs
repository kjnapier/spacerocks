use crate::time::Time;
use nalgebra::Vector3;
use crate::constants::ROTATION_MATRICES;

use crate::spacerock::{SpaceRock, CoordinateFrame, Origin};

use std::sync::Arc;

const DEG_TO_RAD: f64 = std::f64::consts::PI / 180.0;

use serde::{Serialize, Deserialize};

#[derive(PartialEq, Clone, Debug, Default, Serialize, Deserialize)]
pub struct Observer {
    pub position: Vector3<f64>,
    pub velocity: Vector3<f64>,
    pub epoch: Time,
    pub frame: CoordinateFrame,
    pub origin: Origin,
    pub lat: Option<f64>,
    pub lon: Option<f64>,
    pub rho: Option<f64>,
}

impl Observer {

    pub fn new(position: Vector3<f64>, velocity: Vector3<f64>, epoch: Time) -> Self {
        Observer {
            position: position,
            velocity: velocity,
            epoch: epoch,
            frame: CoordinateFrame::J2000,
            origin: Origin::SSB,
            lat: None,
            lon: None,
            rho: None,
        }
    }

    pub fn from_ground(position: Vector3<f64>, velocity: Vector3<f64>, epoch: Time, frame: &CoordinateFrame, origin: &Origin, lat: f64, lon: f64, rho: f64) -> Self {

        Observer {
            position: position,
            velocity: velocity,
            epoch: epoch,
            frame: frame.clone(),
            origin: origin.clone(),
            lat: Some(lat),
            lon: Some(lon),
            rho: Some(rho),
        }
    }

    pub fn from_spacerock(rock: &SpaceRock) -> Self {
        Observer {
            position: rock.position.clone(),
            velocity: rock.velocity.clone(),
            epoch: rock.epoch.clone(),
            frame: rock.frame.clone(),
            origin: rock.origin.clone(),
            lat: None,
            lon: None,
            rho: None,
        }
    }

    pub fn local_sidereal_time(&self) -> Result<f64, Box<dyn std::error::Error>> {

        let lon = match self.lon {
            Some(lon) => lon,
            None => return Err("Longitude not set. Cannot compute local sidereal time".into()),
        };

        let t = (self.epoch.epoch - 2451545.0) / 36525.0;
        let mut theta = 280.46061837 + 360.98564736629 * (self.epoch.epoch - 2451545.0) + (0.000387933 * t * t) - (t * t * t / 38710000.0);
        theta *= DEG_TO_RAD;
        return Ok(theta + lon)
    }

    pub fn change_frame(&mut self, frame: &CoordinateFrame) -> Result<(), Box<dyn std::error::Error>> {

        if frame == &self.frame {
            return Ok(());
        }

        let inv = frame.get_rotation_matrix().try_inverse().ok_or("Could not invert rotation matrix")?;
        let rot = &self.frame.get_rotation_matrix() * inv;

        self.position = rot * self.position;
        self.velocity = rot * self.velocity;
        self.frame = frame.clone();
        
        Ok(())
    }

}