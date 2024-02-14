use crate::time::Time;
use nalgebra::Vector3;
use crate::constants::ROTATION_MATRICES;

const DEG_TO_RAD: f64 = std::f64::consts::PI / 180.0;

#[derive(PartialEq, Clone, Debug)]
pub struct Observer {
    pub position: Vector3<f64>,
    pub velocity: Vector3<f64>,
    pub epoch: Time,
    pub frame: String,
    pub lat: Option<f64>,
    pub lon: Option<f64>,
    pub elevation: Option<f64>,
}

impl Observer {

    pub fn new(position: Vector3<f64>, velocity: Vector3<f64>, epoch: Time) -> Self {
        Observer {
            position: position,
            velocity: velocity,
            epoch: epoch,
            frame: "J2000".to_string(),
            lat: None,
            lon: None,
            elevation: None,
        }
    }

    pub fn from_ground(position: Vector3<f64>, velocity: Vector3<f64>, epoch: Time, frame: &str, lat: f64, lon: f64, elevation: f64) -> Self {
        Observer {
            position: position,
            velocity: velocity,
            epoch: epoch,
            frame: frame.to_string(),
            lat: Some(lat),
            lon: Some(lon),
            elevation: Some(elevation),
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

    pub fn change_frame(&mut self, frame: &str) {
        if frame != self.frame {
            let inv = ROTATION_MATRICES[&self.frame].try_inverse().unwrap();
            let rot = ROTATION_MATRICES[frame] * inv;
            self.position = rot * self.position;
            self.velocity = rot * self.velocity;
            self.frame = frame.to_string();
        }
    }

}