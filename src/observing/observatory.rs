use crate::SpaceRock;
use crate::CoordinateFrame;
use crate::spacerock::Origin;

use crate::observing::observer::Observer;
use crate::constants::{DEG_TO_RAD, M_TO_AU, EQUAT_RAD};
use crate::constants::OBSERVATORIES;
use crate::time::Time;

use nalgebra::Vector3;

pub struct Observatory {
    pub lon: f64,
    pub rho_cos_lat: f64,
    pub rho_sin_lat: f64,
}

impl Observatory {

    pub fn from_obscode(obscode: &str) -> Result<Self, &'static str> {
        let obscode = obscode.to_uppercase();

        if let Some(obs) = OBSERVATORIES.get(&obscode) {
            return Ok(Observatory::from_parallax(obs.0, obs.1, obs.2));
        }

        return Err("Observatory not found")
    }

    pub fn from_parallax(lon: f64, rho_cos_lat: f64, rho_sin_lat: f64) -> Self {
        return Observatory {
            lon: lon,
            rho_cos_lat: rho_cos_lat,
            rho_sin_lat: rho_sin_lat
        }
    }


    pub fn lat(&self) -> f64 {
        return self.rho_sin_lat.atan2(self.rho_cos_lat)
    }

    pub fn rho(&self) -> f64 {
        return (self.rho_sin_lat * self.rho_sin_lat + self.rho_cos_lat * self.rho_cos_lat).sqrt()
    }

    pub fn at(&self, epoch: &Time, origin: &Origin) -> Observer {
        let mut earth = SpaceRock::from_spice("earth", epoch, &CoordinateFrame::J2000, origin);
        let [d_pos, d_vel] = compute_topocentric_correction(self.lon, self.rho_sin_lat, self.rho_cos_lat, epoch.epoch);
        earth.position += d_pos;
        earth.velocity += d_vel;
        Observer::from_ground(earth.position, earth.velocity, earth.epoch, &earth.frame, &earth.origin, self.lat(), self.lon, self.rho())
    }

    pub fn local_sidereal_time(&self, epoch: f64) -> f64 {
        return compute_local_sidereal_time(epoch, self.lon)
    }

}

fn compute_local_sidereal_time(epoch: f64, lon: f64) -> f64 {
    let t = (epoch - 2451545.0) / 36525.0;
    let mut theta = 280.46061837 + 360.98564736629 * (epoch - 2451545.0) + (0.000387933 * t * t) - (t * t * t / 38710000.0);
    theta *= DEG_TO_RAD;
    return theta + lon
}

fn sidereal_rate(epoch: f64) -> f64 {
    let T = (epoch - 2451545.0) / 36525.0;
    let Tprime = 1.0 / 36525.0;
    let theta_dot = 360.98564736629 + 2.0 * 0.000387933 * T * Tprime + 3.0 * T * T * Tprime / 38710000.0;
    return theta_dot * DEG_TO_RAD
}

fn compute_topocentric_correction(lon: f64, rho_sin_lat: f64, rho_cos_lat: f64, epoch: f64) -> [Vector3<f64>; 2] {

    let phi = compute_local_sidereal_time(epoch, lon);
    let phi_rate = sidereal_rate(epoch);
    
    let sin_lon = phi.sin();
    let cos_lon = phi.cos();

    let sin_lon_prime = cos_lon * phi_rate;
    let cos_lon_prime = -sin_lon * phi_rate;
        
    let dx = rho_cos_lat * cos_lon * EQUAT_RAD;
    let dy = rho_cos_lat * sin_lon * EQUAT_RAD;
    let dz = rho_sin_lat * EQUAT_RAD;

    let dvx = cos_lon_prime * rho_cos_lat * EQUAT_RAD;
    let dvy = sin_lon_prime * rho_cos_lat * EQUAT_RAD;
    let dvz = 0.0;

    let d_pos = Vector3::new(dx, dy, dz) * M_TO_AU; // AU
    let d_vel = Vector3::new(dvx, dvy, dvz) * M_TO_AU; // AU/day

    return [d_pos, d_vel];

}

// fn compute_topocentric_correction_lonlatel(lon: f64, lat: f64, elevation: f64, epoch: f64) -> [Vector3<f64>; 2] {
        
//     let observer_lat: f64 = lat * DEG_TO_RAD;
//     let observer_lon = lon * DEG_TO_RAD;
//     let sin_lat = observer_lat.sin();
//     let cos_lat = observer_lat.cos();

//     let phi = compute_local_sidereal_time(epoch, observer_lon);
//     let phi_rate = sidereal_rate(epoch);
    
//     let sin_lon = phi.sin();
//     let cos_lon = phi.cos();

//     let sin_lon_prime = cos_lon * phi_rate;
//     let cos_lon_prime = -sin_lon * phi_rate;
        
//     let mut denom: f64 = O_M_FLATTEN * sin_lat;
//     denom = cos_lat * cos_lat + denom*denom;

//     let mut c_geo: f64 = 1.0 / denom.sqrt();
//     let mut s_geo: f64 = O_M_FLATTEN * O_M_FLATTEN * c_geo;
//     c_geo = c_geo * EQUAT_RAD + elevation;
//     s_geo = s_geo * EQUAT_RAD + elevation;
//     let dx = c_geo * cos_lat * cos_lon;
//     let dy = c_geo * cos_lat * sin_lon;
//     let dz = s_geo * sin_lat;

//     let dvx = cos_lat * (cos_lon_prime * c_geo);
//     let dvy = cos_lat * (sin_lon_prime * c_geo);
//     let dvz = 0.0;

//     let d_pos = Vector3::new(dx, dy, dz) * M_TO_AU; // AU
//     let d_vel = Vector3::new(dvx, dvy, dvz) * M_TO_AU; // AU/day

//     return [d_pos, d_vel];

// }