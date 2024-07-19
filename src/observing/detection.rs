use crate::time::Time;
use crate::observing::observer::Observer;
use nalgebra::Vector3;

use std::sync::Arc;

#[derive(Debug, Clone, PartialEq, Default)]
pub struct Detection {

    pub name: Arc<String>,
    pub epoch: Time,
    pub observer: Observer,

    pub ra: f64,
    pub dec: f64,

    pub ra_rate: Option<f64>,
    pub dec_rate: Option<f64>,

    pub rho: Option<f64>,
    pub rho_rate: Option<f64>,
    
    pub mag: Option<f64>,
    pub filter: Option<String>,

    pub ra_uncertainty: Option<f64>,
    pub dec_uncertainty: Option<f64>,
    pub ra_rate_uncertainty: Option<f64>,
    pub dec_rate_uncertainty: Option<f64>,
    pub rho_uncertainty: Option<f64>,
    pub rho_rate_uncertainty: Option<f64>,
    pub mag_uncertainty: Option<f64>,

}


impl Detection {


    pub fn calc_altaz(&self) -> Result<(f64, f64), Box<dyn std::error::Error>> {

        let lat = match self.observer.lat {
            Some(lat) => lat,
            None => return Err("Observer latitude not set. Cannot compute altaz".into()),
        };

        // let lon = match self.observer.lon {
        //     Some(lon) => lon,
        //     None => return Err("Observer longitude not set. Cannot compute altaz".into()),
        // };

        match self.observer.lon {
            None => return Err("Observer longitude not set. Cannot compute altaz".into()),
            _ => (),
        };

        let local_sidereal_time = self.observer.local_sidereal_time()?;
        let hour_angle = local_sidereal_time - self.ra;
        let sin_alt = lat.sin() * self.dec.sin() + lat.cos() * self.dec.cos() * hour_angle.cos();
        let alt = sin_alt.asin();
        let cos_az = (self.dec.sin() - lat.sin() * alt.sin()) / (lat.cos() * alt.cos());
        let az = cos_az.acos();
        return Ok((alt, (az + std::f64::consts::PI) % (2.0 * std::f64::consts::PI)));
    }

    pub fn proper_motion(&self) -> Option<f64> {
        let ra_rate = match self.ra_rate {
            Some(ra_rate) => ra_rate,
            None => return None,
        };

        let dec_rate = match self.dec_rate {
            Some(dec_rate) => dec_rate,
            None => return None,
        };

        let eta = (ra_rate.powi(2) * (self.dec.cos()).powi(2) + dec_rate.powi(2)).sqrt();
        Some(eta)
    }

    pub fn rho_hat(&self) -> Vector3<f64> {
        let x = self.ra.cos() * self.dec.cos();
        let y = self.ra.sin() * self.dec.cos();
        let z = self.dec.sin();
        Vector3::new(x, y, z)
    }

}
