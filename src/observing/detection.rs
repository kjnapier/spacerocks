use crate::time::Time;
use crate::StateVector;
use crate::SpaceRock;
use crate::observing::observer::Observer;
use nalgebra::Vector3;

#[derive(Debug, Clone, PartialEq)]
pub struct Detection {

    pub name: String,
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

        let lon = match self.observer.lon {
            Some(lon) => lon,
            None => return Err("Observer longitude not set. Cannot compute altaz".into()),
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

    pub fn pointing(&self) -> Vector3<f64> {
        let x = self.ra.cos() * self.dec.cos();
        let y = self.ra.sin() * self.dec.cos();
        let z = self.dec.sin();
        Vector3::new(x, y, z)
    }

    // pub fn mag(&self, sun: &SpaceRock) -> Option<f64> {
       
    //     let H = self.H.unwrap();
        

    //     let Gslope = self.Gslope.unwrap();

    //     let delta = self.topocentric_state.position.norm();
    //     let sun_dist = (self.topocentric_state.position + self.observer.position - sun.position).norm();
    //     let earth_dist = (self.observer.position - sun.position).norm();
    
    //     let q = (sun_dist.powi(2) + delta.powi(2) - earth_dist) / (2.0 * sun_dist * delta);
    
    //     let mut beta = 0.0;
    //     match q {
    //         q if q <= -1.0 => beta = std::f64::consts::PI,
    //         q if q >= 1.0 => beta = 0.0,
    //         _ => beta = q.acos(),
    //     };
    
    //     let psi_1 = (-3.332 * ((beta / 2.0).tan()).powf(0.631)).exp();
    //     let psi_2 = (-1.862 * ((beta / 2.0).tan()).powf(1.218)).exp();
    //     let mut mag = H + 5.0 * (sun_dist * delta).log10();
    
    //     if psi_1 == 0.0 && psi_2 == 0.0 {
    //         return Some(mag);
    //     }
    
    //     mag -= 2.5 * ((1.0 - Gslope) * psi_1 + Gslope * psi_2).log10();
    
    //     return Some(mag);
    // }

}
