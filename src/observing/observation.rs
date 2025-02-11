use crate::{Time, Observer};

use nalgebra::{Vector3, DMatrix};

/// Types of astronomical observations, each containing different measured quantities.
#[derive(Debug, Clone, PartialEq)]
pub enum ObservationType {
    Astrometry { ra: f64, dec: f64 },
    Streak { ra: f64, dec: f64, ra_rate: f64, dec_rate: f64 },
    Radar { ra: f64, dec: f64, range: f64, range_rate: f64 },
    Complete { ra: f64, dec: f64, ra_rate: f64, dec_rate: f64, range: f64, range_rate: f64 },
}

/// An astronomical observation of an object from a specific observer.
#[derive(Debug, Clone, PartialEq)]
pub struct Observation {
    pub epoch: Time,
    pub observation_type: ObservationType,
    pub observer: Observer,
    pub inverse_covariance: Option<DMatrix<f64>>,
    pub mag: Option<f64>,
    pub mag_err: Option<f64>,
    // pub filter: Option<String>,
    // pub obsid: Option<String>,
}

impl Observation {
    /// Creates a new Observation 
    ///
    /// # Arguments
    /// * `epoch` - Time of observation
    /// * `observation_type` - Type of observation (Astrometry, Streak, Radar, or Complete)
    /// * `observer` - Observer making the measurement
    /// * `inverse_covariance` - Optional inverse of measurement covariance matrix
    /// * `mag` - Optional visual magnitude
    /// * `mag_err` - Optional magnitude uncertainty    
    pub fn new(epoch: Time, observation_type: ObservationType, observer: Observer, inverse_covariance: Option<DMatrix<f64>>, mag: Option<f64>, mag_err: Option<f64>) -> Observation {
        Observation { epoch, observation_type, observer, inverse_covariance, mag, mag_err }
    }

    /// Creates an astrometric observation from position measurements
    ///
    /// # Arguments
    /// * `epoch` - Time of observation
    /// * `ra` - Right ascension in radians
    /// * `dec` - Declination in radians
    /// * `observer` - Observer making the measurement
    /// * `covariance` - Optional 2x2 covariance matrix for [ra, dec]
    /// * `mag` - Optional visual magnitude
    /// * `mag_err` - Optional magnitude uncertainty
    ///
    /// # Returns
    /// * `Result<Observation, Box<dyn std::error::Error>>` - The created observation or an error
    pub fn from_astrometry(epoch: Time, ra: f64, dec: f64, observer: Observer, covariance: Option<[[f64; 2]; 2]>, mag: Option<f64>, mag_err: Option<f64>) -> Result<Observation, Box<dyn std::error::Error>> {

        let cov = covariance.map(|cov| {
            DMatrix::from_fn(2, 2, |r, c| cov[r][c])
        });
        let inv_cov = cov.clone().map(|cov| {
            cov.clone().try_inverse().ok_or("Covariance matrix is not invertible").unwrap()
        });
        
        Ok(Observation::new(epoch, ObservationType::Astrometry { ra, dec }, observer, inv_cov, mag, mag_err))
    }

    /// Creates a streak observation including apparent motion
    ///
    /// # Arguments
    /// * `epoch` - Time of observation
    /// * `ra` - Right ascension in radians
    /// * `dec` - Declination in radians
    /// * `ra_rate` - Right ascension rate in rad/day
    /// * `dec_rate` - Declination rate in rad/day
    /// * `observer` - Observer making the measurement
    /// * `covariance` - Optional 4x4 covariance matrix for [ra, dec, ra_rate, dec_rate]
    /// * `mag` - Optional visual magnitude
    /// * `mag_err` - Optional magnitude uncertainty
    ///
    /// # Returns
    /// * `Result<Observation, Box<dyn std::error::Error>>` - The created observation or an error
    pub fn from_streak(epoch: Time, ra: f64, dec: f64, ra_rate: f64, dec_rate: f64, observer: Observer, covariance: Option<[[f64; 4]; 4]>, mag: Option<f64>, mag_err: Option<f64>) -> Result<Observation, Box<dyn std::error::Error>> {
        let cov = covariance.map(|cov| {
            DMatrix::from_fn(4, 4, |r, c| cov[r][c])
        });

        let inv_cov = cov.clone().map(|cov| {
            cov.clone().try_inverse().ok_or("Covariance matrix is not invertible").unwrap()
        });
        Ok(Observation::new(epoch, ObservationType::Streak { ra, dec, ra_rate, dec_rate }, observer, inv_cov, mag, mag_err))
    }

    /// Creates a complete observation with position, motion, and range information
    ///
    /// # Arguments
    /// * `epoch` - Time of observation
    /// * `ra` - Right ascension in radians
    /// * `dec` - Declination in radians
    /// * `ra_rate` - Right ascension rate in rad/day
    /// * `dec_rate` - Declination rate in rad/day
    /// * `range` - Distance to target in AU
    /// * `range_rate` - Range rate in AU/day
    /// * `observer` - Observer making the measurement
    /// * `covariance` - Optional 6x6 covariance matrix for [ra, dec, ra_rate, dec_rate, range, range_rate]
    /// * `mag` - Optional visual magnitude
    /// * `mag_err` - Optional magnitude uncertainty
    ///
    /// # Returns
    /// * `Result<Observation, Box<dyn std::error::Error>>` - The created observation or an error
    pub fn from_complete(epoch: Time, ra: f64, dec: f64, ra_rate: f64, dec_rate: f64, range: f64, range_rate: f64, observer: Observer, covariance: Option<[[f64; 6]; 6]>, mag: Option<f64>, mag_err: Option<f64>) -> Result<Observation, Box<dyn std::error::Error>> {
        let cov = covariance.map(|cov| {
            DMatrix::from_fn(6, 6, |r, c| cov[r][c])
        });

        if cov.is_some() && cov.as_ref().unwrap().determinant() == 0.0 {
            return Err("Covariance matrix is singular".into());
        }

        let inv_cov = cov.clone().map(|cov| {
            cov.clone().try_inverse().ok_or("Covariance matrix is not invertible").unwrap()
        });
        Ok(Observation::new(epoch, ObservationType::Complete { ra, dec, ra_rate, dec_rate, range, range_rate }, observer, inv_cov, mag, mag_err))
    }

    pub fn ra(&self) -> f64 {
        match self.observation_type {
            ObservationType::Astrometry { ra, .. } => ra,
            ObservationType::Streak { ra, .. } => ra,
            ObservationType::Radar { ra, .. } => ra,
            ObservationType::Complete { ra, .. } => ra,
        }
    }

    pub fn dec(&self) -> f64 {
        match self.observation_type {
            ObservationType::Astrometry { dec, .. } => dec,
            ObservationType::Streak { dec, .. } => dec,
            ObservationType::Radar { dec, .. } => dec,
            ObservationType::Complete { dec, .. } => dec,
        }
    }

    pub fn ra_rate(&self) -> Option<f64> {
        match self.observation_type {
            ObservationType::Astrometry { .. } => None,
            ObservationType::Streak { ra_rate, .. } => Some(ra_rate),
            ObservationType::Radar { .. } => None,
            ObservationType::Complete { ra_rate, .. } => Some(ra_rate),
        }
    }

    pub fn dec_rate(&self) -> Option<f64> {
        match self.observation_type {
            ObservationType::Astrometry { .. } => None,
            ObservationType::Streak { dec_rate, .. } => Some(dec_rate),
            ObservationType::Radar { .. } => None,
            ObservationType::Complete { dec_rate, .. } => Some(dec_rate),
        }
    }

    pub fn range(&self) -> Option<f64> {
        match self.observation_type {
            ObservationType::Astrometry { .. } => None,
            ObservationType::Streak { .. } => None,
            ObservationType::Radar { range, .. } => Some(range),
            ObservationType::Complete { range, .. } => Some(range),
        }
    }

    pub fn range_rate(&self) -> Option<f64> {
        match self.observation_type {
            ObservationType::Astrometry { .. } => None,
            ObservationType::Streak { .. } => None,
            ObservationType::Radar { range_rate, .. } => Some(range_rate),
            ObservationType::Complete { range_rate, .. } => Some(range_rate),
        }
    }

    pub fn mag(&self) -> Option<f64> {
        self.mag
    }

    pub fn mag_err(&self) -> Option<f64> {
        self.mag_err
    }

    pub fn set_mag(&mut self, mag: f64) {
        self.mag = Some(mag);
    }

    pub fn set_mag_err(&mut self, mag_err: f64) {
        self.mag_err = Some(mag_err);
    }

    pub fn set_covariance(&mut self, covariance: DMatrix<f64>) {
        self.inverse_covariance = Some(covariance.try_inverse().unwrap());
    }

    pub fn inverse_covariance(&self) -> Option<&DMatrix<f64>> {
        self.inverse_covariance.as_ref()
    }

    pub fn covariance(&self) -> Option<DMatrix<f64>> {
        self.inverse_covariance.as_ref().map(|cov| cov.clone().try_inverse().unwrap())
    }

    pub fn proper_motion(&self) -> Option<f64> {
        let ra_rate = self.ra_rate()?;
        let dec_rate = self.dec_rate()?;
        Some((ra_rate.powi(2) * (self.dec().cos()).powi(2) + dec_rate.powi(2)).sqrt())
    }

    /// Returns unit vector in direction of observation
    pub fn pointing(&self) -> Vector3<f64> {
        let ra = self.ra();
        let dec = self.dec();
        Vector3::new(dec.cos() * ra.cos(), dec.cos() * ra.sin(), dec.sin())
    }
}

// implement a display trait for Observation
impl std::fmt::Display for Observation {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self.observation_type {
            ObservationType::Astrometry { ra, dec, .. } => write!(f, "Astrometric observation at epoch {} with RA: {} and Dec: {}", self.epoch, ra, dec),
            ObservationType::Streak { ra, dec, ra_rate, dec_rate, .. } => write!(f, "Streak observation at epoch {} with RA: {}, Dec: {}, RA rate: {}, Dec rate: {}", self.epoch, ra, dec, ra_rate, dec_rate),
            ObservationType::Radar { ra, dec, range, range_rate, .. } => write!(f, "Radar observation at epoch {} with RA: {}, Dec: {}, Range: {}, Range rate: {}", self.epoch, ra, dec, range, range_rate),
            ObservationType::Complete { ra, dec, ra_rate, dec_rate, range, range_rate, .. } => write!(f, "Complete observation at epoch {} with RA: {}, Dec: {}, RA rate: {}, Dec rate: {}, Range: {}, Range rate: {}", self.epoch, ra, dec, ra_rate, dec_rate, range, range_rate),
        }
    }
}
