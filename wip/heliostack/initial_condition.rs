use crate::transforms::{calc_conic_anomaly_from_mean_anomaly, calc_conic_anomaly_from_true_anomaly, calc_mean_anomaly_from_conic_anomaly, calc_true_anomaly_from_mean_anomaly};
use crate::Time;
use crate::OrbitType;

pub fn nice_acos(x: f64) -> f64 {
    if x > 1.0 {
        0.0
    }
    else if x < -1.0 {
        std::f64::consts::PI
    }
    else {
        x.acos()
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct InitialCondition {
    pub id: String,
    pub q: f64,
    pub e: f64,
    pub inc: f64,
    pub mean_anomaly_0: f64,
    pub true_anomaly_0: f64,
    pub r: f64,
    pub vr: f64,
    pub vo: f64,
    pub epoch: Time,
    pub kappa: f64,
    pub mu: f64,
    pub orbit_type: OrbitType
}

impl InitialCondition {

    pub fn from_kep(id: String, q: f64, e: f64, inc: f64, true_anomaly: f64, epoch: Time, kappa:f64, mu: f64) -> Result<Self, Box<dyn std::error::Error>> {

        let orbit_type = OrbitType::from_eccentricity(e, 1e-10)?;

        // validate here that the true anomaly is consistent with the eccentricity.
        if e > 1.0 {
            let true_anomaly_infinity = nice_acos(-1.0 / e);
            if true_anomaly.abs() > true_anomaly_infinity {
                return Err("True anomaly is out of bounds for the given eccentricity.".to_string().into());
            }
        }

        let p = q * (1.0 + e); // semi-latus rectum
        let h = (p * mu).sqrt(); 
        let r = p / (1.0 + e * true_anomaly.cos());
        let vo = h / r;
        let vr = mu * true_anomaly.sin() / h;

        let eccentric_anomaly = calc_conic_anomaly_from_true_anomaly(e, true_anomaly)?;
        let mean_anomaly = calc_mean_anomaly_from_conic_anomaly(e, eccentric_anomaly)?;
   
        Ok(Self {
            id,
            q,
            e,
            inc,
            mean_anomaly_0: mean_anomaly,
            true_anomaly_0: true_anomaly,
            r,
            vr,
            vo,
            epoch,
            kappa,
            mu,
            orbit_type
        })
    }

    pub fn r(&self, epoch: &Time) -> Result<f64, Box<dyn std::error::Error>> {
        let M = self.M(epoch);
        let true_anomaly = calc_true_anomaly_from_mean_anomaly(self.e, M)?;
        let p = self.q * (1.0 + self.e); // semi-latus rectum
        let r = p / (1.0 + self.e * true_anomaly.cos());
        Ok(r)
    }

    pub fn M(&self, epoch: &Time) -> f64 {
        let dt = epoch.epoch - self.epoch.epoch;
        let M = self.mean_anomaly_0 + self.n() * dt;
        M
    }

    pub fn n(&self) -> f64 {
        match self.a() {
            Some(a) => (self.mu / a.powi(3)).sqrt(),
            None => (self.mu / self.q.powi(3)).sqrt()
        }
    }

    pub fn a(&self) -> Option<f64> {
        match self.orbit_type {
            OrbitType::Hyperbolic => Some(self.q / (self.e - 1.0)),
            OrbitType::Parabolic => None,
            OrbitType::Elliptical => Some(self.q / (1.0 - self.e)),
            OrbitType::Circular => Some(self.q),
            _ => None
        }
    }

}
