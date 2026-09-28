use crate::constants::{MU_BARY, SPEED_OF_LIGHT};
use crate::time::Time;
use spaccrateerocks::Observation;

use nalgebra::Vector3;

pub fn cost_function(observation: &Observation, ic: &InitialCondition, rho: f64) -> f64 {
    let ltt = rho / SPEED_OF_LIGHT;
    let ep = Time::new(observation.epoch.epoch - ltt, "utc", "jd");
    let r = ic.r(&ep);
    r.powi(2) - rho.powi(2) - 2.0 * rho * observation.observer_dot_pointing - observation.observer_distance.powi(2)
}

pub fn gradient(observation: &Observation, ic: &InitialCondition, rho: f64) -> f64 {
    // ignoring the term that goes like recession velocity / speed of light
    - 2.0 * (rho + observation.observer_dot_pointing)
}

pub fn optimize_rho(observation: &Observation, ic: &InitialCondition, rho0: f64) -> f64 {
    // Simple implementation of Newton's method to optimize rho.

    let mut rho = rho0;
    let tol = 1e-8;
    let mut cost = cost_function(&observation, &ic, rho);
    let mut grad = 0.0;

    while cost.abs() > tol {
        grad = gradient(&observation, &ic, rho);
        rho -= cost / grad;
        cost = cost_function(&observation, &ic, rho);
    }

    return rho;
}

pub fn ahat(phi: f64, theta: f64) -> Vector3<f64> {
    return Vector3::new(-phi.sin(), phi.cos(), 0.0);
}

pub fn dhat(phi: f64, theta: f64) -> Vector3<f64> {
    return Vector3::new(-theta.sin() * phi.cos(), -theta.sin() * phi.sin(), theta.cos());
}

pub fn solve_unit_vector_given_invariants(observation: &Observation, ic: &InitialCondition) -> Option<[f64; 3]> {

    // First solve for the distance between the observer and the object. This part corrects for the light travel time.
    let p = ic.q * (1.0 + ic.e); // semi-latus rectum
    let r = p / (1.0 + ic.e * ic.true_anomaly.cos());
    let cos_solar_elongation = -observer.position().dot(&observation.pointing());
    let sin_solar_elongation_sq = 1.0 - cos_solar_elongation.powi(2);
    let observer_distance = observer.position().norm();
    let rho0 = observer_distance * cos_solar_elongation + (r.powi(2) - observer_distance.powi(2) * sin_solar_elongation_sq).sqrt()
    let rho = optimize_rho(&observation, &ic, rho0);

    // The new epoch of the observation is the light-corrected value: t - rho / SPEED_OF_LIGHT
    // let epoch = Time::new(observation.epoch.epoch - rho / SPEED_OF_LIGHT, "utc", "jd");
    let epoch = observation.epoch - rho / SPEED_OF_LIGHT;

    // This gives the position of the rock at the epoch that the light was observed, rather than when it was reflected.
    let r_vec = observation.observer.position + rho * observation.pointing;
    let r = r_vec.norm();
    let phi = r_vec.y.atan2(r_vec.x);
    let theta = (r_vec.z / r).asin(); 
    
    
    let vsq = MU_BARY * (2.0 / r - 1.0 / ic.a);
    let vo = ((ic.e.powi(2) - 1.0) / r.powi(2) * (vsq / MU_BARY.powi(2) - 2.0 / (MU_BARY * r)).powi(-1)).sqrt();
    
    // vo = h / r;
    let vr_sq = vsq - vo.powi(2);
    let mut vr = match vr_sq {
        x if x < 0.0 => 0.0,
        _ => vr_sq.sqrt()    
    };
    vr = if ic.M(&epoch) > std::f64::consts::PI { -vr } else { vr };

    // Now calculate the gauss f and g functions
    let dt = ic.epoch.epoch - epoch.epoch;
    let dE = ic.E0 - ic.E(&epoch);
    let f = 1.0 - ic.a / r * (1.0 - dE.cos());
    let g = dt + 1.0 / ic.n() * (dE.sin() - dE);


    // At this point we have: {rhat, ahat, dhat, vr, vo, r}. We should store these intermediate results.
    // All that's left is psi, which depends on the inclination.

    // If the absolute value of the latitude of the object is greater than the inclination, then the orbit is not possible.
    if ic.inc < theta.abs() {
        return None;
    }
 
    let cos_psi = ic.inc.cos() / theta.cos();
    let sin_psi = ic.kappa * (1.0 - cos_psi.powi(2)).sqrt();

    // construct the velocity vector
    let v_vec = vr * r_vec / r + vo * (cos_psi * ahat(phi, theta) + sin_psi * dhat(phi, theta));
    let new_position = f * r_vec + g * v_vec;
    let new_pointing = new_position.normalize();

    Some(new_pointing.into())
}



use crate::transforms::{calc_conic_anomaly_from_mean_anomaly, calc_conic_anomaly_from_true_anomaly, calc_mean_anomaly_from_conic_anomaly, calc_true_anomaly_from_mean_anomaly};
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
