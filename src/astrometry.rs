use crate::time::Time;
use crate::statevector::StateVector;
use crate::spacerock::SpaceRock;

pub struct Astrometry {
    pub ra: f64,
    pub dec: f64,
    pub epoch: Time,
    pub ra_rate: f64,
    pub dec_rate: f64,
    pub topocentric_state: StateVector,
    pub observer: SpaceRock,
    pub name: String,
    pub H: Option<f64>,
    pub Gslope: Option<f64>,
}

impl Astrometry {

    pub fn mag(&self, sun: &SpaceRock) -> Option<f64> {
       
        let H = self.H.unwrap();
        

        let Gslope = self.Gslope.unwrap();

        let delta = self.topocentric_state.position.norm();
        let sun_dist = (self.topocentric_state.position + self.observer.position - sun.position).norm();
        let earth_dist = (self.observer.position - sun.position).norm();
    
        let q = (sun_dist.powi(2) + delta.powi(2) - earth_dist) / (2.0 * sun_dist * delta);
    
        let mut beta = 0.0;
        match q {
            q if q <= -1.0 => beta = std::f64::consts::PI,
            q if q >= 1.0 => beta = 0.0,
            _ => beta = q.acos(),
        };
    
        let psi_1 = (-3.332 * ((beta / 2.0).tan()).powf(0.631)).exp();
        let psi_2 = (-1.862 * ((beta / 2.0).tan()).powf(1.218)).exp();
        let mut mag = H + 5.0 * (sun_dist * delta).log10();
    
        if psi_1 == 0.0 && psi_2 == 0.0 {
            return Some(mag);
        }
    
        mag -= 2.5 * ((1.0 - Gslope) * psi_1 + Gslope * psi_2).log10();
    
        return Some(mag);
    }

}
