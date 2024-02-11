use crate::constants::*;
use crate::statevector::StateVector;
use crate::observatory::Observatory;
use crate::keplerorbit::KeplerOrbit;
use crate::properties::Properties;
use crate::astrometry::Astrometry;

use crate::time::Time;

use crate::transforms::correct_for_ltt;
use crate::transforms::calc_kep_from_xyz;
use crate::transforms::calc_xyz_from_kepM;
use crate::transforms::calc_f_from_E;
use crate::transforms::calc_E_from_M;

use reqwest::blocking::Client;
use serde_json;
use std::collections::HashMap;
use serde::{Serialize, Deserialize};

use nalgebra::Vector3;

use rand;
use rand::Rng;
use rand::distributions::{Uniform};

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SpaceRock {

    pub name: String,
    pub epoch: Time,

    pub frame: String,
    pub origin: String,

    pub position: Vector3<f64>,
    pub velocity: Vector3<f64>,

    pub orbit: Option<KeplerOrbit>,
    pub mass: Option<f64>,
    pub properties: Option<Properties>,

}

#[allow(dead_code)]
impl SpaceRock {

    // Instantiation Methods

    pub fn from_spice(name: &str, epoch: &Time) -> Self {
        let mut ep = epoch.clone();
        ep.to_utc();
        let et = spice::str2et(&format!("JD{epoch} UTC", epoch=ep.epoch));
        let (state, _) = spice::spkezr(name, et, "J2000", "NONE", "SSB");
        let position = Vector3::new(state[0], state[1], state[2]) * KM_TO_AU;
        let velocity = Vector3::new(state[3], state[4], state[5]) * KM_TO_AU * SECONDS_PER_DAY;
        SpaceRock {
            name: name.to_string(), 
            position: position,
            velocity: velocity,
            epoch: epoch.clone(),
            frame: "J2000".to_string(),
            origin: "SSB".to_string(),
            orbit: Some(KeplerOrbit::from_xyz(StateVector {position: position, velocity: velocity})),
            mass: MASSES.get(name).copied(),
            properties: None,
        }
    }

    pub fn from_xyz(name: &str, x: f64, y: f64, z: f64, vx: f64, vy: f64, vz: f64, epoch: Time, frame: &str, origin: &str) -> Self {
        let position = Vector3::new(x, y, z);
        let velocity = Vector3::new(vx, vy, vz);
        SpaceRock {
            name: name.to_string(),
            position: position,
            velocity: velocity,
            epoch: epoch,
            frame: frame.to_string(),
            origin: origin.to_string(),
            orbit: Some(KeplerOrbit::from_xyz(StateVector {position: position, velocity: velocity})),
            mass: None,
            properties: None,
        }
    }

    pub fn from_state(name: &str, state: StateVector, epoch: Time, frame: &str, origin: &str) -> Self {
        let position = state.position;
        let velocity = state.velocity;
        SpaceRock {
            name: name.to_string(),
            position: position,
            velocity: velocity,
            epoch: epoch,
            frame: frame.to_string(),
            origin: origin.to_string(),
            orbit: Some(KeplerOrbit::from_xyz(StateVector {position: position, velocity: velocity})),
            mass: None,
            properties: None,
        }
    }

    pub fn from_kepler(name: &str, orbit: KeplerOrbit, epoch: Time, frame: &str, origin: &str) -> Self {
        let state = calc_xyz_from_kepM(orbit.a, orbit.e, orbit.inc, orbit.arg, orbit.node, orbit.M());
        SpaceRock {
            name: name.to_string(),
            position: state.position,
            velocity: state.velocity,
            epoch: epoch,
            frame: frame.to_string(),
            origin: origin.to_string(),
            orbit: Some(orbit),
            mass: None,
            properties: None,
        }
    }

    // pub fn from_spherical(name: &str, state: SphericalState, epoch: Time, frame: &str, origin: &str) -> Self {
    //     let s = calc_xyz_from_spherical(state);
    //     SpaceRock {
    //         name: name.to_string(),
    //         position: s.position,
    //         velocity: s.velocity,
    //         epoch: epoch,
    //         frame: frame.to_string(),
    //         origin: origin.to_string(),
    //         orbit: Some(KeplerOrbit::from_xyz(StateVector {position: s.position, velocity: s.velocity})),
    //         mass: None,
    //         properties: None,
    //     }
    // }

    pub fn from_horizons(name: &str) -> Result<Self, Box<dyn std::error::Error>> {

        let client = reqwest::blocking::Client::new();

        let mut params = HashMap::new();
        // make the command the name, but with a single quote around it as a string

        let command_str = format!("'{}'", name);

        params.insert("command", command_str.as_str());
        params.insert("start_time", "'2023-03-14'");
        params.insert("stop_time", "'2023-03-15'");
        params.insert("make_ephem", "'yes'");
        params.insert("ephem_type", "'vectors'");
        params.insert("center", "'@ssb'");
        params.insert("ref_plane", "'ecliptic'");
        params.insert("step_size", "'1d'");
        params.insert("ref_system", "'J2000'");
        params.insert("vec_corr", "'None'");
        params.insert("out_units", "'AU-D'");
        params.insert("csv_format", "'yes'");
        params.insert("vec_delta_t", "'no'");
        params.insert("vec_table", "'2x'");
        params.insert("vec_labels", "'no'");

    
        let response = client.get("https://ssd.jpl.nasa.gov/api/horizons.api?")
            .query(&params)
            .send()?;

        let json: serde_json::Value = response.json()?;
        let text = json["result"].as_str();

        let lines: Vec<&str> = text.ok_or_else(| | "No data")?.split('\n').collect();
        let first_data_line = lines.iter().skip_while(|&line| !line.starts_with("$$SOE")).skip(1).next().ok_or("No data")?;
        
        let data: Vec<f64> = first_data_line.split(',').filter_map(|s| s.trim().parse::<f64>().ok()).collect();
        let epoch = Time { epoch: data[0], timescale: "tdb".to_string(), format: "jd".to_string() };
        let (x, y, z, vx, vy, vz) = (data[1], data[2], data[3], data[4], data[5], data[6]);
        // let (dx, dy, dz, dvx, dvy, dvz) = (data[7], data[8], data[9], data[10], data[11], data[12]);

        let position = Vector3::new(x, y, z);
        let velocity = Vector3::new(vx, vy, vz);

        let rock = SpaceRock {
            name: name.to_string(),
            position: position,
            velocity: velocity,
            epoch: epoch,
            frame: "ECLIPJ2000".to_string(),
            origin: "SSB".to_string(),
            orbit: Some(KeplerOrbit::from_xyz(StateVector {position: position, velocity: velocity})),
            mass: None,
            properties: None,
        };
        return Ok(rock);
    }

    pub fn random() -> Self {
        let mut rng = rand::thread_rng();
        let a = rng.gen_range(0.5..100.0);
        let e = rng.gen_range(0.0..1.0);
        let inc = rng.gen_range(0.0..std::f64::consts::PI);
        let arg = rng.gen_range(0.0..2.0 * std::f64::consts::PI);
        let node = rng.gen_range(0.0..2.0 * std::f64::consts::PI);
        let f = rng.gen_range(0.0..2.0 * std::f64::consts::PI);
        // get today's julian date

        SpaceRock::from_kepler("random", KeplerOrbit::new(a, e, inc, arg, node, f), Time::now(), "ECLIPJ2000", "SSB")
    }

    // Operations

    pub fn calculate_kepler(&mut self) {
        self.orbit = Some(KeplerOrbit::from_xyz(StateVector {position: self.position, velocity: self.velocity}));
    }

    pub fn analytic_propagate(&mut self, epoch: &Time) {

        let timescale = &self.epoch.timescale;
        let mut epoch = epoch.clone();

        epoch.change_timescale(timescale);
        let dt = epoch.epoch - self.epoch.epoch;

        // check that self.orbit is not None
        if let Some(orbit) = &self.orbit {
            let dM = orbit.n() * dt;
            let M_new = orbit.M() + dM;
            let new_state = calc_xyz_from_kepM(orbit.a, orbit.e, orbit.inc, orbit.arg, orbit.node, M_new);
            self.position = Vector3::new(new_state.position[0], new_state.position[1], new_state.position[2]);
            self.velocity = Vector3::new(new_state.velocity[0], new_state.velocity[1], new_state.velocity[2]);
            self.epoch = epoch;
            self.calculate_orbit();
        }

        // let dM = self.orbit.n() * dt;

        // let M_new = self.orbit.M() + dM;

        // let new_state = calc_xyz_from_kepM(self.orbit.a, self.orbit.e, self.orbit.inc, self.orbit.arg, self.orbit.node, M_new);

        // self.position = Vector3::new(new_state.position[0], new_state.position[1], new_state.position[2]);
        // self.velocity = Vector3::new(new_state.velocity[0], new_state.velocity[1], new_state.velocity[2]);
        // self.epoch = epoch;
        // self.calculate_orbit();
        
    }


    pub fn observe(&mut self, observer: &SpaceRock) -> Astrometry {

        self.change_frame("J2000");
        let cr = correct_for_ltt(&self, observer);
        let mut ra = cr.position.y.atan2(cr.position.x);
        if ra < 0.0 {
            ra += 2.0 * std::f64::consts::PI;
        }
        let dec = (cr.position.z / cr.position.norm()).asin();

        let xi = cr.position.x.powi(2) + cr.position.y.powi(2);
        let ra_rate = - (cr.position.y * cr.velocity.x - cr.position.x * cr.velocity.y) / xi;
        let num = -cr.position.z * (cr.position.x * cr.velocity.x + cr.position.y * cr.velocity.y) + xi * cr.velocity.z;
        let denom = xi.sqrt() * cr.position.norm_squared();
        let dec_rate = num / denom;

        let mut H = None;
        let mut Gslope = None;
        let mut properties;
        
        if self.properties.is_some() {
            properties = self.properties.clone().unwrap();

            if properties.H.is_some() {
                H = properties.H;
            }
            if properties.Gslope.is_some() {
                Gslope = properties.Gslope;
            }
        }

        let astro = Astrometry {
            ra: ra,
            dec: dec,
            ra_rate: ra_rate,
            dec_rate: dec_rate,
            epoch: self.epoch.clone(),
            topocentric_state: cr,
            observer: observer.clone(),
            name: self.name.clone(),
            H: H,
            Gslope: Gslope,
        };
        
        return astro;
    }

    pub fn change_frame(&mut self, frame: &str) {
        if frame != self.frame {
            let inv = ROTATION_MATRICES[&self.frame].try_inverse().unwrap();
            let rot = ROTATION_MATRICES[frame] * inv;
            self.position = rot * self.position;
            self.velocity = rot * self.velocity;
            self.frame = frame.to_string();
            self.orbit = Some(KeplerOrbit::from_xyz(StateVector {position: self.position, velocity: self.velocity}));
        }
    }

    pub fn change_origin(&mut self, origin: &SpaceRock) {
        if origin.name == self.origin {
            return;
        }
    }

    pub fn calculate_orbit(&mut self) {
        self.orbit = Some(KeplerOrbit::from_xyz(StateVector {position: self.position, velocity: self.velocity}));
    }

    fn r_squared(&self) -> f64 {
        self.position.dot(&self.position)
    }

    pub fn r(&self) -> f64 {
        self.position.norm()
    }

    fn v_squared(&self) -> f64 {
        self.velocity.dot(&self.velocity)
    }

    fn v(&self) -> f64 {
        self.velocity.norm()
    }

    pub fn hvec(&self) -> Vector3<f64> {
        self.position.cross(&self.velocity)
    }

    pub fn h(&self) -> f64 {
        self.hvec().norm()
    }    

}



#[allow(dead_code)]
fn separation(body1: &SpaceRock, body2: &SpaceRock) -> f64 {
    let d_pos = body1.position - body2.position;
    return d_pos.norm();
}
