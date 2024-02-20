use crate::constants::*;
use crate::StateVector;
use crate::KeplerOrbit;
use crate::Properties;
use crate::observing::Observer;

use crate::observing::{Detection, Observatory};

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
    pub mass: f64,

    pub position: Vector3<f64>,
    pub velocity: Vector3<f64>,
    pub acceleration: Vector3<f64>,

    pub mu: Option<f64>,
    pub orbit: Option<KeplerOrbit>,
    pub properties: Option<Properties>,

}

#[allow(dead_code)]
impl SpaceRock {

    // Instantiation Methods
    pub fn from_spice(name: &str, epoch: &Time, frame: &str, origin: &str) -> Self {
        let mut ep = epoch.clone();
        ep.to_utc();
        let et = spice::str2et(&format!("JD{epoch} UTC", epoch=ep.epoch));
        let (state, _) = spice::spkezr(name, et, frame, "NONE", origin);
        let position = Vector3::new(state[0], state[1], state[2]) * KM_TO_AU;
        let velocity = Vector3::new(state[3], state[4], state[5]) * KM_TO_AU * SECONDS_PER_DAY;
        let acceleration = Vector3::new(0.0, 0.0, 0.0);
        let mu = Some(MU_BARY);

        let mass = match MASSES.get(name) {
            Some(m) => *m,
            None => 0.0,
        };

        if mass == 0.0 {
            println!("Could not find mass for {}. Setting mass to 0.", name);
        }

        SpaceRock {
            name: name.to_string(), 
            position: position,
            velocity: velocity,
            acceleration: acceleration,
            mu: mu,
            epoch: epoch.clone(),
            frame: frame.to_string(),
            origin: origin.to_string(),
            orbit: Some(KeplerOrbit::from_xyz(StateVector {position: position, velocity: velocity})),
            mass: mass,
            properties: None,
        }
    }

    pub fn from_xyz(name: &str, x: f64, y: f64, z: f64, vx: f64, vy: f64, vz: f64, epoch: Time, frame: &str, origin: &str) -> Self {
        let position = Vector3::new(x, y, z);
        let velocity = Vector3::new(vx, vy, vz);
        let acceleration = Vector3::new(0.0, 0.0, 0.0);
        SpaceRock {
            name: name.to_string(),
            position: position,
            velocity: velocity,
            acceleration: acceleration,
            epoch: epoch,
            mu: Some(MU_BARY),
            frame: frame.to_string(),
            origin: origin.to_string(),
            orbit: Some(KeplerOrbit::from_xyz(StateVector {position: position, velocity: velocity})),
            mass: 0.0,
            properties: None,
        }
    }

    pub fn from_state(name: &str, state: StateVector, epoch: Time, frame: &str, origin: &str) -> Self {
        let position = state.position;
        let velocity = state.velocity;
        let acceleration = Vector3::new(0.0, 0.0, 0.0);
        SpaceRock {
            name: name.to_string(),
            position: position,
            velocity: velocity,
            acceleration: acceleration,
            epoch: epoch,
            mu: Some(MU_BARY),
            frame: frame.to_string(),
            origin: origin.to_string(),
            orbit: Some(KeplerOrbit::from_xyz(StateVector {position: position, velocity: velocity})),
            mass: 0.0,
            properties: None,
        }
    }

    pub fn from_kepler(name: &str, orbit: KeplerOrbit, epoch: Time, frame: &str, origin: &str) -> Self {
        let state = calc_xyz_from_kepM(orbit.a, orbit.e, orbit.inc, orbit.arg, orbit.node, orbit.M());
        let acceleration = Vector3::new(0.0, 0.0, 0.0);
        SpaceRock {
            name: name.to_string(),
            position: state.position,
            velocity: state.velocity,
            acceleration: acceleration,
            epoch: epoch,
            mu: Some(MU_BARY),
            frame: frame.to_string(),
            origin: origin.to_string(),
            orbit: Some(orbit),
            mass: 0.0,
            properties: None,
        }
    }

    pub fn from_horizons(name: &str, epoch: &Time, frame: &str, origin: &str) -> Result<Self, Box<dyn std::error::Error>> {

        let client = reqwest::blocking::Client::new();

        let mut params = HashMap::new();
        // make the command the name, but with a single quote around it as a string

        let command_str = format!("'{}'", name);

        let mut ep = epoch.clone();
        ep.to_tdb();

        params.insert("command", command_str.as_str());
        let start_time = format!("'JD {}TDB'", ep.epoch);
        let stop_time = format!("'JD {}'", ep.epoch + 1.0);

        match frame {
            "J2000" => {
                params.insert("ref_system", "'J2000'");
                params.insert("ref_plane", "'frame'");
            },
            "ECLIPJ2000" => {
                params.insert("ref_system", "'J2000'");
                params.insert("ref_plane", "'ecliptic'");
            },
            _ => {
                return Err("Frame not recognized".into());
            }
        }

        let center = format!("'@{}'", origin);

        params.insert("start_time", start_time.as_str());
        params.insert("stop_time", stop_time.as_str());
        params.insert("step_size", "'1d'");

        params.insert("make_ephem", "'yes'");
        params.insert("ephem_type", "'vectors'");
        params.insert("center", center.as_str());

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
        let acceleration = Vector3::new(0.0, 0.0, 0.0);

        let rock = SpaceRock {
            name: name.to_string(),
            position: position,
            velocity: velocity,
            acceleration: acceleration,
            mu: Some(MU_BARY),
            epoch: epoch,
            frame: frame.to_string(),
            origin: origin.to_string(),
            orbit: Some(KeplerOrbit::from_xyz(StateVector {position: position, velocity: velocity})),
            mass: 0.0,
            properties: None,
        };
        return Ok(rock);
    }

    pub fn random() -> Self {
        let mut rng = rand::thread_rng();
        let a = rng.gen_range(40.0..50.0);
        let e = rng.gen_range(0.0..0.3);
        let inc = rng.gen_range(0.0..std::f64::consts::PI/3.0);
        let arg = rng.gen_range(0.0..2.0 * std::f64::consts::PI);
        let node = rng.gen_range(0.0..2.0 * std::f64::consts::PI);
        let f = rng.gen_range(0.0..2.0 * std::f64::consts::PI);

        // uuid for name
        let name = format!("{}", uuid::Uuid::new_v4().simple());

        SpaceRock::from_kepler(&name, KeplerOrbit::new(a, e, inc, arg, node, f), Time::now(), "ECLIPJ2000", "SSB")
    }

    // Methods

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
            // self.orbit = None;
            self.calculate_orbit();
        }        
    }


    pub fn observe(&mut self, observer: &Observer) -> Detection {

        // Make sure the rock is in the J2000 frame
        let original_frame = self.frame.clone();
        self.change_frame("J2000");

        // Calculate the topocentric state, correct for light travel time
        let cr = correct_for_ltt(&self, observer);

        // Calaculate the ra, and dec
        let mut ra = cr.position.y.atan2(cr.position.x);
        if ra < 0.0 {
            ra += 2.0 * std::f64::consts::PI;
        }
        let dec = (cr.position.z / cr.position.norm()).asin();

        // Calculate the ra and dec rates
        let xi = cr.position.x.powi(2) + cr.position.y.powi(2);
        let ra_rate = - (cr.position.y * cr.velocity.x - cr.position.x * cr.velocity.y) / xi;
        let num = -cr.position.z * (cr.position.x * cr.velocity.x + cr.position.y * cr.velocity.y) + xi * cr.velocity.z;
        let denom = xi.sqrt() * cr.position.norm_squared();
        let dec_rate = num / denom;

        // calculate the topocentric range and range rate
        let rho = cr.position.norm();
        let rho_rate = cr.position.dot(&cr.velocity) / rho;

        // construct the detection
        let obs = Detection {
            ra: ra,
            dec: dec,
            ra_rate: Some(ra_rate),
            dec_rate: Some(dec_rate),
            rho: Some(rho),
            rho_rate: Some(rho_rate),
            epoch: self.epoch.clone(),
            observer: observer.clone(),
            name: self.name.clone(),

            mag: None,
            filter: None,
            ra_uncertainty: None,
            dec_uncertainty: None,
            ra_rate_uncertainty: None,
            dec_rate_uncertainty: None,
            rho_uncertainty: None,
            rho_rate_uncertainty: None,
            mag_uncertainty: None,
        };

        // Change the frame back to the original frame
        self.change_frame(&original_frame);
        
        return obs;
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

        // Change the origin to a new body, and set the new mu.

        let origin_position = origin.position;
        let origin_velocity = origin.velocity;

        self.position -= origin_position;
        self.velocity -= origin_velocity;

        self.mu = Some(origin.mass * gravitational_constant);
        self.origin = origin.name.to_string();
        self.orbit = Some(KeplerOrbit::from_xyz(StateVector {position: self.position, velocity: self.velocity}));
    }

    pub fn to_ssb(&mut self) {
        // get the ssb from spice
        let mut ssb = SpaceRock::from_spice("ssb", &self.epoch, self.frame.as_str(), self.origin.as_str());
        ssb.mass = MU_BARY / gravitational_constant;
        self.change_origin(&ssb);
    }

    pub fn to_helio(&mut self) {
        // get the sun from spice
        let sun = SpaceRock::from_spice("sun", &self.epoch, self.frame.as_str(), self.origin.as_str());
        self.change_origin(&sun);
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