//! Fundamental data structure for celestial objects.
use crate::{Origin, ReferencePlane, Time, Properties, Observer, Observation};
use crate::observing::{apparent, Apparent};
use crate::constants::*;
use crate::correct_for_ltt;
use crate::SpiceKernel;
use crate::assist::SpiceSimulation;



use crate::state::{self, Elements, State};

use nalgebra::Vector3;

use rand;
use rand::Rng;

use std::collections::HashMap;
use std::time::Duration;

/// A SpaceRock represents a celestial body with a state vector and optional physical properties.
/// 
/// The state vector consists of position (in AU) and velocity (in AU/day) components in a specified
/// reference frame relative to an origin point. Physical properties like mass, absolute magnitude,
/// and albedo are optional.
/// 
/// # State Components
/// * Position and velocity in Cartesian coordinates
/// * Reference plane defining the orientation
/// * Origin point defining the center
/// * Epoch specifying the time of the state vector
/// 
/// # Physical Properties (Optional)
/// * Mass is in Solar Mass (M☉)
/// * Absolute magnitude (H)
/// * Phase slope parameter (G)
/// * Radius (in km)
/// * Albedo
///
/// A SpaceRock can be instantiated from a spice kernel, random keplerian elements, cartesian coordinates, 
/// spherical coordinates, or the JPL Horizons API. It can be propagated in time, and observed from an observer 
/// on Earth. It can also be transformed to the solar system barycenter or the heliocenter.
/// 
/// # Examples
/// ```no_run
/// use spacerocks::SpaceRock;
/// use spacerocks::Time;
///
/// let epoch = Time::now();
/// let asteroid = SpaceRock::from_horizons("Ceres", &epoch, "ECLIPJ2000", "SSB").unwrap();
/// println!("Semi-major axis: {} AU", asteroid.a());
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct SpaceRock {

    pub name: String,
    pub epoch: Time,

    pub reference_plane: ReferencePlane,
    pub origin: Origin,

    pub position: Vector3<f64>,
    pub velocity: Vector3<f64>,

    pub properties: Option<Properties>,
}


impl SpaceRock {

    /// Instantiate a SpaceRock from a spice kernel. A kernel must be loaded before calling this method.
    ///
    /// # Arguments
    /// * `name` - The name of the object
    /// * `epoch` - The epoch of the ephemeris
    /// * `reference_plane` - The coordinate reference_plane 
    /// * `origin` - The origin of the coordinate system
    ///
    /// # Returns
    /// * [`SpaceRock`] 
    ///
    /// # Example
    /// ```no_run
    /// use spacerocks::{SpaceRock, SpiceKernel, Time};
    ///
    /// let kernel = SpiceKernel::defaults().unwrap();
    /// let epoch = Time::now();
    /// let rock = SpaceRock::from_spice("Earth", &epoch, "J2000", "SSB", &kernel).unwrap();
    /// ```
    pub fn from_spice(name: &str, epoch: &Time, reference_plane: &str, origin: &str, kernel: &SpiceKernel) -> Result<Self, Box<dyn std::error::Error>> {

        // check a priori if the name is in the list of loaded kernels

        let reference_plane = ReferencePlane::from_str(reference_plane)?;
        let origin = Origin::from_str(origin)?;

        let target_id = kernel.body_id(name).ok_or_else(|| crate::spice::SpiceError::UnknownBody(name.to_string()))?;
        let origin_id = kernel.body_id(origin.as_str()).ok_or_else(|| crate::spice::SpiceError::UnknownBody(origin.to_string()))?;

        // J2000 state in AU, AU/day
        let [x, y, z, vx, vy, vz] = kernel.state_au(target_id, origin_id, epoch.tdb().jd())?;
        let position = Vector3::new(x, y, z);
        let velocity = Vector3::new(vx, vy, vz);

        // println!("Position: {:?}", position);
        // let mut ep = epoch.clone();
        // let et = spice::str2et(&format!("JD{epoch} UTC", epoch=epoch.utc().jd()));
        // let (state, _) = spice::spkezr(name, et, reference_plane.as_str(), "NONE", &origin.to_string());
        // let position = Vector3::new(state[0], state[1], state[2]) * KM_TO_AU;
        // let velocity = Vector3::new(state[3], state[4], state[5]) * KM_TO_AU * SECONDS_PER_DAY;

        let mut rock = SpaceRock {
            name: name.to_string(),
            position,
            velocity,
            epoch: epoch.clone(),
            // reference_plane: reference_plane.clone(),
            reference_plane: ReferencePlane::from_str("J2000")?,
            origin,
            properties: None,
        };
        rock.change_reference_plane(reference_plane.as_str())?;

        if let Some(m) = MASSES.get(name.to_lowercase().as_str()) { rock.set_mass(*m) };

        Ok(rock)

    }

    /// Instantiate a SpaceRock with random keplerian elements
    ///
    /// # Arguments
    /// * `epoch` - The epoch of the ephemeris
    /// * `reference_plane` - The coordinate reference_plane
    /// * `origin` - The origin of the coordinate system
    ///
    /// # Returns
    /// * [`SpaceRock`] 

    ///
    /// # Example
    /// ```
    /// use spacerocks::SpaceRock;
    /// use spacerocks::Time;
    ///
    /// let epoch = Time::now();
    /// let rock = SpaceRock::random(&epoch, "J2000", "SSB");
    /// ```
    pub fn random(epoch: &Time, reference_plane: &str, origin: &str) -> Result<Self, Box<dyn std::error::Error>> {

        let mut rng = rand::thread_rng();
        let q = rng.gen_range(2.0..50.0);
        let e = rng.gen_range(0.0..1.5);
        let inc = rng.gen_range(0.0..std::f64::consts::PI);
        let arg = rng.gen_range(0.0..2.0 * std::f64::consts::PI);
        let node = rng.gen_range(0.0..2.0 * std::f64::consts::PI);

        // let max_true_anomaly = ((-1.0 / e) as f64).acos();
        let mut max_true_anomaly = 2.0 * std::f64::consts::PI;
        if e > 1.0 {
            max_true_anomaly = ((-1.0 / e) as f64).acos();
        }

        let true_anomaly = rng.gen_range(-max_true_anomaly..max_true_anomaly);

        // let name = format!("{}", uuid::Uuid::new_v4().simple());
        let name = format!("{}", generate_name(2, 4));


        let rock = SpaceRock::from_kepler(&name, q, e, inc, arg, node, true_anomaly, epoch.clone(), reference_plane, origin)?;
        Ok(rock)
    }

    /// Instantiate a SpaceRock from cartesian coordinates
    ///
    /// # Arguments
    /// * `name` - The name of the object
    /// * `x` - The x-coordinate of the object (au)
    /// * `y` - The y-coordinate of the object (au)
    /// * `z` - The z-coordinate of the object (au)
    /// * `vx` - The x-component of the velocity (au/day)
    /// * `vy` - The y-component of the velocity (au/day)
    /// * `vz` - The z-component of the velocity (au/day)
    /// * `epoch` - The epoch of the ephemeris
    /// * `reference_plane` - The coordinate reference_plane 
    /// * `origin` - The origin of the coordinate system
    ///
    /// # Returns
    /// * [`SpaceRock`] 
    ///
    /// # Example
    /// ```
    /// use spacerocks::SpaceRock;
    /// use spacerocks::Time;
    ///
    /// let epoch = Time::now();
    /// let rock = SpaceRock::from_xyz("Arrokoth", 43.0, 0.0, 0.0, 0.0, 0.0, 0.0, epoch, "J2000", "SSB");
    /// ```
    pub fn from_xyz(name: &str, x: f64, y: f64, z: f64, vx: f64, vy: f64, vz: f64, epoch: Time, reference_plane: &str, origin: &str) -> Result<Self, Box<dyn std::error::Error>> {

        let reference_plane = ReferencePlane::from_str(reference_plane)?;
        let origin = Origin::from_str(origin)?;

        let position = Vector3::new(x, y, z);
        let velocity = Vector3::new(vx, vy, vz);
        let rock = SpaceRock {
                name: name.to_string(),
                position,
                velocity,
                epoch,
                reference_plane: reference_plane.clone(),
                origin: origin.clone(),
                properties: None,
        };

        Ok(rock)
    }

    /// Get a SpaceRock from the JPL Horizons API
    ///
    /// # Arguments
    /// * `name` - The name of the object 
    /// * `epoch` - The epoch of the ephemeris
    /// * `reference_plane` - The coordinate reference_plane 
    /// * `origin` - The origin of the coordinate system
    ///
    /// # Returns
    /// * [`SpaceRock`] 
    ///
    /// # Example
    /// ```
    /// use spacerocks::SpaceRock;
    /// use spacerocks::Time;
    ///
    /// let epoch = Time::now();
    /// let rock = SpaceRock::from_horizons("Arrokoth", &epoch, "J2000", "SSB");
    /// ```
    pub fn from_horizons(name: &str, epoch: &Time, reference_plane: &str, origin: &str) -> Result<Self, Box<dyn std::error::Error>> {

        // let client = reqwest::blocking::Client::new();

        let client = reqwest::blocking::Client::builder()
            .timeout(Duration::from_secs(60)) // wait up to 60s
            .build()?;

        let mut params = HashMap::new();

        let command_str = format!("'{}'", name);
        params.insert("command", command_str.as_str());

        let ep = epoch.clone();

        // let timescale = &ep.timescale.to_str().to_uppercase();
        // let timeformat = &ep.format.to_str().to_uppercase();


        match reference_plane.to_uppercase().as_str() {
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


        // if timescale == "UTC" {
        //     params.insert("TIME_TYPE", "'UT'");
        // } else {
        //     params.insert("TIME_TYPE", timescale);
        // }

        let _timescale = "TDB";
        let timeformat = "JD"; // 'CALENDAR' or 'ISO'

        // ep.to_tdb();
        params.insert("TIME_TYPE", "'TDB'");

        let time_list = format!("'{}'", ep.tdb().jd());
        params.insert("TLIST", time_list.as_str());

        let tf = format!("'{}'", timeformat);
        params.insert("TLIST_TYPE", tf.as_str());

        let center = format!("'@{}'", origin);
        params.insert("center", center.as_str());

        // params.insert("make_ephem", "'yes'");
        params.insert("ephem_type", "'vectors'");
        params.insert("vec_corr", "'None'");
        params.insert("out_units", "'AU-D'");
        params.insert("csv_format", "'yes'");
        params.insert("vec_delta_t", "'no'");
        // params.insert("vec_table", "'2x'");
        params.insert("vec_table", "'2'");
        params.insert("vec_labels", "'no'");

        let response = client.get("https://ssd.jpl.nasa.gov/api/horizons.api?")
            .query(&params)
            .send()?;


        let json: serde_json::Value = response.json()?;
        let text = json["result"].as_str();

        let lines: Vec<&str> = text.ok_or("No data")?.split('\n').collect();

        let first_data_line = lines.iter().skip_while(|&line| !line.starts_with("$$SOE")).nth(1).ok_or("No data")?;
        
        let data: Vec<f64> = first_data_line.split(',').filter_map(|s| s.trim().parse::<f64>().ok()).collect();
        // let given_epoch = Time::new(data[0], "tdb", "jd")?;
        // println!("{:?}", data);
        let (x, y, z, vx, vy, vz) = (data[1], data[2], data[3], data[4], data[5], data[6]);

        let rock = SpaceRock::from_xyz(name, x, y, z, vx, vy, vz, epoch.clone(), reference_plane, origin)?;
        // let rock = SpaceRock::from_xyz(name, x, y, z, vx, vy, vz, given_epoch.clone(), reference_plane, origin)?;
        Ok(rock)
    }

    /// Instantiate a SpaceRock from spherical coordinates (Napier and Holman 2024)
    ///
    /// # Arguments
    /// * `name` - The name of the object
    /// * `phi` - Longitude (radians)
    /// * `theta` - Latutude (radians)
    /// * `r` - Distance from the origin (au)
    /// * `vr` - Radial velocity (au/day)
    /// * `vo` - Tangential velocity (au/day)
    /// * `psi` - Angle between the radial and tangential velocities (radians)
    /// * `epoch` - The epoch of the ephemeris
    /// * `reference_plane` - The coordinate reference_plane 
    /// * `origin` - The origin of the coordinate system
    ///
    /// # Returns
    /// * [`SpaceRock`] 
    pub fn from_spherical(name: &str, phi: f64, theta: f64, r: f64, vr: f64, vo: f64, psi: f64, epoch: Time, reference_plane: &str, origin: &str) -> Result<Self, Box<dyn std::error::Error>> {

        let pointing = Vector3::new(phi.cos() * theta.cos(), phi.sin() * theta.cos(), theta.sin());
        let position = pointing * r;

        let dhat = Vector3::new(-phi.cos() * theta.sin(), -phi.sin() * theta.sin(), theta.cos());
        let ahat = Vector3::new(-phi.sin(), phi.cos(), 0.0);
        let velocity = pointing * vr + vo * (psi.cos() * ahat + psi.sin() * dhat);

        let x = position.x;
        let y = position.y;
        let z = position.z;
        let vx = velocity.x;
        let vy = velocity.y;
        let vz = velocity.z;

        let rock = SpaceRock::from_xyz(name, x, y, z, vx, vy, vz, epoch, reference_plane, origin)?;
        Ok(rock)
    }


    /// Instantiate a SpaceRock from keplerian elements
    ///
    /// # Arguments
    /// * `name` - The name of the object
    /// * `q` - Perihelion distance (au)
    /// * `e` - Eccentricity
    /// * `inc` - Inclination (radians)
    /// * `arg` - Argument of perihelion (radians)
    /// * `node` - Longitude of the ascending node (radians)
    /// * `true_anomaly` - True anomaly (radians)
    /// * `epoch` - The epoch of the ephemeris
    /// * `reference_plane` - The coordinate reference_plane 
    /// * `origin` - The origin of the coordinate system
    ///
    /// # Returns
    /// * [`SpaceRock`] 
    /// # Errors
    /// Returns an error if:
    /// * The true anomaly is not compatible with the eccentricity for hyperbolic orbits
    /// * The reference plane string is invalid
    /// * The origin string is invalid
    pub fn from_kepler(name: &str, q: f64, e: f64, inc: f64, arg: f64, node: f64, true_anomaly: f64, epoch: Time, reference_plane: &str, origin: &str) -> Result<Self, Box<dyn std::error::Error>> {

        let mu = Origin::from_str(origin)?.mu();
        let [x, y, z, vx, vy, vz] = state::from_kepler(q, e, inc, arg, node, true_anomaly, mu)?;

        let rock = SpaceRock::from_xyz(name, x, y, z, vx, vy, vz, epoch, reference_plane, origin)?;
        Ok(rock)
    }


    /// Numerically propagate the SpaceRock to `epoch` with IAS15, in the field of the Sun,
    /// planets, Moon, Pluto and the 16 most massive asteroids (all from `kernel`), with ASSIST's
    /// force model: Newtonian gravity, relativity (Einstein-Infeld-Hoffmann, Sun), Earth J2-J4,
    /// solar J2, and non-gravitational forces if the rock has A1/A2/A3 (`set_nongrav`).
    ///
    /// The integration is done in J2000 about the solar system barycenter; the result is
    /// returned in the rock's original reference plane and origin (SUN or SSB; a custom
    /// origin is returned as SSB).
    pub fn propagate(&mut self, epoch: &Time, kernel: &SpiceKernel) -> Result<(), Box<dyn std::error::Error>> {
        if self.epoch.tdb().jd() == epoch.tdb().jd() {
            return Ok(());
        }
        let original_plane = self.reference_plane.clone();
        let original_origin = self.origin.clone();

        // The simulation's perturber states are barycentric J2000.
        let mut rock = self.clone();
        rock.change_reference_plane("J2000")?;
        rock.to_ssb(kernel)?;

        let mut sim = SpiceSimulation::horizons(&rock.epoch, kernel)?;
        sim.add(rock)?;
        sim.integrate(epoch, kernel)?;
        let p = &sim.state.particles[0];

        self.position = p.position;
        self.velocity = p.velocity;
        self.epoch = epoch.clone();
        self.reference_plane = ReferencePlane::J2000;
        self.origin = Origin::ssb();

        if original_origin == Origin::SUN {
            self.to_helio(kernel)?;
        }
        self.change_reference_plane(original_plane.as_str())?;
        Ok(())
    }

    /// Propagate the SpaceRock in time along a keplerian orbit. The operation is performed in place.
    ///
    /// # Arguments
    /// * `epoch` - The epoch to propagate to
    pub fn analytic_propagate(&mut self, epoch: &Time) -> Result<(), Box<dyn std::error::Error>> {

        let dt = epoch.tdb().jd() - self.epoch.tdb().jd();
        let s = state::kepler_step(&self.state(), self.origin.mu(), dt)?;
        self.set_state(&s);
        self.epoch = epoch.clone();

        Ok(())
    }   

    /// Make a new SpaceRock object with the same properties as the original, but propagated to a new epoch.
    ///
    /// # Arguments
    /// * `epoch` - The epoch to propagate to
    ///
    /// # Returns
    /// * [`SpaceRock`] 
    pub fn analytic_at(&self, epoch: &Time) -> Result<SpaceRock, Box<dyn std::error::Error>> {
        let mut rock = self.clone();
        rock.analytic_propagate(epoch)?;
        Ok(rock)
    }
        

    /// Change the reference plane of the SpaceRock
    ///
    /// # Arguments
    /// * `reference_plane` - The new reference plane
    pub fn change_reference_plane(&mut self, reference_plane: &str) -> Result<(), Box<dyn std::error::Error>> {

        let reference_plane = ReferencePlane::from_str(reference_plane)?;
        if reference_plane == self.reference_plane {
            return Ok(());
        }

        let rot = state::rotation(&self.reference_plane, &reference_plane)?;
        self.position = rot * self.position;
        self.velocity = rot * self.velocity;
        self.reference_plane = reference_plane;

        Ok(())
    }

    /// Change the origin of the SpaceRock
    ///
    /// # Arguments
    /// * `origin` - The SpaceRock object to change the origin to
    pub fn change_origin(&mut self, origin: &SpaceRock) {

        let origin_position = origin.position;
        let origin_velocity = origin.velocity;

        self.position -= origin_position;
        self.velocity -= origin_velocity;

        self.origin = Origin::new_custom(origin.mass() * GRAVITATIONAL_CONSTANT, &origin.name);
    }

    /// Change the origin of the SpaceRock to the solar system barycenter
    /// 
    /// # Example
    /// ```no_run
    /// use spacerocks::{SpaceRock, SpiceKernel, Time};
    ///
    /// let kernel = SpiceKernel::defaults().unwrap();
    /// let epoch = Time::now();
    /// let mut rock = SpaceRock::from_horizons("Arrokoth", &epoch, "J2000", "SSB").unwrap();
    /// rock.to_ssb(&kernel).unwrap();
    /// ```
    pub fn to_ssb(&mut self, kernel: &SpiceKernel) -> Result<(), Box<dyn std::error::Error>> {
        if self.origin == Origin::SSB {
            return Ok(());
        }
        // get the ssb from spice
        // let mut ssb = SpaceRock::from_spice("ssb", &self.epoch, self.reference_plane.as_str(), self.origin.as_str(), &kernel)?;
        let ssb = SpaceRock::from_spice(self.origin.as_str(), &self.epoch, self.reference_plane.as_str(), "ssb", &kernel)?;
        self.position += ssb.position;
        self.velocity += ssb.velocity;

        self.origin = Origin::ssb();

        // set the mass of the ssb to the mass of the bary

        // ssb.set_mass(MU_BARY / GRAVITATIONAL_CONSTANT);
        // self.change_origin(&ssb);
        Ok(())
    }

    /// Change the origin of the SpaceRock to the heliocenter
    ///
    /// # Example
    /// ```no_run
    /// use spacerocks::{SpaceRock, SpiceKernel, Time};
    ///
    /// let kernel = SpiceKernel::defaults().unwrap();
    /// let epoch = Time::now();
    /// let mut rock = SpaceRock::from_horizons("Arrokoth", &epoch, "J2000", "SSB").unwrap();
    /// rock.to_helio(&kernel).unwrap();
    /// ```
    pub fn to_helio(&mut self, kernel: &SpiceKernel) -> Result<(), Box<dyn std::error::Error>> {
        // get the sun from spice
        if self.origin == Origin::SUN {
            return Ok(());
        }
        let sun = SpaceRock::from_spice("sun", &self.epoch, self.reference_plane.as_str(), self.origin.as_str(), &kernel)?;
        self.position -= sun.position;
        self.velocity -= sun.velocity;
        self.origin = Origin::sun();
        Ok(())
    }

    pub fn r_squared(&self) -> f64 {
        self.position.dot(&self.position)
    }

    pub fn v_squared(&self) -> f64 {
        self.velocity.dot(&self.velocity)
    }

    pub fn v(&self) -> f64 {
        self.velocity.norm()
    }

    /// Set the mass in solar masses
    pub fn set_mass(&mut self, mass: f64) {
        if self.properties.is_none() {
            self.properties = Some(Properties::default());
        }
        self.properties.as_mut().unwrap().mass = Some(mass);
    }

    pub fn mass(&self) -> f64 {
        match &self.properties {
            Some(p) => p.mass.unwrap_or(0.0),
            None => 0.0,
        }
    }

    pub fn absolute_magnitude(&self) -> f64 {
        match &self.properties {
            Some(p) => p.absolute_magnitude.unwrap_or(0.0),
            None => 0.0,
        }
    }

    pub fn gslope(&self) -> f64 {
        match &self.properties {
            Some(p) => p.gslope.unwrap_or(0.15),
            None => 0.15,
        }
    }

    pub fn radius(&self) -> f64 {
        match &self.properties {
            Some(p) => p.radius.unwrap_or(0.0),
            None => 0.0,
        }
    }

    pub fn albedo(&self) -> f64 {
        match &self.properties {
            Some(p) => p.albedo.unwrap_or(0.0),
            None => 0.0,
        }
    }

    /// Set the absolute magnitude (H)
    pub fn set_absolute_magnitude(&mut self, absolute_magnitude: f64) {
        if self.properties.is_none() {
            self.properties = Some(Properties::default());
        }
        self.properties.as_mut().unwrap().absolute_magnitude = Some(absolute_magnitude);
        self.properties.as_mut().unwrap().gslope = Some(0.15);
    }

    /// Set the phase slope parameter (G)
    pub fn set_gslope(&mut self, gslope: f64) {
        if self.properties.is_none() {
            self.properties = Some(Properties::default());
        }
        self.properties.as_mut().unwrap().gslope = Some(gslope);
    }

    pub fn set_radius(&mut self, radius: f64) {
        if self.properties.is_none() {
            self.properties = Some(Properties::default());
        }
        self.properties.as_mut().unwrap().radius = Some(radius);
    }

    pub fn set_albedo(&mut self, albedo: f64) {
        if self.properties.is_none() {
            self.properties = Some(Properties::default());
        }
        self.properties.as_mut().unwrap().albedo = Some(albedo);
    }

    /// Set the non-gravitational parameters (A1, A2, A3) of the Marsden model, in AU/day^2
    /// (radial, transverse and normal components at 1 AU), as listed by JPL's small-body
    /// database. They are used by the `NonGravitational` force in numerical propagation.
    pub fn set_nongrav(&mut self, a1: f64, a2: f64, a3: f64) {
        if self.properties.is_none() {
            self.properties = Some(Properties::default());
        }
        self.properties.as_mut().unwrap().nongrav = Some([a1, a2, a3]);
    }

    /// Non-gravitational parameters (A1, A2, A3), if set.
    pub fn nongrav(&self) -> Option<[f64; 3]> {
        self.properties.as_ref().and_then(|p| p.nongrav)
    }

    pub fn r(&self) -> f64 {
        self.position.norm()
    }

    pub fn hvec(&self) -> Vector3<f64> {
        state::angular_momentum(&self.state())
    }

    pub fn nvec(&self) -> Vector3<f64> {
        let hvec = self.hvec();
        Vector3::new(-hvec.y, hvec.x, 0.0)
    }

    pub fn h(&self) -> f64 {
        self.hvec().norm()
    }

    pub fn evec(&self) -> Vector3<f64> {
        state::eccentricity_vector(&self.state(), self.origin.mu())
    }

    /// Calculate the eccentricity (dimensionless)
    pub fn e(&self) -> f64 {
        state::eccentricity(&self.state(), self.origin.mu())
    }

    pub fn specific_energy(&self) -> f64 {
        state::specific_energy(&self.state(), self.origin.mu())
    }

    /// Calculate the semi-major axis in AU
    pub fn a(&self) -> f64 {
        state::semi_major_axis(&self.state(), self.origin.mu())
    }

    /// Calculate the periapsis distance in AU
    pub fn q(&self) -> f64 {
        state::perihelion(&self.state(), self.origin.mu())
    }

    /// Calculate the semi-latus rectum in AU
    pub fn p(&self) -> f64 {
        self.q() * (1.0 + self.e())
    }

    /// Calculate the inclination in radians
    pub fn inc(&self) -> f64 {
        state::inclination(&self.state())
    }

    /// Calculate the argument of perihelion in radians (0 for circular orbits; for equatorial
    /// orbits this is the longitude of perihelion)
    pub fn arg(&self) -> f64 {
        state::argument_of_perihelion(&self.state(), self.origin.mu())
    }

    /// Calculate the longitude of the ascending node in radians (0 for equatorial orbits)
    pub fn node(&self) -> f64 {
        state::node(&self.state())
    }

    /// Calculate the true anomaly in radians. For circular orbits (where perihelion is
    /// undefined) this is the argument of latitude, consistent with `arg() == 0`.
    pub fn true_anomaly(&self) -> f64 {
        state::true_anomaly(&self.state(), self.origin.mu())
    }

    pub fn mean_anomaly(&self) -> f64 {
        state::mean_anomaly(&self.state(), self.origin.mu()).expect("Invalid eccentricity")
    }

    pub fn conic_anomaly(&self) -> f64 {
        state::conic_anomaly(&self.state(), self.origin.mu()).expect("Invalid eccentricity")
    }

    /// All osculating elements in one pass (cheaper than calling the individual methods when
    /// more than one is needed). Anomalies that cannot be computed are NaN.
    pub fn elements(&self) -> Elements {
        state::elements(&self.state(), self.origin.mu())
    }

    /// The state vector `[x, y, z, vx, vy, vz]` (AU, AU/day).
    pub fn state(&self) -> State {
        state::from_pv(&self.position, &self.velocity)
    }

    /// Replace the position and velocity, keeping the epoch, frame and origin.
    pub fn set_state(&mut self, s: &State) {
        let (r, v) = state::pv(s);
        self.position = r;
        self.velocity = v;
    }

    // calculate the osculating elements and return a KeplerOrbit object. This is more expensive than the other 
    // individual methods, but cheaper if you need multiple elements
    // pub fn calculate_orbit(&self) -> KeplerOrbit {
    //     OrbitType::from_eccentricity(e, 1e-10).expect("Invalid eccentricity");
    // }

    /// Compute observational quantities for this object as seen by an observer
    /// 
    /// Calculates topocentric coordinates including light-time correction. The observer
    /// must have the same epoch and reference plane as the SpaceRock.
    /// 
    /// # Arguments
    /// * `observer` - The Observer object representing the viewing location
    /// 
    /// # Returns
    /// * [`Observation`] 
    /// 
    /// # Errors
    /// Returns an error if:
    /// * Observer and SpaceRock have different epochs
    /// * Observer and SpaceRock have different reference planes
    pub fn calc_radec(&self, observer: &Observer) -> Result<(f64, f64), Box<dyn std::error::Error>> {

        // self.change_reference_plane("J2000")?;

        // throw an error if the observer and self have different epochs
        // if self.epoch.utc().jd() != observer.epoch.utc().jd() {
        //     return Err("Observer and SpaceRock have different epochs".into());
        // }

        // if self.reference_plane != observer.reference_plane {
        //     return Err("Observer and SpaceRock have different reference planes".into());
        // }
        // Calculate the topocentric state, correct for light travel time
        let cr = correct_for_ltt(&self, observer);

        // Calaculate the ra, and dec
        let mut ra = cr.position.y.atan2(cr.position.x);
        if ra < 0.0 {
            ra += 2.0 * std::f64::consts::PI;
        }
        let dec = (cr.position.z / cr.position.norm()).asin();

        Ok((ra, dec))
        
    }


    pub fn calc_radec_no_light_time(&self, observer: &Observer) -> Result<(f64, f64), Box<dyn std::error::Error>> {

        // self.change_reference_plane("J2000")?;

        // throw an error if the observer and self have different epochs
        // if self.epoch.utc().jd() != observer.epoch.utc().jd() {
        //     return Err("Observer and SpaceRock have different epochs".into());
        // }

        // if self.reference_plane != observer.reference_plane {
        //     return Err("Observer and SpaceRock have different reference planes".into());
        // }
        // Calculate the topocentric state, correct for light travel time
        let cr = self.position - observer.position;

        // Calaculate the ra, and dec
        let mut ra = cr.y.atan2(cr.x);
        if ra < 0.0 {
            ra += 2.0 * std::f64::consts::PI;
        }
        let dec = (cr.z / cr.norm()).asin();

        Ok((ra, dec))
        
    }


    /// Compute observational quantities for this object as seen by an observer
    /// 
    /// Calculates topocentric coordinates including light-time correction. The observer
    /// must have the same epoch and reference plane as the SpaceRock.
    /// 
    /// # Arguments
    /// * `observer` - The Observer object representing the viewing location
    /// 
    /// # Returns
    /// * [`Observation`] 
    /// 
    /// # Errors
    /// Returns an error if:
    /// * Observer and SpaceRock have different epochs
    /// * Observer and SpaceRock have different reference planes
    pub fn observe(&mut self, observer: &Observer) -> Result<Observation, Box<dyn std::error::Error>> {

        // self.change_reference_plane("J2000")?;

        // throw an error if the observer and self have different epochs
        if self.epoch.utc().jd() != observer.epoch.utc().jd() {
            return Err("Observer and SpaceRock have different epochs".into());
        }

        if self.reference_plane != observer.reference_plane {
            return Err("Observer and SpaceRock have different reference planes".into());
        }
        if observer.velocity.is_none() {
            return Err("Observer velocity is required to compute rates".into());
        }
        let app = self.apparent_unchecked(observer);
        let has_h = self.properties.as_ref().and_then(|p| p.absolute_magnitude).is_some();
        let mag = if has_h { Some(app.mag) } else { None };
        let observation = Observation::from_complete(self.epoch.clone(), app.ra, app.dec, app.ra_rate, app.dec_rate, app.range, app.range_rate, observer.clone(), None, mag, None)?;
        Ok(observation)
    }

    /// Observable quantities (RA/Dec, rates, range, phase, elongation, magnitude) as plain
    /// numbers, with the same light-time correction as [`SpaceRock::observe`] but without
    /// building an [`Observation`]. The observer must be at the rock's epoch (to within 1 µs),
    /// in the same reference plane, and have a velocity. Heliocentric distance, phase and
    /// elongation use the observer's `sun_position` (the origin if it has none).
    pub fn apparent(&self, observer: &Observer) -> Result<Apparent, Box<dyn std::error::Error>> {
        if observer.velocity.is_none() {
            return Err("Observer velocity is required to compute rates".into());
        }
        if self.reference_plane != observer.reference_plane {
            return Err("Observer and SpaceRock have different reference planes".into());
        }
        if (self.epoch.tdb().jd() - observer.epoch.tdb().jd()).abs() > 1e-6 / 86400.0 {
            return Err("Observer and SpaceRock have different epochs".into());
        }
        Ok(self.apparent_unchecked(observer))
    }

    fn apparent_unchecked(&self, observer: &Observer) -> Apparent {
        let obs_vel = observer.velocity.unwrap_or_else(Vector3::zeros);
        let (h, g) = match &self.properties {
            Some(p) => (p.absolute_magnitude, p.gslope.unwrap_or(0.15)),
            None => (None, 0.15),
        };
        apparent(&self.position, &self.velocity, &observer.position, &obs_vel, &observer.sun(), h, g)
    }

}

/// Apparent magnitude in the IAU H-G system (Bowell et al. 1989).
///
/// * `h`, `g` - absolute magnitude and slope parameter
/// * `r`, `delta` - heliocentric and observer distances (AU)
/// * `phase` - Sun-object-observer angle (radians)
pub fn hg_magnitude(h: f64, g: f64, r: f64, delta: f64, phase: f64) -> f64 {
    let t = (phase / 2.0).tan();
    let phi1 = (-3.332 * t.powf(0.631)).exp();
    let phi2 = (-1.862 * t.powf(1.218)).exp();
    let reduced = h + 5.0 * (r * delta).log10();
    let phase_term = (1.0 - g) * phi1 + g * phi2;
    if phase_term > 0.0 {
        reduced - 2.5 * phase_term.log10()
    } else {
        reduced
    }
}

/// Display the SpaceRock object with each field on a new line
impl std::fmt::Display for SpaceRock {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "SpaceRock: {}\nEpoch: {:?}\nReference Plane: {}\nOrigin: {}\nPosition: {:?}\nVelocity: {:?}\nProperties: {:?}", 
        self.name, self.epoch, self.reference_plane, self.origin, self.position, self.velocity, self.properties)
    }
}



/// Generate a "name" made up of random syllables.
/// `min_syllables` - minimum number of syllables
/// `max_syllables` - maximum number of syllables
fn generate_name(min_syllables: usize, max_syllables: usize) -> String {
    // Choose how many syllables this name will have
    let syllables_count = rand::thread_rng().gen_range(min_syllables..=max_syllables);

    // Build the name
    let mut name = String::new();
    for _i in 0..syllables_count {
        let s = generate_syllable();
        name.push_str(&s);
    }

    name
}

/// Generate a single syllable in the form:
/// (optional consonant) + vowel + (optional consonant)
fn generate_syllable() -> String {
    let vowels = ['a', 'e', 'i', 'o', 'u'];
    // You can include more consonants if you like.
    // Some letters (like 'q', 'x', 'z') might produce more unusual results.
    let consonants = [
        'b', 'c', 'd', 'f', 'g', 'h', 'j', 'k', 'l', 'm',
        'n', 'p', 'r', 's', 't', 'v', 'w', 'y', 'z', 'x', 'q',
    ];
    
    let mut rng = rand::thread_rng();
    
    let mut syllable = String::new();
    
    // 50% chance to start with a consonant
    if rng.gen_bool(0.5) {
        syllable.push(consonants[rng.gen_range(0..consonants.len())]);
    }
    
    // Always include a vowel
    syllable.push(vowels[rng.gen_range(0..vowels.len())]);
    
    // 40% chance to add a trailing consonant
    if rng.gen_bool(0.4) {
        syllable.push(consonants[rng.gen_range(0..consonants.len())]);
    }
    
    syllable
}