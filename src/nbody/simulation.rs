//! Set up and manage a collection of gravitationally interacting particles.
use std::collections::HashMap;

use crate::SpaceRock;
use crate::constants::GRAVITATIONAL_CONSTANT;
use crate::time::Time;
use crate::{Properties, ReferencePlane, Origin};
use crate::errors::SimulationError;
use crate::state::{self, State};


use crate::nbody::forces::{Force, NewtonianGravity};
use crate::nbody::integrators::{Integrator, IAS15};
use crate::spice::SpiceKernel;


use nalgebra::Vector3;

/// A simulation maintains:
/// - The states of its particles as one array of `[x, y, z, vx, vy, vz]`, with their masses,
///   names and properties in columns beside it (massive particles first, by decreasing mass)
/// - The current simulation epoch (TDB), shared by every particle, kept as a reference Julian
///   date plus the days since it (as ASSIST does), so steps of any length add up exactly
/// - Reference frame and origin specifications, shared by every particle
/// - Integration method and forces
///
/// The integrator and forces see only the state array and the masses. Particles go in as
/// [`SpaceRock`]s ([`Simulation::add`]) and come out as `SpaceRock`s built on demand
/// ([`Simulation::get_particle`], [`Simulation::particles`]).
#[derive(Clone)]
pub struct Simulation {
    /// TDB Julian date the simulation's clock counts from.
    jd_ref: f64,
    /// Days since `jd_ref` of every particle; the integrators' clock.
    t: f64,
    pub particle_index_map: HashMap<String, usize>,

    pub reference_plane: ReferencePlane,
    pub origin: Origin,

    states: Vec<State>,
    masses: Vec<f64>,
    names: Vec<String>,
    properties: Vec<Option<Properties>>,

    pub integrator: Box<dyn Integrator + Send + Sync>,
    pub forces: Vec<Box<dyn Force + Send + Sync>>,

    /// The integrator's own state while `particles` hold interpolated states.
    synced: Option<SyncedState>,
    /// Whether the integrator's last step describes the current particles.
    dense_valid: bool,
}

/// Where the integrator actually is, kept aside by [`Simulation::integrate_or_interpolate`].
#[derive(Clone)]
struct SyncedState {
    t: f64,
    states: Vec<State>,
}

impl Default for Simulation {
    fn default() -> Self {
        Self::new(&Time::now(), "J2000", "SSB").unwrap()
    }
}

impl Simulation {
    /// Creates a new simulation at the specified epoch (kept as a TDB Julian date) and
    /// reference frame.
    pub fn new(epoch: &Time, reference_plane: &str, origin: &str) -> Result<Simulation, Box<dyn std::error::Error>> {
        let reference_plane = ReferencePlane::from_str(reference_plane)?;
        let origin = Origin::from_str(origin)?;

        Ok(Simulation {
            states: Vec::new(),
            masses: Vec::new(),
            names: Vec::new(),
            properties: Vec::new(),
            jd_ref: epoch.tdb().jd(),
            t: 0.0,
            forces: vec![Box::new(NewtonianGravity)],
            reference_plane: reference_plane,  
            origin: origin,
            integrator: Box::new(IAS15::new(1.0)),
            particle_index_map: HashMap::new(),
            synced: None,
            dense_valid: false,
        })
    }

    /// The simulation's epoch (TDB Julian date).
    pub fn epoch(&self) -> Time {
        Time::new(self.jd_ref + self.t, "tdb", "jd").expect("valid time")
    }

    /// Days since the reference epoch, the integrators' clock.
    pub fn t(&self) -> f64 {
        self.t
    }

    /// Set the epoch without moving the particles (it becomes the new reference epoch).
    pub fn set_epoch(&mut self, epoch: &Time) {
        self.synchronize();
        self.dense_valid = false;
        self.jd_ref = epoch.tdb().jd();
        self.t = 0.0;
    }

    /// `epoch` in days since the reference epoch.
    fn days_since_ref(&self, epoch: &Time) -> f64 {
        epoch.tdb().jd() - self.jd_ref
    }

    /// Instantiate a simulation with the solar system giants.
    ///
    /// # Arguments
    ///
    /// * `epoch` - The epoch of the simulation.
    /// * `reference_plane` - The reference plane of the simulation.
    /// * `origin` - The origin of the simulation.
    ///
    /// # Returns
    ///
    /// * `Result<Simulation, Box<dyn std::error::Error>>` - The simulation with the solar system giants.
    pub fn giants(epoch: &Time, reference_plane: &str, origin: &str, kernel: &SpiceKernel) -> Result<Simulation, Box<dyn std::error::Error>> {

        let mut sim = Simulation::new(epoch, reference_plane, origin)?;
        sim.integrator = Box::new(IAS15::new(1.0));

        // add sun, jupiter barycenter, saturn barycenter, uranus barycenter, neptune barycenter.
        for name in ["sun", "jupiter barycenter", "saturn barycenter", "uranus barycenter", "neptune barycenter"].iter() {
            let particle = SpaceRock::from_spice(name, epoch, reference_plane, origin, kernel)?;
            sim.add(particle)?;
        }
        Ok(sim)
    }

    /// Instantiate a simulation with the solar system planets.
    /// Includes the sun, mercury barycenter, venus barycenter, earth barycenter, mars barycenter, jupiter barycenter, saturn barycenter, uranus barycenter, neptune barycenter.
    ///
    /// # Arguments
    ///
    /// * `epoch` - The epoch of the simulation.
    /// * `reference_plane` - The reference plane of the simulation.
    /// * `origin` - The origin of the simulation.
    ///
    /// # Returns
    ///
    /// * `Result<Simulation, Box<dyn std::error::Error>>` - The simulation with the solar system planets.
    pub fn planets(epoch: &Time, reference_plane: &str, origin: &str, kernel: &SpiceKernel) -> Result<Simulation, Box<dyn std::error::Error>> {
        let mut sim = Simulation::new(epoch, reference_plane, origin)?;
        sim.integrator = Box::new(IAS15::new(1.0));

        let names = ["sun", "mercury barycenter", "venus barycenter", "earth barycenter", "mars barycenter", "jupiter barycenter", 
                     "saturn barycenter", "uranus barycenter", "neptune barycenter"];
        for name in names.iter() {
            let particle = SpaceRock::from_spice(name, epoch, reference_plane, origin, kernel)?;
            sim.add(particle)?;
        }
        Ok(sim)
    }

    /// Instantiate a simulation with the solar system planets and moons.
    /// Includes the sun, mercury barycenter, venus barycenter, earth, moon, mars barycenter, jupiter barycenter, saturn barycenter, uranus barycenter, neptune barycenter, pluto barycenter.
    ///
    /// # Arguments
    ///
    /// * `epoch` - The epoch of the simulation.
    /// * `reference_plane` - The reference plane of the simulation.
    /// * `origin` - The origin of the simulation.
    ///
    /// # Returns
    ///
    /// * `Result<Simulation, Box<dyn std::error::Error>>` - The simulation with the solar system planets and moons.
    pub fn horizons(epoch: &Time, reference_plane: &str, origin: &str, kernel: &SpiceKernel) -> Result<Simulation, Box<dyn std::error::Error>> {
        let mut sim = Simulation::new(epoch, reference_plane, origin)?;
        sim.integrator = Box::new(IAS15::new(0.001));

        let names = ["sun", "mercury barycenter", "venus barycenter", "earth", "moon", "mars barycenter", "jupiter barycenter", 
                     "saturn barycenter", "uranus barycenter", "neptune barycenter", "pluto barycenter", 
                     "2000001", 
                     "2000002", 
                     "2000003", 
                     "2000004", 
                     "2000007",
                     "2000010", 
                     "2000015", 
                     "2000016", 
                     "2000031", 
                    //  "2000048", 
                     "2000052", 
                     "2000065", 
                     "2000087",
                     "2000088", 
                     "2000107",
                    //  "2000451", 
                     "2000511", 
                     "2000704"];
        for name in names.iter() {
            let particle = SpaceRock::from_spice(name, epoch, reference_plane, origin, kernel)?;
            sim.add(particle)?;
        }
        Ok(sim)
    }

    
    /// Add a particle to the simulation.
    ///
    /// # Arguments
    ///
    /// * `particle` - The particle to add to the simulation.
    pub fn add(&mut self, mut particle: SpaceRock) -> Result<(), Box<dyn std::error::Error>> {
        self.synchronize();
        self.dense_valid = false;

        if (self.days_since_ref(&particle.epoch) - self.t).abs() > EPOCH_TOLERANCE {
            let err = SimulationError::EpochMismatch(particle.epoch.clone(), self.epoch(), particle.name.clone());
            return Err(err.into());
        }

        if particle.origin.clone() != self.origin {
            if !self.particle_index_map.contains_key(&particle.origin.to_string()) {
                let err = SimulationError::OriginMismatch(particle.origin.clone(), self.origin.clone(), particle.name.clone());
                return Err(err.into());
            }
            let origin = self.particle(self.particle_index_map[&particle.origin.to_string()]);
            particle.change_origin(&origin);
            println!("Changing origin of {} from {} to {}", particle.name, particle.origin, origin.name);
        }

        particle.change_reference_plane(self.reference_plane.as_str())?;

        let massive = particle.mass() != 0.0;
        self.particle_index_map.insert(particle.name.clone(), self.states.len());
        self.states.push(particle.state());
        self.masses.push(particle.mass());
        self.names.push(particle.name);
        self.properties.push(particle.properties);

        if massive {
            // Keep the particles sorted by decreasing mass (a stable sort).
            let mut order: Vec<usize> = (0..self.states.len()).collect();
            order.sort_by(|&a, &b| self.masses[b].partial_cmp(&self.masses[a]).unwrap());
            self.states = order.iter().map(|&i| self.states[i]).collect();
            self.masses = order.iter().map(|&i| self.masses[i]).collect();
            self.names = order.iter().map(|&i| self.names[i].clone()).collect();
            self.properties = order.iter().map(|&i| self.properties[i].clone()).collect();
            for (idx, name) in self.names.iter().enumerate() {
                self.particle_index_map.insert(name.clone(), idx);
            }
        }

        Ok(())
    }

    /// Remove a particle from the simulation.
    ///
    /// # Arguments
    ///
    /// * `name` - The name of the particle to remove.
    pub fn remove(&mut self, name: &str) -> Result<(), SimulationError> {
        self.synchronize();
        self.dense_valid = false;
        if self.particle_index_map.contains_key(name) {
            let idx = self.particle_index_map[name];
            self.states.remove(idx);
            self.masses.remove(idx);
            self.names.remove(idx);
            self.properties.remove(idx);
            self.particle_index_map.remove(name);
            for value in self.particle_index_map.values_mut() {
                if *value > idx {
                    *value -= 1;
                }
            }
        } else {
            return Err(SimulationError::ParticleNotFound(name.to_string()));
        }
        Ok(())
    }

    /// Move the simulation to the center of mass.
    pub fn move_to_center_of_mass(&mut self) -> Result<(), Box<dyn std::error::Error>> {
        self.synchronize();
        self.dense_valid = false;
        let mut total_mass = 0.0;
        let mut center_of_mass = Vector3::new(0.0, 0.0, 0.0);
        let mut center_of_mass_velocity = Vector3::new(0.0, 0.0, 0.0);

        for (s, &m) in self.states.iter().zip(&self.masses) {
            if m == 0.0 {
                continue;
            }
            let (x, v) = state::pv(s);
            center_of_mass += m * x;
            center_of_mass_velocity += m * v;
            total_mass += m;
        }

        center_of_mass /= total_mass;
        center_of_mass_velocity /= total_mass;

        for s in &mut self.states {
            let (x, v) = state::pv(s);
            *s = state::from_pv(&(x - center_of_mass), &(v - center_of_mass_velocity));
        }

        let origin = Origin::new_custom(total_mass * GRAVITATIONAL_CONSTANT, "simulation_barycenter");
        self.origin = origin;
        Ok(())
    }

    /// Change the origin of the simulation.
    ///
    /// # Arguments
    ///     
    /// * `origin` - The name of the particle to set as the origin.
    pub fn change_origin(&mut self, origin: &str) -> Result<(), String> {
        self.synchronize();
        self.dense_valid = false;

        if !self.particle_index_map.contains_key(origin) {
           return Err(format!("Origin {} not found in perturbers", origin));
        }

        let idx = self.particle_index_map[origin];
        self.origin = Origin::new_custom(self.masses[idx] * GRAVITATIONAL_CONSTANT, origin);

        let (origin_position, origin_velocity) = state::pv(&self.states[idx]);
        for s in &mut self.states {
            let (x, v) = state::pv(s);
            *s = state::from_pv(&(x - origin_position), &(v - origin_velocity));
        }

        Ok(())
    }

    /// Step the simulation forward in time by one timestep.
    pub fn step(&mut self) {
        self.synchronize();
        self.integrator.step(&mut self.states, &self.masses, &mut self.t, &self.forces);
        self.dense_valid = true;
    }

    /// Step the simulation forward by `n` timesteps. Same result as calling [`Simulation::step`]
    /// `n` times, but integrators can share work between the steps: Wisdom–Holman then splits
    /// test particles between threads.
    pub fn steps(&mut self, n: usize) {
        if n == 0 {
            return;
        }
        self.synchronize();
        self.integrator.steps(&mut self.states, &self.masses, &mut self.t, &self.forces, n);
        self.dense_valid = true;
    }

    /// Integrate to a new epoch, or interpolate to it when the integrator supports it.
    ///
    /// The integrator takes its own adaptive steps until one of them brackets `epoch`, and the
    /// particles are then interpolated to `epoch` from that step. The integrator never shortens a
    /// step to land on `epoch`, so calling this for a sorted sequence of epochs costs no more than
    /// integrating to the last one, and the result doesn't depend on how many epochs you ask for.
    /// An epoch inside the last step, before or after the current one, needs no new steps.
    ///
    /// The particles then hold the interpolated states, and the integrator's own state is kept
    /// aside. Any later call that steps, integrates, adds or removes particles, or changes the
    /// origin restores it first, so changes made to interpolated particles are discarded.
    ///
    /// Integrators without dense output (only [`IAS15`] has it) fall back to
    /// [`Simulation::integrate`].
    ///
    /// # Arguments
    ///
    /// * `epoch` - The epoch to integrate or interpolate to.
    pub fn integrate_or_interpolate(&mut self, epoch: &Time) {
        if !self.integrator.has_dense_output() {
            self.integrate(epoch);
            return;
        }
        self.synchronize();

        let target = self.days_since_ref(epoch);
        if (target - self.t).abs() < 1e-16 {
            return;
        }

        loop {
            if self.dense_valid {
                if let Some(states) = self.integrator.interpolate(target) {
                    let synced = SyncedState { t: self.t, states: std::mem::replace(&mut self.states, states) };
                    self.t = target;
                    self.synced = Some(synced);
                    return;
                }
            }

            // Step towards the target.
            let dt = target - self.t;
            if self.dense_valid && dt.abs() < 1e-16 {
                return;
            }
            let timestep = self.integrator.timestep();
            if (dt < 0.0) != (timestep < 0.0) {
                self.integrator.set_timestep(-timestep);
            }
            self.step();
            if (target - self.t) * dt < 0.0 && self.integrator.interpolate(target).is_none() {
                // Stepped past the target without being able to interpolate back; finish exactly.
                self.integrate(epoch);
                return;
            }
        }
    }

    /// Restore the integrator's own state if the particles hold interpolated states.
    fn synchronize(&mut self) {
        if let Some(synced) = self.synced.take() {
            self.states = synced.states;
            self.t = synced.t;
        }
    }

    /// Integrate the simulation to a new epoch.
    ///
    /// Advances the simulation time by taking steps until reaching the target epoch.
    /// The integration direction (forward/backward) is determined automatically based
    /// on the difference between current and target epochs. The timestep is adjusted
    /// automatically when approaching the target epoch to hit it exactly.
    ///
    /// # Arguments
    ///
    /// * `epoch` - The new epoch to integrate to.
    pub fn integrate(&mut self, epoch: &Time) {
        self.synchronize();

        let target = self.days_since_ref(epoch);
        let dt = target - self.t;
        if dt.abs() < 1e-16 {
            return;
        }

        if dt < 0.0 && self.integrator.timestep() > 0.0 {
            self.integrator.set_timestep(-self.integrator.timestep());
        }


        loop {
            let dt = target - self.t;

            // done integrating
            if dt.abs() < 1e-16 {
                break;
            }

            // if we're within a timestep of the epoch, just take a step of that size
            if dt.abs() < self.integrator.timestep().abs() {
                // Shorten one step to land on the epoch, then go back to the timestep the
                // integrator wanted (IAS15 predicts its coefficients again for it). If IAS15
                // rejected the short step and took an even shorter one, keep the timestep it
                // chose instead: a rejected step comes back at most a quarter as long.
                let full_timestep = self.integrator.timestep();
                let start = self.t;
                self.integrator.set_timestep(dt);
                self.step();
                let taken = self.t - start;
                if (taken - dt).abs() < 0.5 * dt.abs() {
                    self.integrator.set_timestep(full_timestep);
                }
                continue;

                // let dt = epoch.tdb().jd() - self.epoch.tdb().jd();
                // if dt.abs() < 1e-16 {
                //     self.integrator.set_timestep(last_timestep);
                //     break;
                // }
                
            }

            // if the timestep is negative, make sure the integrator is set to negative
            if dt < 0.0 {
                if self.integrator.timestep() > 0.0 {
                    self.integrator.set_timestep(-self.integrator.timestep());
                }
            } else if self.integrator.timestep() < 0.0 {
                self.integrator.set_timestep(-self.integrator.timestep());
            }
            // A fixed-step integrator takes every full step this loop would take before the last,
            // shorter one, together; an adaptive one may change its timestep at each step.
            let n = if self.integrator.fixed_timestep() { self.full_steps_to(target) } else { 1 };
            self.steps(n.max(1));
        }
        
        // let dt = &epoch - &self.epoch;
        // let dt = epoch.tdb().jd() - self.epoch.tdb().jd();
        // if dt.abs() < 1e-16 {
        //     return;
        // }
        // // create an exact match for the epoch
        // let old_timestep = self.integrator.timestep();
        // self.integrator.set_timestep(dt);
        // self.step();
        // // reset the timestep
        // self.integrator.set_timestep(old_timestep);
    }

    /// How many steps of the current timestep fit before `target` (days since the reference
    /// epoch), as the loop in [`Simulation::integrate`] counts them, on the same arithmetic.
    fn full_steps_to(&self, target: f64) -> usize {
        let h = self.integrator.timestep();
        let mut t = self.t;
        let mut n = 0;
        loop {
            let dt = target - t;
            if dt.abs() < 1e-16 || dt.abs() < h.abs() || (dt < 0.0) != (h < 0.0) {
                return n;
            }
            t += h;
            n += 1;
        }
    }

    /// The `idx`-th particle as a [`SpaceRock`] (at the simulation's epoch, in its frame and
    /// about its origin).
    pub fn particle(&self, idx: usize) -> SpaceRock {
        let (position, velocity) = state::pv(&self.states[idx]);
        SpaceRock {
            name: self.names[idx].clone(),
            epoch: self.epoch(),
            reference_plane: self.reference_plane.clone(),
            origin: self.origin.clone(),
            position,
            velocity,
            properties: self.properties[idx].clone(),
        }
    }

    /// Every particle as a [`SpaceRock`], in the simulation's order (massive particles first).
    pub fn particles(&self) -> Vec<SpaceRock> {
        (0..self.states.len()).map(|i| self.particle(i)).collect()
    }

    /// Get a particle from the simulation by name.
    ///
    /// # Arguments
    ///
    /// * `name` - The name of the particle to get.
    pub fn get_particle(&self, name: &str) -> Result<SpaceRock, SimulationError> {
        match self.particle_index_map.get(name) {
            Some(&idx) => Ok(self.particle(idx)),
            None => Err(SimulationError::ParticleNotFound(name.to_string())),
        }
    }

    /// Number of particles.
    pub fn len(&self) -> usize {
        self.states.len()
    }

    pub fn is_empty(&self) -> bool {
        self.states.is_empty()
    }

    /// States of the particles `[x, y, z, vx, vy, vz]`, in the simulation's order.
    pub fn states(&self) -> &[State] {
        &self.states
    }

    /// Mutable states of the particles. Changing them restarts the integrator's dense output.
    pub fn states_mut(&mut self) -> &mut [State] {
        self.synchronize();
        self.dense_valid = false;
        &mut self.states
    }

    /// Masses of the particles (solar masses), in the simulation's order.
    pub fn masses(&self) -> &[f64] {
        &self.masses
    }

    /// Names of the particles, in the simulation's order.
    pub fn names(&self) -> &[String] {
        &self.names
    }

    /// Get the energy of the simulation.
    pub fn energy(&self) -> f64 {
        let mut kinetic_energy = 0.0;
        let mut potential_energy = 0.0;

        let n = self.states.len();
        for idx in 0..n {
            let (x_i, v_i) = state::pv(&self.states[idx]);
            kinetic_energy += 0.5 * self.masses[idx] * v_i.norm_squared();
            for jdx in (idx + 1)..n {
                let r = (x_i - state::position(&self.states[jdx])).norm();
                potential_energy -= GRAVITATIONAL_CONSTANT * self.masses[idx] * self.masses[jdx] / r;
            }
        }
        kinetic_energy + potential_energy
    }

    /// Add a force to the simulation.
    ///
    /// # Arguments
    ///
    /// * `force` - The force to add to the simulation.
    pub fn add_force(&mut self, force: Box<dyn Force + Send + Sync>) {
        self.forces.push(force);
    }

}

/// How close (days) a particle's epoch must be to the simulation's to be added: about two
/// rounding steps of a Julian date.
const EPOCH_TOLERANCE: f64 = 1e-9;
