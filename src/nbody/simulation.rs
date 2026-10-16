//! Set up and manage a collection of gravitationally interacting particles.
use std::collections::HashMap;

use crate::SpaceRock;
use crate::constants::GRAVITATIONAL_CONSTANT;
use crate::time::Time;
use crate::{Properties, ReferencePlane, Origin};
use crate::errors::SimulationError;
use crate::state::{self, State};


use crate::nbody::forces::{central_body, Force, NewtonianGravity, SolarGR};
use crate::nbody::integrators::{Integrator, Leapfrog, Trace, WisdomHolman, IAS15};
use crate::nbody::monitor::{Command, ProgressMonitor, Release};
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
    /// The live progress page, if one is attached ([`Simulation::monitor`]).
    pub(crate) progress_monitor: Option<ProgressMonitor>,
    /// Asked between steps whether to stop ([`Simulation::set_interrupt`]), when it was last
    /// asked, and whether it said to.
    interrupt: Option<InterruptCheck>,
    interrupt_asked: Option<std::time::Instant>,
    interrupted: bool,
    /// A new end (days since the reference epoch, infinite for none) for the running
    /// `integrate`, `integrate_or_interpolate` or `run`, set on the progress page.
    end_override: Option<f64>,
    /// Length (days) of the last step taken, for the progress page's Back.
    last_step: f64,
}

/// A check for whether to stop integrating, e.g. whether Ctrl-C was pressed. See
/// [`Simulation::set_interrupt`].
pub type InterruptCheck = std::sync::Arc<dyn Fn() -> bool + Send + Sync>;

/// What the stepping loops do after a [`Simulation::checkpoint`].
enum Go {
    /// Step on (up to [`MONITOR_CHUNK`] steps before the next checkpoint).
    Steps,
    /// Take one step: the progress page's Step button let it through.
    OneStep,
    /// Stop: the interrupt check said so.
    Stop,
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
            progress_monitor: None,
            interrupt: None,
            interrupt_asked: None,
            interrupted: false,
            end_override: None,
            last_step: 0.0,
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

    /// The reference epoch (TDB Julian date) the clock counts from.
    pub(crate) fn jd_ref(&self) -> f64 {
        self.jd_ref
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
        if !matches!(self.checkpoint(), Go::Stop) {
            self.step_unmonitored();
        }
    }

    /// Step the simulation forward by `n` timesteps. Same result as calling [`Simulation::step`]
    /// `n` times, but integrators can share work between the steps: Wisdom–Holman then splits
    /// test particles between threads.
    pub fn steps(&mut self, n: usize) {
        if !self.has_checkpoints() {
            self.steps_unmonitored(n);
            return;
        }
        // Hand the steps over in chunks, so the page updates (and the run can pause or stop)
        // during a long call.
        let mut left = n;
        while left > 0 {
            let chunk = match self.checkpoint() {
                Go::Stop => return,
                Go::OneStep => 1,
                Go::Steps => left.min(MONITOR_CHUNK),
            };
            self.steps_unmonitored(chunk);
            left -= chunk;
        }
    }

    /// One step, with no stop at the progress page's controls.
    fn step_unmonitored(&mut self) {
        self.synchronize();
        let start = self.t;
        self.integrator.step(&mut self.states, &self.masses, &mut self.t, &self.forces);
        self.last_step = self.t - start;
        self.dense_valid = true;
        if let Some(monitor) = &self.progress_monitor {
            monitor.tick(self, 1);
        }
    }

    /// `n` steps, with no stop at the progress page's controls.
    fn steps_unmonitored(&mut self, n: usize) {
        if n == 0 {
            return;
        }
        self.synchronize();
        let start = self.t;
        self.integrator.steps(&mut self.states, &self.masses, &mut self.t, &self.forces, n);
        self.last_step = (self.t - start) / n as f64;
        self.dense_valid = true;
        if let Some(monitor) = &self.progress_monitor {
            monitor.tick(self, n);
        }
    }

    /// Ask `check` between steps (at most every 50 ms) whether to stop. When it says yes,
    /// `step`, `steps`, `integrate` and `integrate_or_interpolate` return at once, between two
    /// steps (also while the progress page holds the run paused), and [`Simulation::interrupted`]
    /// reports it. The Python bindings use this to stop on Ctrl-C. `None` removes the check.
    pub fn set_interrupt(&mut self, check: Option<InterruptCheck>) {
        self.interrupt = check;
        self.interrupt_asked = None;
        self.interrupted = false;
    }

    /// Whether the interrupt check stopped the simulation since it was set.
    pub fn interrupted(&self) -> bool {
        self.interrupted
    }

    /// Whether the stepping loops need to stop between steps (for a progress page or an
    /// interrupt check), rather than hand the integrator every step at once.
    fn has_checkpoints(&self) -> bool {
        self.progress_monitor.is_some() || self.interrupt.is_some()
    }

    /// Ask the interrupt check, if one is set and it wasn't asked in the last 50 ms.
    fn should_stop(&mut self) -> bool {
        if self.interrupted {
            return true;
        }
        let Some(check) = &self.interrupt else { return false };
        let now = std::time::Instant::now();
        if self.interrupt_asked.is_some_and(|last| now - last < INTERRUPT_INTERVAL) {
            return false;
        }
        self.interrupt_asked = Some(now);
        self.interrupted = check();
        self.interrupted
    }

    /// Stop between steps: for the interrupt check, and to wait while the progress page is
    /// paused, applying its commands as they come.
    fn checkpoint(&mut self) -> Go {
        if self.should_stop() {
            return Go::Stop;
        }
        let Some(monitor) = self.progress_monitor.clone() else { return Go::Steps };
        loop {
            let (commands, released) = monitor.wait(self);
            if !commands.is_empty() {
                for command in commands {
                    self.apply_command(command);
                }
                monitor.publish(self);
            }
            match released {
                Some(Release::One) => return Go::OneStep,
                Some(Release::Run) => return Go::Steps,
                Some(Release::Back) => self.step_back(),
                None if self.should_stop() => return Go::Stop,
                None => {}
            }
        }
    }

    /// One step against the direction the simulation is going in, for the progress page's Back:
    /// as long as the last step (so an adaptive integrator goes back about where it was), then
    /// the timestep is as it was before.
    fn step_back(&mut self) {
        let h = self.integrator.timestep();
        let back = if self.last_step != 0.0 { self.last_step.abs() } else { h.abs() };
        self.synchronize();
        self.integrator.set_timestep(-back * h.signum());
        self.step_unmonitored();
        self.integrator.set_timestep(h);
        self.dense_valid = false;
    }

    /// Apply a command from the progress page. Timesteps keep the direction the simulation is
    /// going in.
    fn apply_command(&mut self, command: Command) {
        let sign = if self.integrator.timestep() < 0.0 { -1.0 } else { 1.0 };
        match command {
            Command::Integrator(name, h) => {
                let h = sign * h.abs();
                let integrator: Box<dyn Integrator + Send + Sync> = match name.as_str() {
                    "IAS15" => Box::new(IAS15::new(h)),
                    "WisdomHolman" => Box::new(WisdomHolman::new(h)),
                    "Leapfrog" => Box::new(Leapfrog::new(h)),
                    "Trace" => Box::new(Trace::new(h)),
                    _ => return,
                };
                self.synchronize();
                self.integrator = integrator;
                self.dense_valid = false;
            }
            Command::Timestep(h) => {
                self.synchronize();
                self.integrator.set_timestep(sign * h.abs());
                self.dense_valid = false;
            }
            Command::End(end) => {
                self.end_override = Some(end.unwrap_or(sign * f64::INFINITY));
                if let Some(monitor) = &self.progress_monitor {
                    monitor.set_end(end);
                }
            }
            Command::Stop => {
                self.synchronize();
                self.end_override = Some(self.t);
                if let Some(monitor) = &self.progress_monitor {
                    monitor.note_stop();
                }
            }
            Command::SolarGR(on) => self.set_solar_gr(on),
            Command::Perturber(name, on) => {
                let _ = self.set_perturber(&name, on);
            }
        }
    }

    /// Whether the forces include solar general relativity ([`SolarGR`]).
    pub fn solar_gr(&self) -> bool {
        self.forces.iter().any(|f| f.name() == std::any::type_name::<SolarGR>())
    }

    /// Add solar general relativity ([`SolarGR`]) to the forces, or take it out.
    pub fn set_solar_gr(&mut self, on: bool) {
        if on == self.solar_gr() {
            return;
        }
        self.synchronize();
        self.dense_valid = false;
        if on {
            self.forces.push(Box::new(SolarGR));
        } else {
            self.forces.retain(|f| f.name() != std::any::type_name::<SolarGR>());
        }
    }

    /// Every particle with a mass (from its properties), except the central body, and whether
    /// its gravity is on ([`Simulation::set_perturber`]).
    pub fn perturbers(&self) -> Vec<(String, bool)> {
        let real: Vec<f64> = (0..self.len()).map(|i| self.real_mass(i)).collect();
        let central = central_body(&real);
        (0..self.len())
            .filter(|&i| real[i] != 0.0 && Some(i) != central)
            .map(|i| (self.names[i].clone(), self.masses[i] != 0.0))
            .collect()
    }

    /// Turn the gravity of the particle called `name` off (it then moves as a test particle)
    /// or back on. Its mass is kept in its properties meanwhile. The central body (the most
    /// massive) can't be turned off.
    pub fn set_perturber(&mut self, name: &str, on: bool) -> Result<(), SimulationError> {
        let idx = *self.particle_index_map.get(name).ok_or_else(|| SimulationError::ParticleNotFound(name.to_string()))?;
        let real: Vec<f64> = (0..self.len()).map(|i| self.real_mass(i)).collect();
        if real[idx] == 0.0 || central_body(&real) == Some(idx) {
            return Ok(());
        }
        let mass = if on { real[idx] } else { 0.0 };
        if self.masses[idx] == mass {
            return Ok(());
        }
        self.synchronize();
        self.dense_valid = false;
        self.masses[idx] = mass;
        if let Some(monitor) = &self.progress_monitor {
            monitor.reset_energy();
        }
        Ok(())
    }

    /// The mass in the `idx`-th particle's properties, which [`Simulation::set_perturber`]
    /// leaves alone.
    fn real_mass(&self, idx: usize) -> f64 {
        self.properties[idx].as_ref().and_then(|p| p.mass).unwrap_or(self.masses[idx])
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
        let target = self.days_since_ref(epoch);
        self.end_override = None;
        if let Some(monitor) = self.progress_monitor.clone() {
            monitor.begin_run(self, Some(target), true);
            self.integrate_or_interpolate_to(target);
            monitor.end_run(self);
        } else {
            self.integrate_or_interpolate_to(target);
        }
        self.end_override = None;
    }

    fn integrate_or_interpolate_to(&mut self, mut target: f64) {
        if !self.integrator.has_dense_output() {
            self.integrate_to(target);
            return;
        }
        self.synchronize();

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
            if matches!(self.checkpoint(), Go::Stop) {
                return;
            }
            if let Some(end) = self.end_override.take() {
                // The page moved the end: go there instead.
                self.synchronize();
                target = end;
                if (target - self.t).abs() < 1e-16 {
                    return;
                }
                continue;
            }
            if !self.integrator.has_dense_output() {
                // The page switched to an integrator without dense output.
                self.integrate_to(target);
                return;
            }
            let timestep = self.integrator.timestep();
            if (dt < 0.0) != (timestep < 0.0) {
                self.integrator.set_timestep(-timestep);
            }
            self.step_unmonitored();
            if (target - self.t) * dt < 0.0 && self.integrator.interpolate(target).is_none() {
                // Stepped past the target without being able to interpolate back; finish exactly.
                self.integrate_to(target);
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
        let target = self.days_since_ref(epoch);
        self.run_to(Some(target));
    }

    /// Integrate with no end: forwards (or backwards, with a negative timestep) until the
    /// progress page's Stop, an end set there, or the interrupt check
    /// ([`Simulation::set_interrupt`]) ends it. Without a progress page or an interrupt check it
    /// never returns.
    pub fn run(&mut self) {
        self.run_to(None);
    }

    /// Integrate to `target` (days since the reference epoch), or with no end.
    fn run_to(&mut self, target: Option<f64>) {
        self.end_override = None;
        let sign = if self.integrator.timestep() < 0.0 { -1.0 } else { 1.0 };
        let end = target.unwrap_or(sign * f64::INFINITY);
        if let Some(monitor) = self.progress_monitor.clone() {
            monitor.begin_run(self, target, true);
            self.integrate_to(end);
            monitor.end_run(self);
        } else {
            self.integrate_to(end);
        }
        self.end_override = None;
    }

    fn integrate_to(&mut self, mut target: f64) {
        self.synchronize();

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
            let single = match self.checkpoint() {
                Go::Stop => return,
                Go::OneStep => true,
                Go::Steps => false,
            };
            if let Some(end) = self.end_override.take() {
                // The page moved the end.
                self.synchronize();
                target = end;
            }
            let dt = target - self.t;
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
                self.step_unmonitored();
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
            // With a progress page or an interrupt check, at most a chunk of them before the next
            // checkpoint (so counting them stops there too).
            let cap = match (self.has_checkpoints(), single) {
                (false, _) if target.is_infinite() => MONITOR_CHUNK,
                (false, _) => usize::MAX,
                (true, true) => 1,
                (true, false) => MONITOR_CHUNK,
            };
            let n = if self.integrator.fixed_timestep() { self.full_steps_to(target, cap) } else { 1 };
            self.steps_unmonitored(n.max(1));
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
    /// Counts at most `cap`.
    fn full_steps_to(&self, target: f64, cap: usize) -> usize {
        let h = self.integrator.timestep();
        let mut t = self.t;
        let mut n = 0;
        loop {
            let dt = target - t;
            if n == cap || dt.abs() < 1e-16 || dt.abs() < h.abs() || (dt < 0.0) != (h < 0.0) {
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

/// Most steps [`Simulation::steps`] hands the integrator at once while a progress page or an
/// interrupt check is attached.
const MONITOR_CHUNK: usize = 1000;

/// Least time between two calls of the interrupt check.
const INTERRUPT_INTERVAL: std::time::Duration = std::time::Duration::from_millis(50);
