//! Set up and manage a collection of gravitationally interacting particles.
use std::collections::HashMap;
use std::sync::Arc;

use crate::assist::PerturberCache;

use crate::SpaceRock;
use crate::time::Time;
use crate::ReferencePlane;

use crate::assist::forces::{assist_default_forces, Force, NewtonianGravity};
use crate::assist::EphemerisConstants;
use crate::assist::integrators::{Integrator, IAS15};

use crate::spice::{SpiceBody, SpiceKernel};

use nalgebra::Vector3;

#[derive(Debug, Clone)]
pub struct VariationalParticle {
    pub name: String,
    pub initial_state: [f64; 6],
    pub epoch: f64,
    pub position: Vector3<f64>,
    pub velocity: Vector3<f64>,
    pub acceleration: Vector3<f64>,
    pub parent: usize,
    /// Variations of the parent's non-gravitational parameters (A1, A2, A3) carried by this
    /// variational particle (constant in time; usually 0 or 1).
    pub nongrav: [f64; 3],
}

#[derive(Debug, Clone)]
pub struct SimulationParticle {
    pub mass: f64,
    pub epoch: f64,
    pub position: Vector3<f64>,
    pub velocity: Vector3<f64>,
    pub acceleration: Vector3<f64>,
    /// Non-gravitational parameters (A1, A2, A3) of the Marsden model, in AU/day^2.
    pub nongrav: [f64; 3],
}

/// Partial derivatives of one particle's acceleration, kept in
/// [`SimulationState::partials`] only while variational particles exist.
#[derive(Debug, Clone, Copy, Default)]
pub struct Partials {
    /// Jacobian of the acceleration: rows 3..6 hold d(acceleration)/d(position) in columns 0..3
    /// and d(acceleration)/d(velocity) in columns 3..6 (rows 0..3 are unused).
    pub stm: [[f64; 6]; 6],
    /// d(acceleration)/d(A1, A2, A3).
    pub nongrav: [[f64; 3]; 3],
}

#[derive(Debug, Clone)]
pub struct SimulationSpiceParticle {
    pub mass: f64,
    pub epoch: f64,
    pub position: Vector3<f64>,
    pub velocity: Vector3<f64>,
    pub acceleration: Vector3<f64>,
}

#[derive(Debug, Clone)]
pub struct SimulationState {
    /// Current epoch (TDB Julian date) of `particles`.
    pub epoch: f64,
    /// Reference epoch (TDB Julian date, the simulation's starting epoch). The `epoch` fields of
    /// the integrator's particles (`particles_0/1`, `variational_particles_0/1`,
    /// `spice_particles`) are days since `jd_ref`, as in ASSIST: a Julian date near 2.46e6 only
    /// resolves ~40 µs, which would misplace a fast object by tens of centimetres per step.
    pub jd_ref: f64,
    
    pub particles_0: Vec<SimulationParticle>,
    pub particles_1: Vec<SimulationParticle>,
    /// Acceleration partials of `particles_1`, one per particle when variational particles
    /// exist and empty otherwise.
    pub partials: Vec<Partials>,

    pub variational_particles_0: Vec<VariationalParticle>,
    pub variational_particles_1: Vec<VariationalParticle>,

    pub spice_particles: Vec<SimulationSpiceParticle>,
    pub spice_bodies: Vec<SpiceBody>,

    pub particle_index_map: HashMap<String, usize>,
    pub spice_particle_index_map: HashMap<String, usize>,

    pub particles: Vec<SpaceRock>,
    pub variational_particles: Vec<VariationalParticle>,
    pub reference_plane: ReferencePlane,
    pub origin: SpiceBody,

    /// Optional precomputed perturber ephemeris (see [`PerturberCache`]); epochs it does not
    /// cover are read from the kernel.
    pub perturber_cache: Option<Arc<PerturberCache>>,
}

impl SimulationState {
    /// Barycentric states of the perturbers `ids` at `jd`, from the cache when possible.
    #[inline]
    /// `t` is in days since `jd_ref`.
    pub fn perturber_states(&self, kernel: &SpiceKernel, ids: &[i32], t: f64, out: &mut [[f64; 6]]) -> Result<(), crate::spice::SpiceError> {
        if let Some(c) = &self.perturber_cache {
            if c.states_rel(ids, self.jd_ref, t, out) {
                return Ok(());
            }
        }
        kernel.barycentric_states_au_rel(ids, self.jd_ref, t, out)
    }
}

/// Perturbers used by [`SpiceSimulation::horizons`] (Sun, planets, Moon, Pluto and 16 asteroids).
pub const HORIZONS_PERTURBERS: [&str; 27] = [
    "sun",
    "mercury barycenter",
    "venus barycenter",
    "earth",
    "moon",
    "mars barycenter",
    "jupiter barycenter",
    "saturn barycenter",
    "uranus barycenter",
    "neptune barycenter",
    "pluto barycenter",
    "2000001",
    "2000002",
    "2000003",
    "2000004",
    "2000007",
    "2000010",
    "2000015",
    "2000016",
    "2000031",
    "2000052",
    "2000065",
    "2000087",
    "2000088",
    "2000107",
    "2000511",
    "2000704",
];

pub struct SpiceSimulation {
    pub state: SimulationState,
    pub integrator: Box<dyn Integrator + Send + Sync>,
    pub forces: Vec<Box<dyn Force + Send + Sync>>,
}

impl SpiceSimulation {

    /// A simulation at `epoch` with the Sun, planets, Moon, Pluto and 16 massive asteroids
    /// ([`HORIZONS_PERTURBERS`]) read from `kernel`, and ASSIST's default force model
    /// ([`crate::assist::forces::assist_default_forces`]: Newtonian gravity, Einstein-Infeld-
    /// Hoffmann relativity from the Sun, Earth J2-J4, solar J2 and non-gravitational forces for
    /// particles that have A1/A2/A3), with constants read from the kernel's planetary
    /// ephemeris. Use [`SpiceSimulation::set_forces`] to change the model.
    pub fn horizons(epoch: &Time, kernel: &SpiceKernel) -> Result<SpiceSimulation, Box<dyn std::error::Error>> {

        let reference_plane = ReferencePlane::from_str("J2000")?;
        let origin = SpiceBody::from_name("SSB")?;

        let state = SimulationState {
            epoch: epoch.tdb().jd(),
            jd_ref: epoch.tdb().jd(),
            particles: Vec::new(),
            particles_0: Vec::new(),
            particles_1: Vec::new(),
            partials: Vec::new(),
            variational_particles: Vec::new(),
            variational_particles_0: Vec::new(),
            variational_particles_1: Vec::new(),
            spice_particles: Vec::new(),
            particle_index_map: HashMap::new(),
            spice_particle_index_map: HashMap::new(),
            spice_bodies: Vec::new(),
            reference_plane,
            origin,
            perturber_cache: None,
        };

        let mut sim = SpiceSimulation {
            state,
            integrator: Box::new(IAS15::new(0.001)),
            forces: assist_default_forces(&EphemerisConstants::from_kernel(kernel)),
        };

        for name in HORIZONS_PERTURBERS.iter() {
            // uppercase the name
            let name = name.to_uppercase();
            let particle = SpiceBody::from_name(&name)?;
            sim.add_spice_body(particle, kernel)?;
        }

        Ok(sim)
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
    // pub fn giants(epoch: &Time, kernel: &SpiceKernel) -> Result<SpiceSimulation, Box<dyn std::error::Error>> {

    //     let reference_plane = ReferencePlane::from_str("J2000")?;
    //     let origin = SpiceBody::from_name("SSB")?;

    //     let mut state = SimulationState {
    //         epoch: epoch.tdb().jd(),
    //         particles: Vec::new(),
    //         particles_0: Vec::new(),
    //         particles_1: Vec::new(),
    //         spice_particles: Vec::new(),
    //         particle_index_map: HashMap::new(),
    //         spice_particle_index_map: HashMap::new(),
    //         spice_bodies: Vec::new(),
    //         reference_plane,
    //         origin,
    //     };

    //     let mut sim = SpiceSimulation {
    //         state,
    //         integrator: Box::new(IAS15::new(0.001)),
    //         forces: vec![Box::new(NewtonianGravity)],
    //     };

    //     // add sun, jupiter barycenter, saturn barycenter, uranus barycenter, neptune barycenter.
    //     for name in ["sun", "jupiter barycenter", "saturn barycenter", "uranus barycenter", "neptune barycenter"].iter() {
    //         // uppercase the name
    //         let name = name.to_uppercase();
    //         let particle = SpiceBody::from_name(&name)?;
    //         sim.add_spice_body(particle, kernel);
    //     }
    //     Ok(sim)
    // }

    /// NAIF IDs of [`HORIZONS_PERTURBERS`], in order.
    pub fn horizons_body_ids() -> Result<Vec<i32>, Box<dyn std::error::Error>> {
        HORIZONS_PERTURBERS
            .iter()
            .map(|n| Ok(SpiceBody::from_name(&n.to_uppercase())?.code))
            .collect()
    }

    /// Use a precomputed perturber ephemeris (shared, e.g. between parallel simulations).
    pub fn set_perturber_cache(&mut self, cache: Arc<PerturberCache>) {
        self.state.perturber_cache = Some(cache);
    }

    pub fn add(&mut self, body: SpaceRock) -> Result<(), Box<dyn std::error::Error>> {
        // Check if the body is already in the simulation
        // if self.state.particle_index_map.contains_key(&body.name) {
        //     return Err(SimulationError::BodyAlreadyExists(body.name.clone()));
        // }

        // Create a new simulation particle
        let simulation_particle = SimulationParticle {
            mass: body.mass(),
            position: body.position,
            velocity: body.velocity,
            acceleration: Vector3::zeros(),
            epoch: body.epoch.tdb().jd() - self.state.jd_ref,
            nongrav: body.properties.as_ref().and_then(|p| p.nongrav).unwrap_or([0.0; 3]),
        };

        // Add the simulation particle to the simulation
        self.state.particles.push(body.clone());
        self.state.particles_0.push(simulation_particle.clone());
        self.state.particles_1.push(simulation_particle);
        self.state.particle_index_map.insert(body.name.clone(), self.state.particles.len() - 1);

        Ok(())
    }

    /// Add the six variational particles (x, y, z, vx, vy, vz) of `parent`, which together
    /// integrate its state transition matrix.
    pub fn add_full_variation(&mut self, parent: &str) -> Result<(), Box<dyn std::error::Error>> {
        for s in ["x", "y", "z", "vx", "vy", "vz"].iter() {
            self.add_variation(s, parent)?;
        }
        Ok(())
    }

    /// Add a variational particle for one initial-state component of `parent`
    /// ("x", "y", "z", "vx", "vy", "vz") or one of its non-gravitational parameters
    /// ("A1", "A2", "A3"). The latter gives the partial derivatives of the trajectory with
    /// respect to that parameter.
    pub fn add_variation(&mut self, dimension: &str, parent: &str) -> Result<(), Box<dyn std::error::Error>> {
        let mut state = [0.0; 6];
        let mut nongrav = [0.0; 3];
        match dimension {
            "x" => state[0] = 1.0,
            "y" => state[1] = 1.0,
            "z" => state[2] = 1.0,
            "vx" => state[3] = 1.0,
            "vy" => state[4] = 1.0,
            "vz" => state[5] = 1.0,
            "A1" | "a1" => nongrav[0] = 1.0,
            "A2" | "a2" => nongrav[1] = 1.0,
            "A3" | "a3" => nongrav[2] = 1.0,
            _ => return Err(format!("invalid variation '{}' (expected x, y, z, vx, vy, vz, A1, A2 or A3)", dimension).into()),
        }
        self.add_variation_with_nongrav(&format!("dd{}", dimension), state, nongrav, parent)
    }

    /// Add a variational particle with an arbitrary initial state variation.
    pub fn add_variation_from_state(&mut self, name: &str, state: [f64; 6], parent: &str) -> Result<(), Box<dyn std::error::Error>> {
        self.add_variation_with_nongrav(name, state, [0.0; 3], parent)
    }

    /// Add a variational particle with an initial state variation and a variation of the
    /// parent's non-gravitational parameters (A1, A2, A3).
    pub fn add_variation_with_nongrav(&mut self, name: &str, state: [f64; 6], nongrav: [f64; 3], parent: &str) -> Result<(), Box<dyn std::error::Error>> {
        let parent_idx = *self
            .state
            .particle_index_map
            .get(parent)
            .ok_or_else(|| format!("no particle named '{}' in the simulation", parent))?;
        let variational_particle = VariationalParticle {
            name: name.to_string(),
            initial_state: state,
            epoch: self.state.epoch - self.state.jd_ref,
            position: Vector3::new(state[0], state[1], state[2]),
            velocity: Vector3::new(state[3], state[4], state[5]),
            acceleration: Vector3::zeros(),
            parent: parent_idx,
            nongrav,
        };
        self.state.variational_particles_0.push(variational_particle.clone());
        self.state.variational_particles_1.push(variational_particle.clone());
        self.state.variational_particles.push(variational_particle);
        Ok(())
    }

    pub fn add_spice_body(&mut self, body: SpiceBody, kernel: &SpiceKernel) -> Result<(), Box<dyn std::error::Error>> {
        // Check if the body is already in the simulation
        // if self.state.spice_particle_index_map.contains_key(&body.name) {
        //     return Err(SimulationError::BodyAlreadyExists(body.name.clone()));
        // }

        // Get the state of the body at the current epoch
        let [x, y, z, vx, vy, vz] = kernel.state_au(body.code, self.state.origin.code, self.state.epoch)?;
        let position = Vector3::new(x, y, z);
        let velocity = Vector3::new(vx, vy, vz);

        // Create a new spice particle
        let spice_particle = SimulationSpiceParticle {
            mass: body.mass,
            position: position,
            velocity: velocity,
            acceleration: Vector3::zeros(),
            epoch: self.state.epoch - self.state.jd_ref,
        };

        // Add the spice particle to the simulation
        self.state.spice_bodies.push(body.clone());
        self.state.spice_particles.push(spice_particle);
        self.state.spice_particle_index_map.insert(body.name.clone(), self.state.spice_particles.len() - 1);

        Ok(())
    }

    /// Replace the force model (see [`crate::assist::forces`]).
    pub fn set_forces(&mut self, forces: Vec<Box<dyn Force + Send + Sync>>) {
        self.forces = forces;
    }

    /// Use Newtonian point-mass gravity only (the force model before ASSIST's corrections were
    /// added).
    pub fn newtonian_only(&mut self) {
        self.forces = vec![Box::new(NewtonianGravity)];
    }

    /// Add a force to the simulation.
    ///
    /// # Arguments
    ///
    /// * `force` - The force to add to the simulation.
    pub fn add_force(&mut self, force: Box<dyn Force + Send + Sync>) {
        self.forces.push(force);
    }


    pub fn step(&mut self, kernel: &SpiceKernel) -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
        self.integrator.step(&mut self.state, &self.forces, kernel)
    }


    pub fn integrate(&mut self, epoch: &Time, kernel: &SpiceKernel) -> Result<(), Box<dyn std::error::Error>> {
        self.integrate_jd(epoch.tdb().jd(), kernel)
    }

    /// Integrate to a TDB Julian date. The simulation's `particles` then hold the states at
    /// exactly `t` (interpolated from the last step), while the integrator keeps its own state,
    /// so integrating through a sorted sequence of epochs costs no more than integrating to the
    /// last one.
    pub fn integrate_jd(&mut self, t: f64, kernel: &SpiceKernel) -> Result<(), Box<dyn std::error::Error>> {
        // Early return if the epoch is the same as the current epoch
        if t == self.state.epoch {
            return Ok(());
        }
        self.integrate_rel(t - self.state.jd_ref, kernel)
    }

    /// Like [`SpiceSimulation::integrate_jd`], with the target given in days since
    /// `state.jd_ref`. Use this when the target is computed (e.g. an observation epoch minus a
    /// light time): an absolute Julian date only resolves ~40 µs.
    pub fn integrate_rel(&mut self, t: f64, kernel: &SpiceKernel) -> Result<(), Box<dyn std::error::Error>> {
        if t == self.state.epoch - self.state.jd_ref && self.state.particles_0[0].epoch == self.state.particles_1[0].epoch {
            return Ok(());
        }
        let t_abs = self.state.jd_ref + t;

        // Set the timestep to be in the correct direction
        if (t < self.state.particles_0[0].epoch) && (t < self.state.particles_1[0].epoch) {
            self.integrator.set_timestep(-1.0 * self.integrator.timestep().abs());
        } else {
            self.integrator.set_timestep(self.integrator.timestep().abs());
        }

        // Step until the last step brackets the target epoch, then interpolate.
        loop {
            let a = self.state.particles_0[0].epoch;
            let b = self.state.particles_1[0].epoch;
            if a != b && (t - a) * (t - b) <= 0.0 {
                break;
            }
            if let Err(e) = self.step(kernel) {
                // A full adaptive step can overshoot the target and need ephemeris data past the
                // end of the loaded kernels. Retry once with a step ending exactly on the target.
                let now = self.state.particles_1[0].epoch;
                let remaining = t - now;
                if remaining != 0.0 && remaining.abs() < self.integrator.timestep().abs() {
                    self.integrator.set_timestep(remaining);
                    self.step(kernel).map_err(|e2| format!("{} (retry ending at the target epoch also failed: {})", e, e2))?;
                } else {
                    return Err(e.to_string().into());
                }
            }
        }

        self.interpolate_simulation(t)?;
        self.state.epoch = t_abs;

        Ok(())
    }

    /// Interpolate the particles to `time` (days since `jd_ref`) within the last step.
    pub fn interpolate_simulation(&mut self, time: f64) -> Result<(), Box<dyn std::error::Error>> {
        let bs_vector = &self.integrator.bs_last();

        // let h = 1.0 - (self.state.epoch - time) / self.integrator.last_timestep(); 
        // let h = 1.0 - (self.state.particles_0[0].epoch - time) / self.integrator.last_timestep(); 

        let h = (time - self.state.particles_0[0].epoch) / self.integrator.last_timestep(); 

        let mut s = vec![0.0; 9];
        s[0] = self.integrator.last_timestep() * h;
        s[1] =       s[0] * s[0] / 2.0;
        s[2] =       s[1] * h / 3.0;
        s[3] =       s[2] * h / 2.0;
        s[4] = 3.0 * s[3] * h / 5.0;
        s[5] = 2.0 * s[4] * h / 3.0;
        s[6] = 5.0 * s[5] * h / 7.0;
        s[7] = 3.0 * s[6] * h / 4.0;
        s[8] = 7.0 * s[7] * h / 9.0;

        let mut u = vec![0.0; 8];
        u[0] = self.integrator.last_timestep() * h;
        u[1] =      u[0] * h / 2.;
        u[2] = 2. * u[1] * h / 3.;
        u[3] = 3. * u[2] * h / 4.;
        u[4] = 4. * u[3] * h / 5.;
        u[5] = 5. * u[4] * h / 6.;
        u[6] = 6. * u[5] * h / 7.;
        u[7] = 7. * u[6] * h / 8.;

        // let mut z = vec![0.0; 7];
        // z[0] = h;
        // z[1] = h * h;
        // z[2] = z[1] * h;
        // z[3] = z[2] * h;
        // z[4] = z[3] * h;
        // z[5] = z[4] * h;
        // z[6] = z[5] * h;


        let new_time = self.state.particles_0[0].epoch + s[0];
        let new_time_object = Time::new(self.state.jd_ref + new_time, "tdb", "jd")?;

        for idx in 0..self.state.particles_0.len() {
            let p = &mut self.state.particles_0[idx];
            let bs = &bs_vector[idx];

            let new_pos = p.position + (s[8] * bs.p6 + s[7] * bs.p5 + s[6] * bs.p4 + s[5] * bs.p3 + s[4] * bs.p2 + s[3] * bs.p1 + s[2] * bs.p0 + s[1] * p.acceleration + s[0] * p.velocity);
            self.state.particles[idx].position = new_pos;

            let new_vel = p.velocity + (u[7] * bs.p6 + u[6] * bs.p5 + u[5] * bs.p4 + u[4] * bs.p3 + u[3] * bs.p2 + u[2] * bs.p1 + u[1] * bs.p0 + u[0] * p.acceleration);
            self.state.particles[idx].velocity = new_vel;

            // let new_acc = p.acceleration + bs.p0 * z[0] + bs.p1 * z[1] + bs.p2 * z[2] + bs.p3 * z[3] + bs.p4 * z[4] + bs.p5 * z[5] + bs.p6 * z[6];
            // self.state.particles[idx].acceleration = new_acc;

            self.state.particles[idx].epoch = new_time_object.clone();
        }
        // now do the same for the variational particles
        for idx in 0..self.state.variational_particles_0.len() {
            let p = &mut self.state.variational_particles_0[idx];
            let bs = &bs_vector[idx + self.state.particles_0.len()];

            let new_pos = p.position + (s[8] * bs.p6 + s[7] * bs.p5 + s[6] * bs.p4 + s[5] * bs.p3 + s[4] * bs.p2 + s[3] * bs.p1 + s[2] * bs.p0 + s[1] * p.acceleration + s[0] * p.velocity);
            self.state.variational_particles[idx].position = new_pos;

            let new_vel = p.velocity + (u[7] * bs.p6 + u[6] * bs.p5 + u[5] * bs.p4 + u[4] * bs.p3 + u[3] * bs.p2 + u[2] * bs.p1 + u[1] * bs.p0 + u[0] * p.acceleration);
            self.state.variational_particles[idx].velocity = new_vel;

            // let new_acc = p.acceleration + bs.p0 * z[0] + bs.p1 * z[1] + bs.p2 * z[2] + bs.p3 * z[3] + bs.p4 * z[4] + bs.p5 * z[5] + bs.p6 * z[6];
            // self.state.variational_particles[idx].acceleration = new_acc;

            self.state.variational_particles[idx].epoch = new_time;
        }

        self.state.epoch = self.state.jd_ref + new_time;

        Ok(())
    }

}
