use crate::SpaceRock;
use crate::time::Time;
use crate::constants::GRAVITATIONAL_CONSTANT;
use crate::nbody::integrators::{Integrator, Leapfrog};
use crate::nbody::forces::Force;
use crate::transforms::universal_kepler_step;

use nalgebra::Vector3;

/// A second-order Wisdom–Holman symplectic integrator in democratic heliocentric coordinates
/// (Duncan, Levison & Lee 1998).
///
/// The Hamiltonian is split into a Keplerian part about the central body, an interaction
/// part, and a "jump" part (from the central body's momentum), and each step is the
/// symmetric composition
///
/// kick(h/2) jump(h/2) kepler(h) jump(h/2) kick(h/2)
///
/// Positions are heliocentric (relative to the central body) and velocities barycentric.
/// The central body is the most massive particle. The Kepler drift is solved exactly with
/// universal variables, so for a test particle orbiting the central body alone the step is
/// exact, and with perturbers the energy error stays bounded instead of growing.
///
/// The kick uses the simulation's forces: the interaction acceleration on each body is its
/// total acceleration minus the central body's Keplerian pull, so extra forces (GR, J2, ...)
/// enter as perturbations. Velocity-dependent forces are evaluated at the start of each half
/// kick, which keeps the step explicit but not strictly symplectic.
///
/// The timestep is fixed. Close encounters between bodies other than the central one are not
/// handled, so choose a timestep well below the shortest orbital period (about 1/20 of it).
#[derive(PartialEq, Debug, Clone)]
pub struct WisdomHolman {
    /// Current timestep in simulation time units
    pub timestep: f64,
    /// Buffers reused from step to step
    state: Democratic,
}

impl WisdomHolman {
    /// Creates a new Wisdom–Holman integrator with the specified timestep.
    ///
    /// # Arguments
    ///
    /// * `timestep` - Fixed timestep to use for integration
    pub fn new(timestep: f64) -> WisdomHolman {
        WisdomHolman { timestep, state: Democratic::default() }
    }
}

/// Democratic heliocentric state: heliocentric positions `q`, barycentric velocities `u`,
/// and the barycenter's position and velocity. The central body's own `q` and `u` are unused.
#[derive(PartialEq, Debug, Clone, Default)]
pub(crate) struct Democratic {
    pub(crate) central: usize,
    pub(crate) m_central: f64,
    pub(crate) m_total: f64,
    /// Masses, read once per step.
    pub(crate) m: Vec<f64>,
    pub(crate) q: Vec<Vector3<f64>>,
    pub(crate) u: Vec<Vector3<f64>>,
    pub(crate) x_cm: Vector3<f64>,
    pub(crate) v_cm: Vector3<f64>,
    /// The massive bodies other than the central one.
    pub(crate) massive: Vec<usize>,
    /// Scratch for the kick's accelerations.
    acc: Vec<Vector3<f64>>,
}

impl Democratic {
    /// Load the state of `particles`, reusing this state's buffers.
    pub(crate) fn load(&mut self, particles: &[SpaceRock], central: usize) {
        self.m.clear();
        self.m.extend(particles.iter().map(|p| p.mass()));
        self.m_total = self.m.iter().sum();
        let mut x_cm = Vector3::zeros();
        let mut v_cm = Vector3::zeros();
        for (p, &mi) in particles.iter().zip(&self.m) {
            if mi != 0.0 {
                x_cm += mi * p.position;
                v_cm += mi * p.velocity;
            }
        }
        x_cm /= self.m_total;
        v_cm /= self.m_total;
        let xc = particles[central].position;
        self.q.clear();
        self.q.extend(particles.iter().map(|p| p.position - xc));
        self.u.clear();
        self.u.extend(particles.iter().map(|p| p.velocity - v_cm));
        self.massive.clear();
        self.massive.extend((0..particles.len()).filter(|&i| i != central && self.m[i] != 0.0));
        self.central = central;
        self.m_central = self.m[central];
        self.x_cm = x_cm;
        self.v_cm = v_cm;
    }

    /// Write inertial positions and velocities back into `particles`.
    pub(crate) fn to_particles(&self, particles: &mut [SpaceRock]) {
        let mut mq = Vector3::zeros();
        let mut mu = Vector3::zeros();
        for &i in &self.massive {
            mq += self.m[i] * self.q[i];
            mu += self.m[i] * self.u[i];
        }
        let xc = self.x_cm - mq / self.m_total;
        let vc = self.v_cm - mu / self.m_central;
        for (i, p) in particles.iter_mut().enumerate() {
            if i == self.central {
                p.position = xc;
                p.velocity = vc;
            } else {
                p.position = self.q[i] + xc;
                p.velocity = self.u[i] + self.v_cm;
            }
        }
    }

    /// Kick the barycentric velocities by the forces minus the central body's Keplerian pull.
    /// The Newtonian pull between each pair in `skip` is left out too (TRACE moves those pairs
    /// into the drift), which assumes the forces include Newtonian gravity.
    pub(crate) fn kick(&mut self, particles: &mut Vec<SpaceRock>, forces: &Vec<Box<dyn Force + Send + Sync>>, h: f64, skip: &[(usize, usize)]) {
        let n = self.q.len();
        let mut acc = std::mem::take(&mut self.acc);
        acc.clear();
        acc.resize(n, Vector3::zeros());

        // Newtonian gravity is the central body's Keplerian pull, which the drift handles, plus
        // the interactions, computed here from the heliocentric positions. Other forces see the
        // inertial state.
        let mut n_gravity = 0;
        let mut synced = false;
        for force in forces {
            if force.is_newtonian_gravity() {
                self.add_interactions(&mut acc);
                n_gravity += 1;
            } else {
                if !synced {
                    self.to_particles(particles);
                    synced = true;
                }
                force.add_acceleration(particles, &mut acc);
            }
        }
        // Without exactly one gravity force, the Keplerian pull the drift adds is still taken
        // out of (or the extra ones added to) the kick.
        if n_gravity != 1 {
            let gm = GRAVITATIONAL_CONSTANT * self.m_central * (n_gravity as f64 - 1.0);
            for i in 0..n {
                if i != self.central {
                    let r2 = self.q[i].norm_squared();
                    acc[i] -= gm / (r2 * r2.sqrt()) * self.q[i];
                }
            }
        }

        for i in 0..n {
            if i != self.central {
                self.u[i] += h * acc[i];
            }
        }
        for &(i, j) in skip {
            let d = self.q[j] - self.q[i];
            let r = d.norm();
            let g = GRAVITATIONAL_CONSTANT * d / (r * r * r);
            self.u[i] -= h * self.m[j] * g;
            self.u[j] += h * self.m[i] * g;
        }
        self.acc = acc;
    }

    /// Add the Newtonian pull between every pair of bodies other than the central one (test
    /// particles don't pull on each other).
    fn add_interactions(&self, acc: &mut [Vector3<f64>]) {
        for &i in &self.massive {
            let (qi, mi) = (self.q[i], self.m[i]);
            let mut acc_i = Vector3::zeros();
            for (j, ((qj, &mj), aj)) in self.q.iter().zip(&self.m).zip(acc.iter_mut()).enumerate() {
                // Each pair of massive bodies once, from the first of the two.
                if j == i || j == self.central || (mj != 0.0 && j < i) {
                    continue;
                }
                let d = qj - qi;
                let r2 = d.norm_squared();
                let g = (GRAVITATIONAL_CONSTANT / (r2 * r2.sqrt())) * d;
                acc_i += mj * g;
                *aj -= mi * g;
            }
            acc[i] += acc_i;
        }
    }

    /// Copy `other` into this state, reusing this state's buffers.
    pub(crate) fn copy_from(&mut self, other: &Democratic) {
        self.central = other.central;
        self.m_central = other.m_central;
        self.m_total = other.m_total;
        self.m.clone_from(&other.m);
        self.q.clone_from(&other.q);
        self.u.clone_from(&other.u);
        self.x_cm = other.x_cm;
        self.v_cm = other.v_cm;
        self.massive.clone_from(&other.massive);
    }

    pub(crate) fn jump(&mut self, h: f64) {
        let mut p = Vector3::zeros();
        for &i in &self.massive {
            p += self.m[i] * self.u[i];
        }
        let dq = h * p / self.m_central;
        for i in 0..self.q.len() {
            if i != self.central {
                self.q[i] += dq;
            }
        }
    }

    pub(crate) fn kepler(&mut self, h: f64) {
        let gm = GRAVITATIONAL_CONSTANT * self.m_central;
        for i in 0..self.q.len() {
            if i != self.central {
                let (q, u) = kepler_drift(&self.q[i], &self.u[i], gm, h);
                self.q[i] = q;
                self.u[i] = u;
            }
        }
        self.x_cm += h * self.v_cm;
    }
}

/// Advance a two-body orbit about a body of gravitational parameter `mu` by `dt`.
pub(crate) fn kepler_drift(r0: &Vector3<f64>, v0: &Vector3<f64>, mu: f64, dt: f64) -> (Vector3<f64>, Vector3<f64>) {
    universal_kepler_step(r0, v0, mu, dt).expect("WisdomHolman: universal Kepler solve failed")
}

/// The most massive particle (the first, on ties), if any has mass.
pub(crate) fn central_body(particles: &[SpaceRock]) -> Option<usize> {
    let best = particles.iter().enumerate().fold(None, |best: Option<(usize, f64)>, (i, p)| match best {
        Some((_, m)) if m >= p.mass() => best,
        _ => Some((i, p.mass())),
    });
    best.filter(|&(_, m)| m > 0.0).map(|(i, _)| i)
}

impl Integrator for WisdomHolman {
    fn step(&mut self, particles: &mut Vec<SpaceRock>, epoch: &mut Time, forces: &Vec<Box<dyn Force + Send + Sync>>) {
        // The central body is the most massive particle (the first, on ties).
        let central = match central_body(particles) {
            Some(i) => i,
            // Nothing to orbit: fall back to a plain drift-kick-drift step.
            _ => return Leapfrog::new(self.timestep).step(particles, epoch, forces),
        };

        let h = self.timestep;
        let s = &mut self.state;
        s.load(particles, central);
        s.kick(particles, forces, 0.5 * h, &[]);
        s.jump(0.5 * h);
        s.kepler(h);
        s.jump(0.5 * h);
        s.kick(particles, forces, 0.5 * h, &[]);
        s.to_particles(particles);

        *epoch += h;
        for particle in particles.iter_mut() {
            particle.epoch = epoch.clone();
        }
    }

    fn timestep(&self) -> f64 {
        self.timestep
    }

    fn set_timestep(&mut self, timestep: f64) {
        self.timestep = timestep;
    }
}
