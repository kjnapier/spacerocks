use crate::constants::GRAVITATIONAL_CONSTANT;
use crate::nbody::integrators::{Integrator, Leapfrog};
use crate::nbody::forces::{central_body, Force};
use crate::state::{from_pv, position, pv, velocity, State};
use crate::transforms::universal_kepler_step;

use nalgebra::Vector3;
use rayon::prelude::*;

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
    /// Load inertial `states` with `masses`, reusing this state's buffers.
    pub(crate) fn load(&mut self, states: &[State], masses: &[f64], central: usize) {
        self.m.clear();
        self.m.extend_from_slice(masses);
        self.m_total = self.m.iter().sum();
        let mut x_cm = Vector3::zeros();
        let mut v_cm = Vector3::zeros();
        for (p, &mi) in states.iter().zip(&self.m) {
            if mi != 0.0 {
                x_cm += mi * position(p);
                v_cm += mi * velocity(p);
            }
        }
        x_cm /= self.m_total;
        v_cm /= self.m_total;
        let xc = position(&states[central]);
        self.q.clear();
        self.q.extend(states.iter().map(|p| position(p) - xc));
        self.u.clear();
        self.u.extend(states.iter().map(|p| velocity(p) - v_cm));
        self.massive.clear();
        self.massive.extend((0..states.len()).filter(|&i| i != central && self.m[i] != 0.0));
        self.central = central;
        self.m_central = self.m[central];
        self.x_cm = x_cm;
        self.v_cm = v_cm;
    }

    /// Write inertial states back into `states`.
    pub(crate) fn to_inertial(&self, states: &mut [State]) {
        let mut mq = Vector3::zeros();
        let mut mu = Vector3::zeros();
        for &i in &self.massive {
            mq += self.m[i] * self.q[i];
            mu += self.m[i] * self.u[i];
        }
        let xc = self.x_cm - mq / self.m_total;
        let vc = self.v_cm - mu / self.m_central;
        for (i, p) in states.iter_mut().enumerate() {
            if i == self.central {
                *p = from_pv(&xc, &vc);
            } else {
                *p = from_pv(&(self.q[i] + xc), &(self.u[i] + self.v_cm));
            }
        }
    }

    /// Kick the barycentric velocities by the forces minus the central body's Keplerian pull.
    /// The Newtonian pull between each pair in `skip` is left out too (TRACE moves those pairs
    /// into the drift), which assumes the forces include Newtonian gravity. Forces other than
    /// Newtonian gravity see the inertial state, written into `states` first.
    pub(crate) fn kick(&mut self, states: &mut [State], forces: &[Box<dyn Force + Send + Sync>], h: f64, skip: &[(usize, usize)]) {
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
                    self.to_inertial(states);
                    synced = true;
                }
                force.add_acceleration(states, &self.m, &mut acc);
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
        // Between massive bodies, each pair once, from the first of the two.
        for (k, &i) in self.massive.iter().enumerate() {
            let (qi, mi) = (self.q[i], self.m[i]);
            let mut acc_i = Vector3::zeros();
            for &j in &self.massive[k + 1..] {
                let g = pull(&self.q[j], &qi);
                acc_i += self.m[j] * g;
                acc[j] -= mi * g;
            }
            acc[i] += acc_i;
        }
        // On each test particle, from every massive body.
        for (j, (aj, qj)) in acc.iter_mut().zip(&self.q).enumerate() {
            if j == self.central || self.m[j] != 0.0 {
                continue;
            }
            for &i in &self.massive {
                *aj -= self.m[i] * pull(qj, &self.q[i]);
            }
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

    /// How far a jump of `h` moves every body: the massive bodies' momentum over the central
    /// body's mass, times `h`.
    fn jump_shift(&self, h: f64) -> Vector3<f64> {
        let mut p = Vector3::zeros();
        for &i in &self.massive {
            p += self.m[i] * self.u[i];
        }
        h * p / self.m_central
    }

    pub(crate) fn jump(&mut self, h: f64) {
        let dq = self.jump_shift(h);
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

/// Most steps per block in [`WisdomHolman::steps`] and [`Trace`](super::Trace)'s: the record
/// of the massive bodies for a block takes a few hundred bytes per massive body per step.
pub(crate) const BLOCK: usize = 512;

/// Test-particle steps below which blocks aren't worth it ([`WisdomHolman::steps`] and
/// [`Trace`](super::Trace)'s take their steps one at a time) or a block's test particles aren't
/// worth splitting between threads (handing work to the thread pool and waiting for it costs
/// tens of microseconds; a test-particle step takes about 0.2).
pub(crate) const MIN_WORK: usize = 4096;

/// Run `f` on every element of `items`, split between threads if `work` is at least
/// [`MIN_WORK`].
pub(crate) fn for_each_maybe_par<T: Send, F: Fn(&mut T) + Sync + Send>(items: &mut [T], work: usize, f: F) {
    if work >= MIN_WORK {
        items.par_iter_mut().for_each(f);
    } else {
        items.iter_mut().for_each(f);
    }
}

/// What a test particle needs from the massive bodies for one step (see [`Democratic`]: `q`
/// heliocentric positions, `u` barycentric velocities), besides their positions at the kicks.
#[derive(Clone, Copy, Default)]
pub(crate) struct Frame {
    /// Central body's position and the barycenter's velocity when the step loads the state.
    xc_load: Vector3<f64>,
    v_cm_load: Vector3<f64>,
    /// The two jumps.
    dq1: Vector3<f64>,
    dq2: Vector3<f64>,
    /// Central body's position and the barycenter's velocity when the step writes it back.
    xc_out: Vector3<f64>,
    v_cm_out: Vector3<f64>,
}

/// The massive bodies of a simulation stepped on their own, for integrators that then take the
/// test particles through the same steps separately. Test particles act on nothing, so a step
/// of the massive bodies alone does exactly the arithmetic it does with them present.
pub(crate) struct MassiveRun {
    /// Index of each massive body in the simulation's particles.
    pub(crate) idx: Vec<usize>,
    /// Inertial states and masses of the massive bodies.
    pub(crate) bodies: Vec<State>,
    pub(crate) masses: Vec<f64>,
    pub(crate) central: usize,
    pub(crate) s: Democratic,
    /// Masses of the massive bodies other than the central one, in order.
    pub(crate) pullers: Vec<f64>,
    pub(crate) gm: f64,
    /// Per recorded step: the frame, and the positions of the bodies in `pullers` at the first
    /// and at the second kick.
    pub(crate) frames: Vec<Frame>,
    pub(crate) kicks: Vec<Vector3<f64>>,
}

impl MassiveRun {
    pub(crate) fn new(states: &[State], masses: &[f64], central: usize) -> MassiveRun {
        let idx: Vec<usize> = (0..states.len()).filter(|&i| masses[i] != 0.0).collect();
        let central_m = idx.iter().position(|&i| i == central).unwrap();
        let bodies: Vec<State> = idx.iter().map(|&i| states[i]).collect();
        let masses: Vec<f64> = idx.iter().map(|&i| masses[i]).collect();
        let gm = GRAVITATIONAL_CONSTANT * masses[central_m];
        let mut s = Democratic::default();
        s.load(&bodies, &masses, central_m);
        let pullers = s.massive.iter().map(|&i| s.m[i]).collect();
        MassiveRun { idx, bodies, masses, central: central_m, s, pullers, gm, frames: Vec::new(), kicks: Vec::new() }
    }

    /// Forget the recorded steps.
    pub(crate) fn clear(&mut self) {
        self.frames.clear();
        self.kicks.clear();
    }

    /// One Wisdom–Holman step of the massive bodies (Newtonian gravity only), recorded.
    /// `check` sees the loaded state at the start of the step (`false`) and the state after
    /// the second kick (`true`); if it returns false the step is abandoned, nothing is recorded
    /// and the bodies keep their states.
    pub(crate) fn step(&mut self, forces: &[Box<dyn Force + Send + Sync>], h: f64, mut check: impl FnMut(&Democratic, bool) -> bool) -> bool {
        let half = 0.5 * h;
        let s = &mut self.s;
        s.load(&self.bodies, &self.masses, self.central);
        if !check(s, false) {
            return false;
        }
        let mut f = Frame { xc_load: position(&self.bodies[self.central]), v_cm_load: s.v_cm, ..Frame::default() };
        let mark = self.kicks.len();
        self.kicks.extend(s.massive.iter().map(|&i| s.q[i]));
        s.kick(&mut self.bodies, forces, half, &[]);
        f.dq1 = s.jump_shift(half);
        s.jump(half);
        s.kepler(h);
        f.dq2 = s.jump_shift(half);
        s.jump(half);
        self.kicks.extend(s.massive.iter().map(|&i| s.q[i]));
        s.kick(&mut self.bodies, forces, half, &[]);
        if !check(s, true) {
            self.kicks.truncate(mark);
            return false;
        }
        s.to_inertial(&mut self.bodies);
        f.xc_out = position(&self.bodies[self.central]);
        f.v_cm_out = s.v_cm;
        self.frames.push(f);
        true
    }

    /// A test particle's democratic heliocentric state (`q`, `u`) as recorded step `k` loads it
    /// from inertial `x`, `v`.
    #[inline]
    pub(crate) fn test_load(&self, k: usize, x: &Vector3<f64>, v: &Vector3<f64>) -> (Vector3<f64>, Vector3<f64>) {
        let f = &self.frames[k];
        (x - f.xc_load, v - f.v_cm_load)
    }

    /// Take a test particle (inertial `x`, `v`) through recorded step `k`, with the arithmetic
    /// of a full step. Returns its democratic heliocentric state (`q`, `u`) as loaded at the
    /// start of the step and after the second kick.
    #[inline]
    pub(crate) fn test_step(&self, k: usize, h: f64, x: &mut Vector3<f64>, v: &mut Vector3<f64>) -> [(Vector3<f64>, Vector3<f64>); 2] {
        let half = 0.5 * h;
        let np = self.pullers.len();
        let f = &self.frames[k];
        let at = &self.kicks[2 * np * k..2 * np * (k + 1)];
        let kick = |q: &Vector3<f64>, u: &mut Vector3<f64>, at: &[Vector3<f64>]| {
            let mut acc = Vector3::zeros();
            for (qi, &mi) in at.iter().zip(&self.pullers) {
                acc -= mi * pull(q, qi);
            }
            *u += half * acc;
        };
        let (mut q, mut u) = self.test_load(k, x, v);
        let start = (q, u);
        kick(&q, &mut u, &at[..np]);
        q += f.dq1;
        (q, u) = kepler_drift(&q, &u, self.gm, h);
        q += f.dq2;
        kick(&q, &mut u, &at[np..]);
        *x = q + f.xc_out;
        *v = u + f.v_cm_out;
        [start, (q, u)]
    }

    /// Write the massive bodies' states back into `states`.
    pub(crate) fn write_back(&self, states: &mut [State]) {
        for (&i, body) in self.idx.iter().zip(&self.bodies) {
            states[i] = *body;
        }
    }
}

impl WisdomHolman {
    fn steps_blocked(&mut self, states: &mut [State], masses: &[f64], t: &mut f64, central: usize, forces: &[Box<dyn Force + Send + Sync>], n: usize) {
        let h = self.timestep;
        let mut run = MassiveRun::new(states, masses, central);
        let test_idx: Vec<usize> = (0..states.len()).filter(|&i| masses[i] == 0.0).collect();
        let mut tests: Vec<(Vector3<f64>, Vector3<f64>)> = test_idx.iter().map(|&i| pv(&states[i])).collect();
        let n_tests = tests.len();
        let mut done = 0;
        while done < n {
            let len = BLOCK.min(n - done);
            run.clear();
            for _ in 0..len {
                run.step(forces, h, |_, _| true);
                *t += h;
            }
            // Each test particle through the block on its own.
            for_each_maybe_par(&mut tests, n_tests * len, |(x, v)| {
                for k in 0..len {
                    run.test_step(k, h, x, v);
                }
            });
            done += len;
        }
        run.write_back(states);
        for (&i, (x, v)) in test_idx.iter().zip(tests) {
            states[i] = from_pv(&x, &v);
        }
    }
}

/// Newtonian pull per unit mass towards `from` felt at `at`, negated: `G (at - from) / r^3`.
#[inline]
fn pull(at: &Vector3<f64>, from: &Vector3<f64>) -> Vector3<f64> {
    let d = at - from;
    let r2 = d.norm_squared();
    (GRAVITATIONAL_CONSTANT / (r2 * r2.sqrt())) * d
}

/// Advance a two-body orbit about a body of gravitational parameter `mu` by `dt`.
pub(crate) fn kepler_drift(r0: &Vector3<f64>, v0: &Vector3<f64>, mu: f64, dt: f64) -> (Vector3<f64>, Vector3<f64>) {
    universal_kepler_step(r0, v0, mu, dt).expect("WisdomHolman: universal Kepler solve failed")
}

impl Integrator for WisdomHolman {
    fn step(&mut self, states: &mut [State], masses: &[f64], t: &mut f64, forces: &[Box<dyn Force + Send + Sync>]) {
        // The central body is the most massive particle (the first, on ties).
        let central = match central_body(masses) {
            Some(i) => i,
            // Nothing to orbit: fall back to a plain drift-kick-drift step.
            _ => return Leapfrog::new(self.timestep).step(states, masses, t, forces),
        };

        let h = self.timestep;
        let s = &mut self.state;
        s.load(states, masses, central);
        s.kick(states, forces, 0.5 * h, &[]);
        s.jump(0.5 * h);
        s.kepler(h);
        s.jump(0.5 * h);
        s.kick(states, forces, 0.5 * h, &[]);
        s.to_inertial(states);

        *t += h;
    }

    /// Takes `n` steps. With Newtonian gravity as the only force and some test particles, the
    /// massive bodies are stepped first, a block of steps at a time, recording what the test
    /// particles need from them; then each test particle is taken through the block on its own,
    /// the particles split between threads. Test particles don't act on anything, so this does
    /// exactly the arithmetic of `n` calls to [`Integrator::step`] and gives identical results,
    /// with one hand-off to the thread pool per block instead of several per step.
    fn steps(&mut self, states: &mut [State], masses: &[f64], t: &mut f64, forces: &[Box<dyn Force + Send + Sync>], n: usize) {
        let central = central_body(masses);
        let newtonian_only = forces.len() == 1 && forces[0].is_newtonian_gravity();
        let n_tests = masses.iter().filter(|&&m| m == 0.0).count();
        match central {
            Some(central) if newtonian_only && n_tests * n >= MIN_WORK => self.steps_blocked(states, masses, t, central, forces, n),
            _ => {
                for _ in 0..n {
                    self.step(states, masses, t, forces);
                }
            }
        }
    }

    fn timestep(&self) -> f64 {
        self.timestep
    }

    fn set_timestep(&mut self, timestep: f64) {
        self.timestep = timestep;
    }

    fn fixed_timestep(&self) -> bool {
        true
    }
}
