use crate::SpaceRock;
use crate::time::Time;
use crate::constants::GRAVITATIONAL_CONSTANT;
use crate::nbody::integrators::{Integrator, Leapfrog};
use crate::nbody::integrators::wisdom_holman::{central_body, kepler_drift, Democratic};
use crate::nbody::forces::Force;

use nalgebra::Vector3;

/// TRACE: a time-reversible hybrid integrator (Lu, Hernandez & Rein 2024).
///
/// Away from close encounters a step is exactly a [`WisdomHolman`](super::WisdomHolman) step in
/// democratic heliocentric coordinates. Two kinds of encounter switch parts of the step to an
/// adaptive Bulirsch–Stoer integration:
///
/// * **Close pairs.** When two bodies other than the central one come within `r_crit_hill`
///   Hill radii, now or within half a step along their relative motion, their
///   mutual pull is taken out of the kicks, and the bodies in such pairs are drifted together
///   under the central body's Kepler potential plus their mutual pull with Bulirsch–Stoer.
///   Everything else still drifts analytically.
/// * **Pericenter passages.** When the step is longer than `peri_crit_eta` times a body's
///   Keplerian timescale about the central body (Pham, Rein & Spiegel 2024), where the
///   Kepler–jump splitting loses accuracy, the whole step is integrated with Bulirsch–Stoer
///   using the simulation's forces, with no splitting.
///
/// The switching is decided reversibly: the criteria are evaluated at the start of the step,
/// the step is taken, and they are evaluated again at the end. If the end flags anything the
/// start did not, the step is redone from the start with both sets flagged. A step and its
/// reverse therefore make the same choice, so the integrator stays time-reversible.
///
/// The switching criteria, the composition of the step and the pericenter handling follow
/// REBOUND's TRACE (its default `FULL_BS` pericenter mode), and with no encounters a step is
/// REBOUND's TRACE step to round-off. Unlike REBOUND, every particle is checked for pericenter
/// passages, test particles included, as REBOUND does when `N_active` is left unset.
///
/// Close-pair handling assumes the simulation's forces include Newtonian gravity (the default).
#[derive(PartialEq, Debug, Clone, Copy)]
pub struct Trace {
    /// Current timestep in simulation time units
    pub timestep: f64,
    /// Pairs closer than this many Hill radii (of the larger body) are close encounters.
    pub r_crit_hill: f64,
    /// A step longer than this many Pham–Rein–Spiegel timescales of any body's Keplerian
    /// motion triggers a pericenter step.
    pub peri_crit_eta: f64,
    /// Relative tolerance of the Bulirsch–Stoer integrations.
    pub bs_epsilon: f64,
}

impl Trace {
    /// Creates a TRACE integrator with the specified timestep and the default criteria
    /// (REBOUND's: 3 Hill radii and `peri_crit_eta` 1; Bulirsch–Stoer tolerance 1e-12).
    ///
    /// # Arguments
    ///
    /// * `timestep` - Fixed timestep to use for integration
    pub fn new(timestep: f64) -> Trace {
        Trace { timestep, r_crit_hill: 3.0, peri_crit_eta: 1.0, bs_epsilon: 1e-12 }
    }
}

/// What the switching criteria flagged.
#[derive(PartialEq, Debug, Clone, Default)]
struct Flags {
    peri: bool,
    pairs: Vec<(usize, usize)>,
}

impl Flags {
    fn union(&self, other: &Flags) -> Flags {
        let mut pairs = self.pairs.clone();
        for p in &other.pairs {
            if !pairs.contains(p) {
                pairs.push(*p);
            }
        }
        pairs.sort();
        Flags { peri: self.peri || other.peri, pairs }
    }

    fn contains(&self, other: &Flags) -> bool {
        (self.peri || !other.peri) && other.pairs.iter().all(|p| self.pairs.contains(p))
    }
}

impl Trace {
    /// The switching criteria of REBOUND's TRACE (`reb_integrator_trace_switch_default` and
    /// `reb_integrator_trace_switch_peri_default`), evaluated in democratic heliocentric
    /// coordinates: heliocentric positions, barycentric velocities.
    fn flags(&self, particles: &[SpaceRock], central: usize, h: f64) -> Flags {
        let m0 = particles[central].mass();
        let gm0 = GRAVITATIONAL_CONSTANT * m0;
        let n = particles.len();
        let s = Democratic::from_particles(particles, central);
        let mut flags = Flags::default();

        // Pericenter: the step is longer than peri_crit_eta times the Pham, Rein & Spiegel
        // (2024) timescale of the Keplerian motion about the central body.
        for i in 0..n {
            if i != central && h * h > self.peri_crit_eta * self.peri_crit_eta * prs_timescale2(&s.q[i], &s.u[i], gm0) {
                flags.peri = true;
            }
        }

        // Pairs: closer than r_crit_hill Hill radii (of the larger body, with heliocentric
        // distance for the semimajor axis) now, or at the closest point of their straight-line
        // relative motion within half a step towards the approach.
        let hill6: Vec<f64> = (0..n)
            .map(|i| {
                let d2 = s.q[i].norm_squared();
                let mr = if i == central { 0.0 } else { particles[i].mass() / (3.0 * m0) };
                d2 * d2 * d2 * mr * mr
            })
            .collect();
        let rc2 = self.r_crit_hill * self.r_crit_hill;
        let h2 = 0.5 * h.abs();
        // Each pair with at least one massive member (test particles never meet each other).
        for i in (0..n).filter(|&i| hill6[i] > 0.0) {
            for j in (0..n).filter(|&j| j != i && j != central && (j > i || hill6[j] == 0.0)) {
                let dcrit6 = rc2 * rc2 * rc2 * hill6[i].max(hill6[j]);
                let dx = s.q[i] - s.q[j];
                let rp = dx.norm_squared();
                let close = if rp * rp * rp < dcrit6 {
                    true
                } else {
                    let dv = s.u[i] - s.u[j];
                    let v2 = dv.norm_squared();
                    let qv = dx.dot(&dv);
                    if qv == 0.0 {
                        false
                    } else {
                        let d = if qv < 0.0 { 1.0 } else { -1.0 };
                        let tmin = -d * qv / v2;
                        let dmin2 = if tmin < h2 { rp - qv * qv / v2 } else { rp + 2.0 * d * qv * h2 + v2 * h2 * h2 };
                        dmin2 * dmin2 * dmin2 < dcrit6
                    }
                };
                if close {
                    flags.pairs.push((i.min(j), i.max(j)));
                }
            }
        }
        flags.pairs.sort();
        flags
    }

    fn try_step(&self, particles: &mut Vec<SpaceRock>, forces: &Vec<Box<dyn Force + Send + Sync>>, central: usize, flags: &Flags) {
        let h = self.timestep;
        if flags.peri {
            bs_full(particles, forces, h, self.bs_epsilon);
            return;
        }
        let mut s = Democratic::from_particles(particles, central);
        s.kick(particles, forces, 0.5 * h, &flags.pairs);
        s.jump(particles, 0.5 * h);

        // Drift: bodies in close pairs together with Bulirsch–Stoer, the rest analytically.
        let mut in_encounter = vec![false; particles.len()];
        for &(i, j) in &flags.pairs {
            in_encounter[i] = true;
            in_encounter[j] = true;
        }
        // When every flagged pair has a test particle, the massive bodies feel nothing extra
        // and keep their analytic drift (REBOUND's `tponly_encounter`).
        let tp_only = flags.pairs.iter().all(|&(i, j)| particles[i].mass() == 0.0 || particles[j].mass() == 0.0);
        let analytic = |i: usize| !in_encounter[i] || (tp_only && particles[i].mass() > 0.0);
        let gm0 = GRAVITATIONAL_CONSTANT * s.m_central;
        let before: Vec<_> = s.q.iter().zip(&s.u).map(|(q, u)| (*q, *u)).collect();
        for i in 0..particles.len() {
            if i != central && analytic(i) {
                let (q, u) = kepler_drift(&s.q[i], &s.u[i], gm0, h);
                s.q[i] = q;
                s.u[i] = u;
            }
        }
        if !flags.pairs.is_empty() {
            let after: Vec<_> = s.q.iter().zip(&s.u).map(|(q, u)| (*q, *u)).collect();
            for i in 0..particles.len() {
                if in_encounter[i] {
                    (s.q[i], s.u[i]) = before[i];
                }
            }
            bs_encounter(&mut s, particles, &in_encounter, &flags.pairs, h, self.bs_epsilon);
            for i in 0..particles.len() {
                if in_encounter[i] && analytic(i) {
                    (s.q[i], s.u[i]) = after[i];
                }
            }
        }
        s.x_cm += h * s.v_cm;

        s.jump(particles, 0.5 * h);
        s.kick(particles, forces, 0.5 * h, &flags.pairs);
        s.to_particles(particles);
    }
}

impl Integrator for Trace {
    fn step(&mut self, particles: &mut Vec<SpaceRock>, epoch: &mut Time, forces: &Vec<Box<dyn Force + Send + Sync>>) {
        let central = match central_body(particles) {
            Some(i) => i,
            // Nothing to orbit: fall back to a plain drift-kick-drift step.
            None => return Leapfrog::new(self.timestep).step(particles, epoch, forces),
        };

        let start: Vec<_> = particles.iter().map(|p| (p.position, p.velocity)).collect();
        let mut flags = self.flags(particles, central, self.timestep);
        loop {
            self.try_step(particles, forces, central, &flags);
            let end = self.flags(particles, central, -self.timestep);
            if flags.contains(&end) {
                break;
            }
            flags = flags.union(&end);
            for (p, &(x, v)) in particles.iter_mut().zip(&start) {
                p.position = x;
                p.velocity = v;
            }
        }

        *epoch += self.timestep;
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

/// Squared Pham, Rein & Spiegel (2024) timescale (their eq. 16) of Keplerian motion about a
/// body of `gm`, from the second to fourth time derivatives of the position.
fn prs_timescale2(x: &Vector3<f64>, v: &Vector3<f64>, gm: f64) -> f64 {
    let d2 = x.norm_squared();
    let d = d2.sqrt();
    let a = -gm / (d2 * d) * x;
    let jerk = gm / (d2 * d2 * d) * (-d2 * v + 3.0 * x.dot(v) * x);
    let xv = x.dot(v);
    // Fourth derivative, component by component as REBOUND writes it.
    let snap = |k: usize| {
        let (o1, o2) = ((k + 1) % 3, (k + 2) % 3);
        let others2 = x[o1] * x[o1] + x[o2] * x[o2];
        let inner = -a[k] * others2 + 2.0 * x[k] * x[k] * a[k] + v[k] * (x[o1] * v[o1] + x[o2] * v[o2])
            + x[k] * (4.0 * v[k] * v[k] + 3.0 * (x[o1] * a[o1] + v[o1] * v[o1] + x[o2] * a[o2] + v[o2] * v[o2]));
        let j = -v[k] * others2 + 2.0 * x[k] * x[k] * v[k] + 3.0 * x[k] * (x[o1] * v[o1] + x[o2] * v[o2]);
        gm / (d2 * d2 * d2 * d) * (d2 * inner - 5.0 * xv * j)
    };
    let s = Vector3::new(snap(0), snap(1), snap(2));
    let (an, sn) = (a.norm(), s.norm());
    2.0 * an * an / (jerk.norm_squared() + an * sn)
}

/// Drift the bodies in close pairs under the central body's Kepler potential plus the pull
/// within each flagged pair (the Hamiltonian the kicks left out), in democratic heliocentric
/// coordinates.
fn bs_encounter(s: &mut Democratic, particles: &[SpaceRock], in_encounter: &[bool], pairs: &[(usize, usize)], h: f64, eps: f64) {
    let idx: Vec<usize> = (0..particles.len()).filter(|&i| in_encounter[i]).collect();
    let slot = |i: usize| idx.iter().position(|&k| k == i).unwrap();
    let local_pairs: Vec<(usize, usize, f64, f64)> = pairs.iter()
        .map(|&(i, j)| (slot(i), slot(j), particles[i].mass(), particles[j].mass()))
        .collect();
    let gm0 = GRAVITATIONAL_CONSTANT * s.m_central;

    let mut y = Vec::with_capacity(6 * idx.len());
    for &i in &idx {
        y.extend_from_slice(s.q[i].as_slice());
        y.extend_from_slice(s.u[i].as_slice());
    }
    let rhs = |y: &[f64], dy: &mut [f64]| {
        let q = |k: usize| Vector3::new(y[6 * k], y[6 * k + 1], y[6 * k + 2]);
        for k in 0..idx.len() {
            let qk = q(k);
            let r = qk.norm();
            let a = -gm0 * qk / (r * r * r);
            dy[6 * k..6 * k + 3].copy_from_slice(&y[6 * k + 3..6 * k + 6]);
            dy[6 * k + 3..6 * k + 6].copy_from_slice(a.as_slice());
        }
        for &(a, b, ma, mb) in &local_pairs {
            let d = q(b) - q(a);
            let r = d.norm();
            let g = GRAVITATIONAL_CONSTANT * d / (r * r * r);
            for c in 0..3 {
                dy[6 * a + 3 + c] += mb * g[c];
                dy[6 * b + 3 + c] -= ma * g[c];
            }
        }
    };
    bulirsch_stoer(&mut y, rhs, h, eps);
    for (k, &i) in idx.iter().enumerate() {
        s.q[i] = Vector3::new(y[6 * k], y[6 * k + 1], y[6 * k + 2]);
        s.u[i] = Vector3::new(y[6 * k + 3], y[6 * k + 4], y[6 * k + 5]);
    }
}

/// Integrate every particle in inertial coordinates under the simulation's forces.
fn bs_full(particles: &mut Vec<SpaceRock>, forces: &Vec<Box<dyn Force + Send + Sync>>, h: f64, eps: f64) {
    let n = particles.len();
    let mut y = Vec::with_capacity(6 * n);
    for p in particles.iter() {
        y.extend_from_slice(p.position.as_slice());
        y.extend_from_slice(p.velocity.as_slice());
    }
    let mut scratch = particles.clone();
    let rhs = |y: &[f64], dy: &mut [f64]| {
        for (k, p) in scratch.iter_mut().enumerate() {
            p.position = Vector3::new(y[6 * k], y[6 * k + 1], y[6 * k + 2]);
            p.velocity = Vector3::new(y[6 * k + 3], y[6 * k + 4], y[6 * k + 5]);
        }
        let mut acc = vec![Vector3::zeros(); n];
        for force in forces {
            for (a, da) in acc.iter_mut().zip(force.calculate_acceleration(&mut scratch)) {
                *a += da;
            }
        }
        for k in 0..n {
            dy[6 * k..6 * k + 3].copy_from_slice(&y[6 * k + 3..6 * k + 6]);
            dy[6 * k + 3..6 * k + 6].copy_from_slice(acc[k].as_slice());
        }
    };
    bulirsch_stoer(&mut y, rhs, h, eps);
    for (k, p) in particles.iter_mut().enumerate() {
        p.position = Vector3::new(y[6 * k], y[6 * k + 1], y[6 * k + 2]);
        p.velocity = Vector3::new(y[6 * k + 3], y[6 * k + 4], y[6 * k + 5]);
    }
}

/// Advance `y' = f(y)` by `span` with adaptive Bulirsch–Stoer (Gragg's modified midpoint and
/// polynomial extrapolation in the squared substep), to relative tolerance `eps`.
fn bulirsch_stoer<F: FnMut(&[f64], &mut [f64])>(y: &mut Vec<f64>, mut f: F, span: f64, eps: f64) {
    const SEQ: [usize; 8] = [2, 4, 6, 8, 10, 12, 14, 16];
    let dim = y.len();
    let mut done = 0.0;
    let mut big_h = span;
    let mut f0 = vec![0.0; dim];
    let mut z0 = vec![0.0; dim];
    let mut z1 = vec![0.0; dim];
    let mut fz = vec![0.0; dim];

    while (span - done).abs() > 1e-15 * span.abs() {
        if (big_h.abs()) > (span - done).abs() {
            big_h = span - done;
        }
        f(y, &mut f0);
        let mut table: Vec<Vec<Vec<f64>>> = Vec::with_capacity(SEQ.len());
        let mut accepted = None;
        for (k, &n) in SEQ.iter().enumerate() {
            // Modified midpoint with n substeps.
            let hs = big_h / n as f64;
            z0.copy_from_slice(y);
            for d in 0..dim {
                z1[d] = y[d] + hs * f0[d];
            }
            for _ in 1..n {
                f(&z1, &mut fz);
                for d in 0..dim {
                    let z2 = z0[d] + 2.0 * hs * fz[d];
                    z0[d] = z1[d];
                    z1[d] = z2;
                }
            }
            f(&z1, &mut fz);
            let est: Vec<f64> = (0..dim).map(|d| 0.5 * (z1[d] + z0[d] + hs * fz[d])).collect();

            // Neville extrapolation to zero substep, row k.
            let mut row = vec![est];
            for j in 1..=k {
                let ratio = (n as f64 / SEQ[k - j] as f64).powi(2) - 1.0;
                let prev_row = &table[k - 1];
                let cur = &row[j - 1];
                let next: Vec<f64> = (0..dim).map(|d| cur[d] + (cur[d] - prev_row[j - 1][d]) / ratio).collect();
                row.push(next);
            }
            if k > 0 {
                let best = &row[k];
                let prev = &row[k - 1];
                // Relative error of each body's position and velocity vectors.
                let err = (0..dim / 3)
                    .map(|g| {
                        let (b, p) = (&best[3 * g..3 * g + 3], &prev[3 * g..3 * g + 3]);
                        let diff = ((b[0] - p[0]).powi(2) + (b[1] - p[1]).powi(2) + (b[2] - p[2]).powi(2)).sqrt();
                        let size = (b[0] * b[0] + b[1] * b[1] + b[2] * b[2]).sqrt();
                        if size > 0.0 { diff / size } else { diff }
                    })
                    .fold(0.0, f64::max);
                if err < eps {
                    accepted = Some((row[k].clone(), k));
                    break;
                }
            }
            table.push(row);
        }
        match accepted {
            Some((y_new, k)) => {
                y.copy_from_slice(&y_new);
                done += big_h;
                if k <= 4 {
                    big_h *= 1.5;
                } else if k >= 7 {
                    big_h *= 0.7;
                }
            }
            None => {
                big_h *= 0.25;
                assert!(big_h.abs() > 1e-12 * span.abs(), "Trace: Bulirsch-Stoer step size underflow");
            }
        }
    }
}
