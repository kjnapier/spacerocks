//! Astrometric residuals of a trial orbit and their partial derivatives.
//!
//! Follows layup's `compute_residuals`: the trial state is integrated with the ASSIST force
//! model and its variational equations, once forward through the detections after the epoch
//! and once backward through those before it, reading states off the integrator's dense
//! output. Each detection is light-time corrected by four fixed-point iterations.
//!
//! Optical residuals are observed minus computed, projected on the tangent plane of the
//! observed direction (`compute_optical_residuals`). Rates follow `compute_streak_residuals`,
//! and radar delay and Doppler follow `compute_radar_residuals`: a two-leg light time with the
//! transmitter placed at the transmit time, the exact round-trip Doppler, and the Shapiro
//! delay. Where the transmitting antenna's Earth-fixed position is known, its state at the
//! transmit time is computed exactly rather than extrapolated (see [`Astrometry`]).

use nalgebra::Vector3;

use crate::assist::forces::NonGravitational;
use crate::assist::{AdaptiveMode, SpiceSimulation};
use crate::data::SPEED_OF_LIGHT;
use crate::observing::earth_fixed_state;
use crate::orbfit::{Astrometry, FitOptions, RowKind};
use crate::spice::SpiceKernel;
use crate::time::Time;
use crate::SpaceRock;

type BoxError = Box<dyn std::error::Error + Send + Sync>;

const NONGRAV_NAMES: [&str; 3] = ["A1", "A2", "A3"];

/// Residuals of a trial orbit at every detection, and optionally their partial derivatives.
///
/// Residual rows are packed: each detection contributes the rows [`Astrometry::rows`] lists
/// (RA, Dec, then RA rate, Dec rate, then delay, Doppler, as measured), detections in input
/// order. Residuals are observed minus computed: along RA (times cos Dec) and Dec in radians,
/// rates in radians/day, delay in days and Doppler in AU/day.
#[derive(Debug, Clone)]
pub struct Residuals {
    pub resid: Vec<f64>,
    /// Row-major `rows x npar` matrix of d(resid)/d(parameters): the state at the epoch (x, y,
    /// z, vx, vy, vz), then the fitted non-gravitational parameters. Empty unless requested.
    pub jacobian: Vec<f64>,
    pub npar: usize,
    /// Detection and measurement behind each row.
    pub row_detection: Vec<usize>,
    pub row_kind: Vec<RowKind>,
}

impl Residuals {
    /// Residuals per detection as `[ra, dec, ra_rate, dec_rate, delay, doppler]`, NaN where not
    /// measured.
    pub fn per_detection(&self, n: usize) -> Vec<[f64; 6]> {
        let mut out = vec![[f64::NAN; 6]; n];
        for ((&i, &k), &r) in self.row_detection.iter().zip(&self.row_kind).zip(&self.resid) {
            out[i][k as usize] = r;
        }
        out
    }
}

/// Indices of the detections after `epoch` (time order) and at or before it (reverse time
/// order), as in layup's `create_sequences`.
pub(crate) fn sequences(astrometry: &Astrometry, epoch: f64) -> (Vec<usize>, Vec<usize>) {
    let order = astrometry.time_order();
    let forward: Vec<usize> = order.iter().copied().filter(|&i| astrometry.epoch[i] > epoch).collect();
    let mut reverse: Vec<usize> = order.iter().copied().filter(|&i| astrometry.epoch[i] <= epoch).collect();
    reverse.reverse();
    (forward, reverse)
}

/// Residuals of the barycentric J2000 `state` (AU, AU/day) at TDB Julian date `epoch`, with
/// non-gravitational parameters `nongrav` (AU/day^2). With `partials`, also the Jacobian with
/// respect to the state and the non-gravitational parameters selected in `opts.fit_nongrav`.
pub fn residuals(
    astrometry: &Astrometry,
    epoch: f64,
    state: &[f64; 6],
    nongrav: &[f64; 3],
    kernel: &SpiceKernel,
    opts: &FitOptions,
    partials: bool,
) -> Result<Residuals, BoxError> {
    residuals_split(astrometry, epoch, state, nongrav, None, kernel, opts, partials)
}

/// Residuals with piecewise-constant non-gravitational parameters (layup's `per_arc`): the
/// detections before `epoch` (arc A) feel `nongrav`, those after it (arc B) `nongrav_arc2`. With
/// `partials`, the Jacobian has the state, then arc A's parameters selected in
/// `opts.fit_nongrav`, then arc B's (`6 + 2 * nactive` columns).
#[allow(clippy::too_many_arguments)]
pub fn residuals_per_arc(
    astrometry: &Astrometry,
    epoch: f64,
    state: &[f64; 6],
    nongrav: &[f64; 3],
    nongrav_arc2: &[f64; 3],
    kernel: &SpiceKernel,
    opts: &FitOptions,
    partials: bool,
) -> Result<Residuals, BoxError> {
    residuals_split(astrometry, epoch, state, nongrav, Some(nongrav_arc2), kernel, opts, partials)
}

#[allow(clippy::too_many_arguments)]
fn residuals_split(
    astrometry: &Astrometry,
    epoch: f64,
    state: &[f64; 6],
    nongrav: &[f64; 3],
    nongrav_arc2: Option<&[f64; 3]>,
    kernel: &SpiceKernel,
    opts: &FitOptions,
    partials: bool,
) -> Result<Residuals, BoxError> {
    let active: Vec<usize> = if partials { (0..3).filter(|&k| opts.fit_nongrav[k]).collect() } else { Vec::new() };
    let nact = active.len();
    // Columns per pass, and in the output (arc B's amplitudes get their own block).
    let npass = 6 + nact;
    let npar = if nongrav_arc2.is_some() { 6 + 2 * nact } else { npass };
    let n = astrometry.len();
    let (forward, reverse) = sequences(astrometry, epoch);

    let after = nongrav_arc2.unwrap_or(nongrav);
    let run = |seq: &[usize], ng: &[f64; 3]| pass(astrometry, seq, epoch, state, ng, &active, kernel, opts, partials);
    let (f, r) = if opts.parallel {
        rayon::join(|| run(&forward, after), || run(&reverse, nongrav))
    } else {
        (run(&forward, after), run(&reverse, nongrav))
    };

    // Row offsets of each detection in the packed output.
    let mut offset = Vec::with_capacity(n + 1);
    let mut out = Residuals { resid: Vec::new(), jacobian: Vec::new(), npar, row_detection: Vec::new(), row_kind: Vec::new() };
    offset.push(0);
    for i in 0..n {
        for k in astrometry.rows(i) {
            out.row_detection.push(i);
            out.row_kind.push(k);
        }
        offset.push(out.row_kind.len());
    }
    let nrows = out.row_kind.len();
    out.resid = vec![0.0; nrows];
    if partials {
        out.jacobian = vec![0.0; nrows * npar];
    }
    // Where a pass's amplitude columns go: forward (arc B) after arc A's block when split.
    let split = nongrav_arc2.is_some();
    for (seq, rows, amp_at) in [(&forward, f?, if split { 6 + nact } else { 6 }), (&reverse, r?, 6)] {
        let mut k = 0;
        for &i in seq.iter() {
            let (a, b) = (offset[i], offset[i + 1]);
            out.resid[a..b].copy_from_slice(&rows.resid[k..k + b - a]);
            if partials {
                for (row, src) in (a..b).zip(k..k + b - a) {
                    let (dst, from) = (&mut out.jacobian[row * npar..(row + 1) * npar], &rows.jacobian[src * npass..(src + 1) * npass]);
                    dst[..6].copy_from_slice(&from[..6]);
                    dst[amp_at..amp_at + nact].copy_from_slice(&from[6..]);
                }
            }
            k += b - a;
        }
    }
    Ok(out)
}

/// A simulation holding the trial orbit (and its variational particles) at `epoch`.
pub(crate) fn trial_simulation(
    epoch: f64,
    state: &[f64; 6],
    nongrav: &[f64; 3],
    active: &[usize],
    kernel: &SpiceKernel,
    opts: &FitOptions,
    partials: bool,
) -> Result<SpiceSimulation, BoxError> {
    let t0 = Time::new(epoch, "tdb", "jd")?;
    let mut sim = SpiceSimulation::horizons(&t0, kernel).map_err(|e| e.to_string())?;
    if let Some([alpha, nm, nn, nk, r0]) = opts.gofr {
        sim.forces[0] = Box::new(NonGravitational { alpha, r0, m: nm, n: nn, k: nk });
    }
    sim.integrator.set_epsilon(opts.epsilon);
    sim.integrator.set_adaptive_mode(AdaptiveMode::Prs23);
    let mut rock = SpaceRock::from_xyz("rock", state[0], state[1], state[2], state[3], state[4], state[5], t0, "J2000", "SSB")
        .map_err(|e| e.to_string())?;
    if nongrav.iter().any(|&a| a != 0.0) {
        rock.set_nongrav(nongrav[0], nongrav[1], nongrav[2]);
    }
    sim.add(rock).map_err(|e| e.to_string())?;
    if partials {
        sim.add_full_variation("rock").map_err(|e| e.to_string())?;
        for &k in active {
            sim.add_variation(NONGRAV_NAMES[k], "rock").map_err(|e| e.to_string())?;
        }
    }
    Ok(sim)
}

/// Move `sim` to the emission time of light received at `t` (days since the simulation's
/// reference epoch) by an observer at `r_obs`: four fixed-point iterations starting from no
/// delay, as layup's `integrate_light_time`.
pub(crate) fn integrate_light_time(sim: &mut SpiceSimulation, t: f64, r_obs: &Vector3<f64>, kernel: &SpiceKernel) -> Result<(), BoxError> {
    let mut lt = 0.0;
    for _ in 0..4 {
        sim.integrate_rel(t - lt, kernel).map_err(|e| e.to_string())?;
        lt = (sim.state.particles[0].position - r_obs).norm() / SPEED_OF_LIGHT;
    }
    Ok(())
}

/// GM of the Sun (AU^3/day^2), for the Shapiro delay.
const GM_SUN: f64 = 2.9591220828559115e-4;

#[allow(clippy::too_many_arguments)]
fn pass(
    astrometry: &Astrometry,
    seq: &[usize],
    epoch: f64,
    state: &[f64; 6],
    nongrav: &[f64; 3],
    active: &[usize],
    kernel: &SpiceKernel,
    opts: &FitOptions,
    partials: bool,
) -> Result<Residuals, BoxError> {
    let npar = 6 + active.len();
    let mut out = Residuals { resid: Vec::with_capacity(2 * seq.len()), jacobian: Vec::new(), npar, row_detection: Vec::new(), row_kind: Vec::new() };
    if seq.is_empty() {
        return Ok(out);
    }
    if partials {
        out.jacobian.reserve(2 * seq.len() * npar);
    }
    let mut sim = trial_simulation(epoch, state, nongrav, active, kernel, opts, partials)?;
    let jd_ref = sim.state.jd_ref;
    let c = SPEED_OF_LIGHT;

    for &i in seq {
        let r_obs = astrometry.observer(i);
        integrate_light_time(&mut sim, astrometry.epoch[i] - jd_ref, &r_obs, kernel)?;

        // Geometry shared by every observable: the unit line of sight, the distance, and the
        // variational position partials with their projection on the line of sight.
        let p = &sim.state.particles[0];
        let (ra, va) = (p.position, p.velocity);
        let d = ra - r_obs;
        let dist = d.norm();
        let rho = d / dist;
        let invd = 1.0 / dist;
        let mut dk = [Vector3::zeros(); 9];
        let mut dvk = [Vector3::zeros(); 9];
        let mut ddist = [0.0; 9];
        if partials {
            for (j, vp) in sim.state.variational_particles.iter().take(npar).enumerate() {
                dk[j] = vp.position;
                dvk[j] = vp.velocity;
                ddist[j] = rho.dot(&vp.position);
            }
        }
        let mut row = |resid: f64, partial: &dyn Fn(usize) -> f64| {
            out.resid.push(resid);
            if partials {
                for j in 0..npar {
                    out.jacobian.push(partial(j));
                }
            }
        };

        if astrometry.has_optical(i) {
            let (a, dv) = astrometry.tangent_basis(i);
            let drho = |j: usize| dk[j] * invd - rho * (invd * ddist[j]) - va * (ddist[j] * invd / c);
            row(-rho.dot(&a), &|j| -drho(j).dot(&a));
            row(-rho.dot(&dv), &|j| -drho(j).dot(&dv));

            if astrometry.has_rates(i) {
                // Great-circle rates with the (1 - d(light time)/dt) factor, as Sorcha and
                // layup's compute_streak_residuals; the partials use the geometric Jacobian.
                let vo = astrometry.observer_velocity(i);
                let vr = va - vo;
                let q = rho.dot(&vr);
                let dr = va * (1.0 - q / c) - vo;
                let w = (dr - rho * q) * invd;
                let invd2 = invd * invd;
                let jac = |x: &Vector3<f64>| {
                    let (px, sx) = (x.dot(&rho), x.dot(&vr));
                    (-(x * q + vr * px + rho * (sx - 3.0 * px * q)) * invd2, (x - rho * px) * invd)
                };
                for (x, obs) in [(a, astrometry.ra_rate[i]), (dv, astrometry.dec_rate[i])] {
                    let (dr_dr, dr_dv) = jac(&x);
                    row(obs - x.dot(&w), &|j| {
                        let drk = dk[j] - va * (ddist[j] / c);
                        -(dr_dr.dot(&drk) + dr_dv.dot(&dvk[j]))
                    });
                }
            }
        }

        if astrometry.has_delay(i) || astrometry.has_doppler(i) {
            // Two-leg radar model (layup's compute_radar_residuals). The sim is at the bounce
            // (down-leg emission) time; place the transmitter at the transmit time.
            let vo = astrometry.observer_velocity(i);
            let ao = astrometry.observer_acceleration(i);
            let vr = va - vo;
            let q = rho.dot(&vr);
            let ltdenom = 1.0 + rho.dot(&va) / c;
            let tau_d = dist / c;
            let t_receive = astrometry.epoch[i] - jd_ref;
            let (rtx, vtx, rho_u, rhu) = if let Some((rt, vt)) = astrometry.transmitter(i) {
                let u = ra - rt;
                (rt, vt, u.norm(), u)
            } else if let Some(site) = astrometry.transmitter_site(i) {
                // The antenna's exact state at the transmit time, iterated with the up-leg
                // light time.
                let (mut rt, mut vt, mut u, mut rho_u) = (r_obs, vo, d, dist);
                let mut tau_u = tau_d;
                for _ in 0..3 {
                    let s = earth_fixed_state(&site, jd_ref + (t_receive - tau_d - tau_u), kernel).map_err(|e| e.to_string())?;
                    rt = Vector3::new(s[0], s[1], s[2]);
                    vt = Vector3::new(s[3], s[4], s[5]);
                    u = ra - rt;
                    rho_u = u.norm();
                    tau_u = rho_u / c;
                }
                (rt, vt, rho_u, u)
            } else {
                // Taylor-extrapolate the receiving station back to the transmit time. (layup has
                // -a tau^2 / 2 here; the expansion of x(t - tau) is +a tau^2 / 2.)
                let (mut rt, mut u, mut rho_u) = (r_obs, d, dist);
                let mut tau_u = tau_d;
                for _ in 0..3 {
                    let tau = tau_d + tau_u;
                    rt = r_obs - vo * tau + ao * (0.5 * tau * tau);
                    u = ra - rt;
                    rho_u = u.norm();
                    tau_u = rho_u / c;
                }
                (rt, vo - ao * (tau_d + tau_u), rho_u, u)
            };
            let tau_u = rho_u / c;
            let rhu = rhu / rho_u;

            if astrometry.has_delay(i) {
                let mut model = tau_d + tau_u;
                // Shapiro delay on each leg, with the Sun at the integrator's current time.
                let mut sun = [[0.0; 6]];
                let t_sim = sim.state.particles_1[0].epoch;
                if kernel.barycentric_states_au_rel(&[10], jd_ref, t_sim, &mut sun).is_ok() {
                    let s = Vector3::new(sun[0][0], sun[0][1], sun[0][2]);
                    let (rb, rr, rt) = ((ra - s).norm(), (r_obs - s).norm(), (rtx - s).norm());
                    let up = (rt + rb + rho_u) / (rt + rb - rho_u);
                    let down = (rb + rr + dist) / (rb + rr - dist);
                    if up > 0.0 && down > 0.0 {
                        model += 2.0 * GM_SUN / (c * c * c) * (up.ln() + down.ln());
                    }
                }
                row(astrometry.delay[i] - model, &|j| -2.0 * (ddist[j] / ltdenom) / c);
            }
            if astrometry.has_doppler(i) {
                let dt_bounce = (c + rho.dot(&vo)) / (c + rho.dot(&va));
                let dt_transmit = dt_bounce * (c - rhu.dot(&va)) / (c - rhu.dot(&vtx));
                let model = c * (1.0 - dt_transmit);
                row(astrometry.doppler[i] - model, &|j| {
                    let range_partial = ddist[j] / ltdenom;
                    let drk = dk[j] - va * (range_partial / c);
                    let dq = rho.dot(&dvk[j]) + (vr.dot(&drk) - q * rho.dot(&drk)) * invd;
                    -2.0 * dq
                });
            }
        }
    }
    Ok(out)
}
