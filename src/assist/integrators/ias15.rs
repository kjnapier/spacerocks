use crate::spice::SpiceKernel;

use crate::assist::{Partials, SimulationParticle, SimulationState};
use crate::assist::spice_simulation::VariationalParticle;
use crate::assist::integrators::Integrator;
use crate::assist::forces::Force;

use nalgebra::Vector3;

// Gauss Radau spacings
const H: [f64; 8] = [0.0, 
                     0.056_262_560_536_922_15, 
                     0.180_240_691_736_892_36, 
                     0.352_624_717_113_169_6, 
                     0.547_153_626_330_555_4, 
                     0.734_210_177_215_410_5, 
                     0.885_320_946_839_095_8, 
                     0.977_520_613_561_287_5];
// Other constants
const RR: [f64; 28] = [0.056_262_560_536_922_15, 
                       0.180_240_691_736_892_36, 
                       0.123_978_131_199_970_21, 
                       0.352_624_717_113_169_6, 
                       0.296_362_156_576_247_5, 
                       0.172_384_025_376_277_28, 
                       0.547_153_626_330_555_4, 
                       0.490_891_065_793_633_23, 
                       0.366_912_934_593_663_03, 
                       0.194_528_909_217_385_75, 
                       0.734_210_177_215_410_5, 
                       0.677_947_616_678_488_4, 
                       0.553_969_485_478_518_2, 
                       0.381_585_460_102_240_87, 
                       0.187_056_550_884_855_15, 
                       0.885_320_946_839_095_8, 
                       0.829_058_386_302_173_7, 
                       0.705_080_255_102_203_4, 
                       0.532_696_229_725_926_1, 
                       0.338_167_320_508_540_37, 
                       0.151_110_769_623_685_25, 
                       0.977_520_613_561_287_5, 
                       0.921_258_053_024_365_4, 
                       0.797_279_921_824_395_1, 
                       0.624_895_896_448_117_8, 
                       0.430_366_987_230_732_13, 
                       0.243_310_436_345_876_96, 
                       0.092_199_666_722_191_74];

const C: [f64; 21] = [-0.056_262_560_536_922_15, 
                       0.010_140_802_830_063_63, 
                      -0.236_503_252_273_814_52, 
                      -0.003_575_897_729_251_617_6, 
                       0.093_537_695_259_462_07, 
                      -0.589_127_969_386_984_2, 
                       0.001_956_565_409_947_221, 
                      -0.054_755_386_889_068_69, 
                       0.415_881_200_082_306_83, 
                      -1.136_281_595_717_539_6, 
                      -0.001_436_530_236_370_891_5, 
                       0.042_158_527_721_268_706, 
                      -0.360_099_596_502_056_8, 
                       1.250_150_711_840_691, 
                      -1.870_491_772_932_95, 
                       0.001_271_790_309_026_867_8, 
                      -0.038_760_357_915_906_77, 
                       0.360_962_243_452_846, 
                      -1.466_884_208_400_427, 
                       2.906_136_259_308_429_4, 
                      -2.755_812_719_772_045_7];

const D: [f64; 21] = [0.056_262_560_536_922_15, 
                      0.003_165_475_718_170_829_3, 
                      0.236_503_252_273_814_52, 
                      0.000_178_097_769_221_743_38, 
                      0.045_792_985_506_027_92, 
                      0.589_127_969_386_984_2, 
                      0.000_010_020_236_522_329_128, 
                      0.008_431_857_153_525_702, 
                      0.253_534_069_054_569_27, 
                      1.136_281_595_717_539_6, 
                      0.000_000_563_764_163_931_820_8, 
                      0.001_529_784_002_500_465_7, 
                      0.097_834_236_532_444_01, 
                      0.875_254_664_684_091_1, 
                      1.870_491_772_932_95, 
                      0.000_000_031_718_815_401_761_364, 
                      0.000_276_293_090_982_647_7, 
                      0.036_028_553_983_736_46, 
                      0.576_733_000_277_078_7, 
                      2.248_588_760_769_16, 
                      2.755_812_719_772_045_7];

const SAFETY_FACTOR: f64 = 0.25;

#[derive(Debug, Clone)]
/// IAS15 (Implicit integrator with Adaptive Step size control, 15th order) numerical integrator.
/// 
/// This implements a high-precision integrator based on the Gauss-Radau quadrature.
/// It features:
/// - 15th order accuracy
/// - Adaptive timestep control
/// - Iterative refinement at each step using carefully chosen substep positions
/// 
/// The method uses a predictor-corrector scheme with Gauss-Radau spacings to achieve
/// high accuracy while maintaining reasonable performance.
pub struct IAS15 {
    /// Current timestep in simulation time units
    pub timestep: f64,
    /// Desired precision of the integrator
    pub epsilon: f64,
    /// Last timestep used by the integrator
    pub last_timestep: f64,
    /// Current step coefficients
    pub bs: Vec<CoefficientSeptet>,
    /// Intermediate step coefficients
    pub gs: Vec<CoefficientSeptet>,
    /// Error estimate coefficients
    pub es: Vec<CoefficientSeptet>,
    /// Previous step coefficients
    pub bs_last: Vec<CoefficientSeptet>,
    /// Previous error estimate coefficients
    pub es_last: Vec<CoefficientSeptet>,
    /// Smallest allowed |timestep| (days). The adaptive controller never goes below this; at
    /// this step size the predictor-corrector must still converge or `step` returns an error.
    pub min_timestep: f64,
    /// How the next timestep is chosen (see [`AdaptiveMode`]).
    pub adaptive_mode: AdaptiveMode,
    /// How round-off is controlled when summing (see [`Summation`]).
    pub summation: Summation,
    /// Compensation terms carried from step to step for each particle's position and velocity
    /// (the part of the value that did not fit in the f64).
    pub csx: Vec<Vector3<f64>>,
    pub csv: Vec<Vector3<f64>>,
}

/// Round-off control in IAS15's sums.
///
/// Positions and velocities are accumulated over many steps, so without compensation their
/// round-off grows with the number of steps. Measured on round trips of 64 Holman clones
/// (validation/assist), compensating them cuts the round-off-limited error 2-10 times
/// (e.g. 1.7 m to 0.15 m over 10^5 days at epsilon 1e-11) for 0-3% more time per step.
/// Also compensating the corrector's sums, as REBOUND does, cost about 8% more and was no
/// more accurate in those tests.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Summation {
    /// Plain floating-point sums.
    None,
    /// Kahan-compensated positions and velocities: the compensation terms are carried from step
    /// to step and fed back into the predictor (the default).
    Kahan,
    /// As `Kahan`, and the corrector's acceleration differences and series coefficients are
    /// compensated too (REBOUND's and ASSIST's scheme).
    Full,
}

impl Summation {
    pub fn from_str(s: &str) -> Result<Summation, String> {
        match s.to_lowercase().replace(['-', ' '], "_").as_str() {
            "none" | "plain" | "off" => Ok(Summation::None),
            "kahan" | "compensated" => Ok(Summation::Kahan),
            "full" | "rebound" => Ok(Summation::Full),
            _ => Err(format!("unknown summation '{}' (expected none, kahan or full)", s)),
        }
    }

    pub fn as_str(&self) -> &'static str {
        match self {
            Summation::None => "none",
            Summation::Kahan => "kahan",
            Summation::Full => "full",
        }
    }
}

/// Timestep criteria, as in REBOUND's `ri_ias15.adaptive_mode`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AdaptiveMode {
    /// Largest last-term coefficient over largest acceleration, over all particles
    /// (REBOUND's original criterion; the one ASSIST uses).
    Global,
    /// Largest per-component ratio of last-term coefficient to acceleration.
    Individual,
    /// Pham, Rein & Spiegel (2024): timescale from acceleration, jerk and snap
    /// (REBOUND's default since 2024).
    Prs23,
    /// The `Global` criterion evaluated for each particle on its own (largest last-term
    /// coefficient over largest acceleration of that particle), taking the most demanding
    /// particle. Identical to `Global` for a single particle; with many particles it cannot be
    /// masked by another particle's larger acceleration.
    PerParticle,
}

impl AdaptiveMode {
    pub fn from_str(s: &str) -> Result<AdaptiveMode, String> {
        match s.to_lowercase().replace(['-', ' '], "_").as_str() {
            "prs23" | "prs" => Ok(AdaptiveMode::Prs23),
            "global" | "assist" | "legacy" => Ok(AdaptiveMode::Global),
            "individual" => Ok(AdaptiveMode::Individual),
            "per_particle" | "perparticle" => Ok(AdaptiveMode::PerParticle),
            _ => Err(format!("unknown adaptive mode '{}' (expected prs23, global, individual or per_particle)", s)),
        }
    }

    pub fn as_str(&self) -> &'static str {
        match self {
            AdaptiveMode::Prs23 => "prs23",
            AdaptiveMode::Global => "global",
            AdaptiveMode::Individual => "individual",
            AdaptiveMode::PerParticle => "per_particle",
        }
    }
}

impl IAS15 {
    /// Creates a new IAS15 integrator with the specified initial timestep
    ///
    /// # Arguments
    ///
    /// * `timestep` - Initial integration timestep in simulation time units
    ///
    /// The integrator will automatically adjust this timestep based on the 
    /// local truncation error to maintain the specified precision (epsilon).
    pub fn new(timestep: f64) -> IAS15 {
        IAS15 { timestep, epsilon: 1e-9, last_timestep: 0.0, bs: vec![], gs: vec![], es: vec![], bs_last: vec![], es_last: vec![], min_timestep: 1e-8, adaptive_mode: AdaptiveMode::Prs23, summation: Summation::Kahan, csx: vec![], csv: vec![] }
    }

    /// Resets all coefficient vectors to zero for the specified number of particles
    ///
    /// # Arguments
    ///
    /// * `n` - Number of particles in the simulation
    ///
    /// This is called when the number of particles changes or at the start of integration
    /// to ensure proper sizing of all coefficient vectors.
    /// Discard the current attempt and retry the step with `new_timestep`: restore the
    /// particles to the start of the step and re-predict the b coefficients for the new size.
    fn retry_with_timestep(
        &mut self,
        new_timestep: f64,
        state: &mut SimulationState,
        initial_particles: &[SimulationParticle],
        initial_variational_particles: &[VariationalParticle],
    ) {
        self.timestep = new_timestep;
        state.particles_1.clone_from_slice(initial_particles);
        state.variational_particles_1.clone_from_slice(initial_variational_particles);
        if self.last_timestep != 0.0 {
            let ratio = self.timestep / self.last_timestep;
            predict_next_coefficients(&ratio, &self.es_last, &self.bs_last, &mut self.es, &mut self.bs);
        } else {
            let n = self.bs.len();
            self.bs = vec![CoefficientSeptet::zeros(); n];
            self.es = vec![CoefficientSeptet::zeros(); n];
        }
    }

    pub fn reset_coefficients(&mut self, n: usize) {
        self.csx = vec![Vector3::zeros(); n];
        self.csv = vec![Vector3::zeros(); n];
        self.bs = vec![CoefficientSeptet::zeros(); n];
        self.gs = vec![CoefficientSeptet::zeros(); n];
        self.es = vec![CoefficientSeptet::zeros(); n];
        self.bs_last = vec![CoefficientSeptet::zeros(); n];
        self.es_last = vec![CoefficientSeptet::zeros(); n];
    }

}

impl Integrator for IAS15 {
    fn step(&mut self, state: &mut SimulationState, forces: &Vec<Box<dyn Force + Send + Sync>>, kernel: &SpiceKernel) -> StepResult {
        let nreal = state.particles_0.len();
        let nvar = state.variational_particles_0.len();
        let n = nreal + nvar;
        let with_stm = nvar > 0;
        let sum = self.summation;
        let full = sum == Summation::Full;

        let spice_ids: Vec<i32> = state.spice_bodies.iter().map(|b| b.code).collect();
        let nb = spice_ids.len();
        let mut buf = vec![[0.0f64; 6]; nb];

        // Forces at the start of the step.
        let t_begin = state.particles_1[0].epoch;
        state.perturber_states(kernel, &spice_ids, t_begin, &mut buf)?;
        for (sp, s) in state.spice_particles.iter_mut().zip(&buf) {
            sp.position = Vector3::new(s[0], s[1], s[2]);
            sp.velocity = Vector3::new(s[3], s[4], s[5]);
            sp.epoch = t_begin;
        }
        update_forces(state, forces, with_stm);

        let initial_particles = state.particles_1.clone();
        let initial_variational_particles = state.variational_particles_1.clone();

        // Start-of-step state of every particle (real, then variational).
        let mut x0: Vec<Vector3<f64>> = Vec::with_capacity(n);
        let mut v0: Vec<Vector3<f64>> = Vec::with_capacity(n);
        let mut a0: Vec<Vector3<f64>> = Vec::with_capacity(n);
        for p in &initial_particles {
            x0.push(p.position);
            v0.push(p.velocity);
            a0.push(p.acceleration);
        }
        for p in &initial_variational_particles {
            x0.push(p.position);
            v0.push(p.velocity);
            a0.push(p.acceleration);
        }

        if self.bs.len() != n || self.gs.len() != n || self.csx.len() != n {
            self.reset_coefficients(n);
        }

        // After a change of direction the predicted coefficients (extrapolated from the last
        // step, in the other direction) are meaningless: start the predictor from zero.
        if self.last_timestep * self.timestep < 0.0 {
            for c in self.bs.iter_mut().chain(self.es.iter_mut()) {
                *c = CoefficientSeptet::zeros();
            }
        }
        if self.timestep.abs() < self.min_timestep {
            self.timestep = self.min_timestep.copysign(self.timestep);
        }

        let mut rri = [0.0f64; 28];
        for (r, x) in rri.iter_mut().zip(RR.iter()) {
            *r = 1.0 / x;
        }
        let mut at = vec![Vector3::zeros(); n];
        let mut csb = vec![CoefficientSeptet::zeros(); n];
        let mut substeps = vec![[0.0f64; 6]; 8 * nb];

        'integration_loop: loop {
            let dt = self.timestep;
            for c in csb.iter_mut() {
                *c = CoefficientSeptet::zeros();
            }
            for (g, b) in self.gs.iter_mut().zip(&self.bs) {
                g.p0 = b.p6 * D[15] + b.p5 * D[10] + b.p4 * D[6] + b.p3 * D[3] + b.p2 * D[1] + b.p1 * D[0] + b.p0;
                g.p1 = b.p6 * D[16] + b.p5 * D[11] + b.p4 * D[7] + b.p3 * D[4] + b.p2 * D[2] + b.p1;
                g.p2 = b.p6 * D[17] + b.p5 * D[12] + b.p4 * D[8] + b.p3 * D[5] + b.p2;
                g.p3 = b.p6 * D[18] + b.p5 * D[13] + b.p4 * D[9] + b.p3;
                g.p4 = b.p6 * D[19] + b.p5 * D[14] + b.p4;
                g.p5 = b.p6 * D[20] + b.p5;
                g.p6 = b.p6;
            }
            // Predictor tables for this step size. Position at substep s:
            //   x0 + v0 dt h + dt^2 h^2 (a0/2 + sum_j b_j h^(j+1) / ((j+2)(j+3))),
            // velocity: v0 + dt h (a0 + sum_j b_j h^(j+1) / (j+2)).
            // px[s] = [dt^2 h^(j+3)/((j+2)(j+3)) for j in 0..7, dt^2 h^2/2, dt h]; pv likewise.
            let mut px = [[0.0f64; 9]; 8];
            let mut pv = [[0.0f64; 8]; 8];
            for s in 1..8 {
                let h = H[s];
                let mut hp = h * h; // h^(j+2)
                for j in 0..7 {
                    let jf = j as f64;
                    pv[s][j] = hp / (jf + 2.0) * dt;
                    px[s][j] = hp * h / ((jf + 2.0) * (jf + 3.0)) * dt * dt;
                    hp *= h;
                }
                px[s][7] = h * h / 2.0 * dt * dt;
                px[s][8] = h * dt;
                pv[s][7] = h * dt;
            }
            // Perturber states at the seven substeps (they do not change between iterations).
            for s in 1..8 {
                state.perturber_states(kernel, &spice_ids, t_begin + dt * H[s], &mut substeps[s * nb..(s + 1) * nb])?;
            }

            let mut predictor_corrector_error = 1e300;
            let mut predictor_corrector_error_last = 2.0;
            let mut iterations = 0;
            loop {
                if predictor_corrector_error < 1e-16 {
                    break;
                }
                if iterations > 2 && predictor_corrector_error_last <= predictor_corrector_error {
                    break;
                }
                if iterations >= 12 {
                    if dt.abs() <= self.min_timestep {
                        return Err(format!(
                            "IAS15 predictor-corrector did not converge at the minimum timestep ({:e} days) at epoch {} (error {:e})",
                            self.min_timestep, state.jd_ref + t_begin, predictor_corrector_error
                        )
                        .into());
                    }
                    let new_timestep = (dt.abs() / 2.0).max(self.min_timestep).copysign(dt);
                    self.retry_with_timestep(new_timestep, state, &initial_particles, &initial_variational_particles);
                    continue 'integration_loop;
                }
                predictor_corrector_error_last = predictor_corrector_error;
                predictor_corrector_error = 0.0;
                iterations += 1;

                for s in 1..8 {
                    let h = H[s];
                    let ts = t_begin + dt * h;
                    // Predict positions and velocities at the substep: the series in h with its
                    // coefficients (which depend only on h and dt) taken from the tables.
                    let (cx, cv) = (&px[s], &pv[s]);
                    for k in 0..n {
                        let b = self.bs[k].as_array();
                        let mut dx = b[6] * cx[6];
                        let mut dv = b[6] * cv[6];
                        for j in (0..6).rev() {
                            dx += b[j] * cx[j];
                            dv += b[j] * cv[j];
                        }
                        let dx = dx + a0[k] * cx[7] + v0[k] * cx[8];
                        let dv = dv + a0[k] * cv[7];
                        let (x, v) = match sum {
                            Summation::None => (x0[k] + dx, v0[k] + dv),
                            _ => ((dx - self.csx[k]) + x0[k], (dv - self.csv[k]) + v0[k]),
                        };
                        if k < nreal {
                            let p = &mut state.particles_1[k];
                            p.position = x;
                            p.velocity = v;
                            p.epoch = ts;
                        } else {
                            let p = &mut state.variational_particles_1[k - nreal];
                            p.position = x;
                            p.velocity = v;
                        }
                    }
                    for (i, sp) in state.spice_particles.iter_mut().enumerate() {
                        let st = &substeps[s * nb + i];
                        sp.position = Vector3::new(st[0], st[1], st[2]);
                        sp.velocity = Vector3::new(st[3], st[4], st[5]);
                        sp.epoch = ts;
                    }
                    update_forces(state, forces, with_stm);
                    for k in 0..nreal {
                        at[k] = state.particles_1[k].acceleration;
                    }
                    for k in nreal..n {
                        at[k] = state.variational_particles_1[k - nreal].acceleration;
                    }

                    // Corrector: update g[s-1] and the b coefficients.
                    let ro = s * (s - 1) / 2; // offset into RR
                    let co = if s >= 2 { (s - 1) * (s - 2) / 2 } else { 0 }; // offset into C
                    let mut max_ak = 0.0f64;
                    let mut max_db6 = 0.0f64;
                    for k in 0..n {
                        // at - a0, compensated.
                        let gk = if full {
                            let mut gk = at[k];
                            let mut cs = Vector3::zeros();
                            add_comp(true, &mut gk, &mut cs, -a0[k]);
                            gk - cs
                        } else {
                            at[k] - a0[k]
                        };
                        let g = self.gs[k].as_array_mut();
                        let mut val = gk * rri[ro];
                        for j in 1..s {
                            val = (val - g[j - 1]) * rri[ro + j];
                        }
                        let tmp = val - g[s - 1];
                        g[s - 1] = val;
                        let b = self.bs[k].as_array_mut();
                        let cb = csb[k].as_array_mut();
                        for j in 0..s - 1 {
                            add_comp(full, &mut b[j], &mut cb[j], tmp * C[co + j]);
                        }
                        add_comp(full, &mut b[s - 1], &mut cb[s - 1], tmp);
                        if s == 7 && k < nreal {
                            for c in 0..3 {
                                let ak = at[k][c].abs();
                                if ak.is_normal() && ak > max_ak {
                                    max_ak = ak;
                                }
                                let db = tmp[c].abs();
                                if db.is_normal() && db > max_db6 {
                                    max_db6 = db;
                                }
                            }
                        }
                    }
                    if s == 7 {
                        let e = max_db6 / max_ak;
                        predictor_corrector_error = if e.is_finite() { e } else { 0.0 };
                    }
                }
            }

            // Step-size control.
            let old_timestep = dt;
            let particles_for_error = &state.particles_1;
            let prs = || calculate_new_timestep(particles_for_error, &initial_particles, &self.bs, &old_timestep, &self.epsilon);
            let mut new_timestep = match self.adaptive_mode {
                AdaptiveMode::Prs23 => prs(),
                // The last-term criteria are limited by round-off when the step is tiny and would
                // take many steps to grow out of a small initial step: take the first step size
                // from the PRS23 estimate (the next step's own criterion still checks it).
                _ if self.last_timestep == 0.0 => prs(),
                mode => legacy_new_timestep(mode, &state.particles_1, &initial_particles, &self.bs, old_timestep, self.epsilon),
            };
            if new_timestep.abs() < self.min_timestep {
                new_timestep = self.min_timestep.copysign(new_timestep);
            }
            if (new_timestep / old_timestep).abs() < SAFETY_FACTOR && old_timestep.abs() > self.min_timestep {
                self.retry_with_timestep(new_timestep.copysign(old_timestep), state, &initial_particles, &initial_variational_particles);
                continue 'integration_loop;
            }
            // Do not let the step grow by more than 1/SAFETY_FACTOR at once (as in REBOUND).
            if (new_timestep / old_timestep).abs() > 1.0 / SAFETY_FACTOR {
                new_timestep = old_timestep / SAFETY_FACTOR;
            }
            let new_timestep = new_timestep.abs().copysign(old_timestep);

            // Accept the step: final positions and velocities, smallest terms first.
            let t_end = t_begin + dt;
            for k in 0..n {
                let b = &self.bs[k];
                let (mut x, mut v) = (x0[k], v0[k]);
                let (cx, cv) = (&mut self.csx[k], &mut self.csv[k]);
                if sum == Summation::None {
                    x = x0[k] + dt * v0[k]
                        + dt * dt * (a0[k] / 2.0 + b.p0 / 6.0 + b.p1 / 12.0 + b.p2 / 20.0 + b.p3 / 30.0 + b.p4 / 42.0 + b.p5 / 56.0 + b.p6 / 72.0);
                    v = v0[k] + dt * (a0[k] + b.p0 / 2.0 + b.p1 / 3.0 + b.p2 / 4.0 + b.p3 / 5.0 + b.p4 / 6.0 + b.p5 / 7.0 + b.p6 / 8.0);
                } else {
                    // dt*dt is not precomputed, to avoid biased round-off with a fixed step.
                    for term in [b.p6 / 72.0 * dt * dt, b.p5 / 56.0 * dt * dt, b.p4 / 42.0 * dt * dt, b.p3 / 30.0 * dt * dt,
                                 b.p2 / 20.0 * dt * dt, b.p1 / 12.0 * dt * dt, b.p0 / 6.0 * dt * dt, a0[k] / 2.0 * dt * dt, v0[k] * dt] {
                        add_comp(true, &mut x, cx, term);
                    }
                    for term in [b.p6 / 8.0 * dt, b.p5 / 7.0 * dt, b.p4 / 6.0 * dt, b.p3 / 5.0 * dt, b.p2 / 4.0 * dt,
                                 b.p1 / 3.0 * dt, b.p0 / 2.0 * dt, a0[k] * dt] {
                        add_comp(true, &mut v, cv, term);
                    }
                }
                if k < nreal {
                    let p1 = &mut state.particles_1[k];
                    p1.position = x;
                    p1.velocity = v;
                    p1.epoch = t_end;
                    let p0 = &mut state.particles_0[k];
                    p0.position = initial_particles[k].position;
                    p0.velocity = initial_particles[k].velocity;
                    p0.acceleration = initial_particles[k].acceleration;
                    p0.epoch = t_begin;
                } else {
                    let j = k - nreal;
                    let p1 = &mut state.variational_particles_1[j];
                    p1.position = x;
                    p1.velocity = v;
                    p1.epoch = t_end;
                    let p0 = &mut state.variational_particles_0[j];
                    p0.position = initial_variational_particles[j].position;
                    p0.velocity = initial_variational_particles[j].velocity;
                    p0.acceleration = initial_variational_particles[j].acceleration;
                    p0.epoch = t_begin;
                }
            }
            for (i, sp) in state.spice_particles.iter_mut().enumerate() {
                let st = &substeps[7 * nb + i];
                sp.position = Vector3::new(st[0], st[1], st[2]);
                sp.velocity = Vector3::new(st[3], st[4], st[5]);
                sp.epoch = t_end;
            }

            self.last_timestep = dt;
            self.timestep = new_timestep;
            let ratio = self.timestep / self.last_timestep;
            self.es_last = self.es.clone();
            self.bs_last = self.bs.clone();
            predict_next_coefficients(&ratio, &self.es_last, &self.bs_last, &mut self.es, &mut self.bs);
            break 'integration_loop;
        }
        Ok(())
    }

    fn timestep(&self) -> f64 {
        self.timestep
    }

    fn set_timestep(&mut self, timestep: f64) {
        self.timestep = timestep;
    }

    fn last_timestep(&self) -> f64 {
        self.last_timestep
    }

    fn bs_last(&self) -> &Vec<CoefficientSeptet> {
        &self.bs_last
    }

    fn epsilon(&self) -> f64 {
        self.epsilon
    }

    fn set_epsilon(&mut self, epsilon: f64) {
        self.epsilon = epsilon;
    }

    fn min_timestep(&self) -> f64 {
        self.min_timestep
    }

    fn set_min_timestep(&mut self, min_timestep: f64) {
        self.min_timestep = min_timestep.abs();
    }

    fn adaptive_mode(&self) -> AdaptiveMode {
        self.adaptive_mode
    }

    fn set_adaptive_mode(&mut self, mode: AdaptiveMode) {
        self.adaptive_mode = mode;
    }

    fn summation(&self) -> Summation {
        self.summation
    }

    fn set_summation(&mut self, summation: Summation) {
        self.summation = summation;
    }
}

/// A 7-element coefficient set used by the IAS15 integrator to represent series expansions
/// of position, velocity and acceleration for each body. Each component (p0 through p6) 
/// represents a term in the series approximation (third through ninth order derivatives of position).
#[derive(Clone, Debug)]
#[repr(C)]
pub struct CoefficientSeptet {
    pub p0: Vector3<f64>,
    pub p1: Vector3<f64>,
    pub p2: Vector3<f64>,
    pub p3: Vector3<f64>,
    pub p4: Vector3<f64>,
    pub p5: Vector3<f64>,
    pub p6: Vector3<f64>,
}

/// Add `inp` to `p`; with `compensated`, Kahan's update carries the rounding error in `cs` (the
/// true sum is `p - cs`), as REBOUND's `add_cs`.
#[inline(always)]
fn add_comp(compensated: bool, p: &mut Vector3<f64>, cs: &mut Vector3<f64>, inp: Vector3<f64>) {
    if compensated {
        let y = inp - *cs;
        let t = *p + y;
        *cs = (t - *p) - y;
        *p = t;
    } else {
        *p += inp;
    }
}

impl CoefficientSeptet {
    /// The seven coefficients as an array (the struct is `repr(C)` with seven fields of the
    /// same type, so it has the array's layout).
    #[inline(always)]
    pub fn as_array(&self) -> &[Vector3<f64>; 7] {
        // SAFETY: repr(C) struct of seven Vector3<f64> fields: same size, alignment and field
        // order as [Vector3<f64>; 7], with no padding.
        unsafe { &*(self as *const CoefficientSeptet as *const [Vector3<f64>; 7]) }
    }

    #[inline(always)]
    pub fn as_array_mut(&mut self) -> &mut [Vector3<f64>; 7] {
        // SAFETY: as in `as_array`.
        unsafe { &mut *(self as *mut CoefficientSeptet as *mut [Vector3<f64>; 7]) }
    }

    // fn new(p0: Vector3<f64>, p1: Vector3<f64>, p2: Vector3<f64>, p3: Vector3<f64>, p4: Vector3<f64>, p5: Vector3<f64>, p6: Vector3<f64>) -> CoefficientSeptet {
    //     CoefficientSeptet { p0, p1, p2, p3, p4, p5, p6 }
    // }

    fn zeros() -> CoefficientSeptet {
        CoefficientSeptet { p0: Vector3::zeros(), p1: Vector3::zeros(), p2: Vector3::zeros(), p3: Vector3::zeros(), p4: Vector3::zeros(), p5: Vector3::zeros(), p6: Vector3::zeros() }
    }
}


/// Predicts coefficients for the next integration step based on the current solution
///
/// # Arguments
///
/// * `ratio` - Ratio of new timestep to current timestep
/// * `es_last` - Previous error estimate coefficients
/// * `bs_last` - Previous solution coefficient\]
/// * `es` - Current error estimate coefficients (modified in-place)
/// * `bs` - Current solution coefficients (modified in-place)
///
/// Uses polynomial extrapolation to predict initial values for the next step's
/// coefficients, improving convergence of the predictor-corrector iteration.
fn predict_next_coefficients(ratio: &f64, es_last: &Vec<CoefficientSeptet>, bs_last: &Vec<CoefficientSeptet>, es: &mut Vec<CoefficientSeptet>, bs: &mut Vec<CoefficientSeptet>) {

    let rat = *ratio;

    if rat > 20.0 {
        for e in es.iter_mut() {
            e.p0 = Vector3::zeros();
            e.p1 = Vector3::zeros();
            e.p2 = Vector3::zeros();
            e.p3 = Vector3::zeros();
            e.p4 = Vector3::zeros();
            e.p5 = Vector3::zeros();
            e.p6 = Vector3::zeros();
        }
        for b in bs.iter_mut() {
            b.p0 = Vector3::zeros();
            b.p1 = Vector3::zeros();
            b.p2 = Vector3::zeros();
            b.p3 = Vector3::zeros();
            b.p4 = Vector3::zeros();
            b.p5 = Vector3::zeros();
            b.p6 = Vector3::zeros();
        }
    } else {
        let q1 = rat;
        let q2 = q1.powi(2);
        let q3 = q1 * q2;
        let q4 = q2.powi(2);
        let q5 = q2 * q3;
        let q6 = q3.powi(2);
        let q7 = q3 * q4;

        for idx in 0..es.len() {
            let e = &mut es[idx];
            let b = &mut bs[idx];
            let e_last = &es_last[idx];
            let b_last = &bs_last[idx];

            let be0 = b_last.p0 - e_last.p0;
            let be1 = b_last.p1 - e_last.p1;
            let be2 = b_last.p2 - e_last.p2;
            let be3 = b_last.p3 - e_last.p3;
            let be4 = b_last.p4 - e_last.p4;
            let be5 = b_last.p5 - e_last.p5;
            let be6 = b_last.p6 - e_last.p6;

            e.p0 = q1 * (b_last.p6 * 7.0 + b_last.p5 * 6.0 + b_last.p4 * 5.0 + b_last.p3 * 4.0 + b_last.p2 * 3.0 + b_last.p1 * 2.0 + b_last.p0);
            e.p1 = q2 * (b_last.p6 * 21.0 + b_last.p5 * 15.0 + b_last.p4 * 10.0 + b_last.p3 * 6.0 + b_last.p2 * 3.0 + b_last.p1);
            e.p2 = q3 * (b_last.p6 * 35.0 + b_last.p5 * 20.0 + b_last.p4 * 10.0 + b_last.p3 * 4.0 + b_last.p2);
            e.p3 = q4 * (b_last.p6 * 35.0 + b_last.p5 * 15.0 + b_last.p4 * 5.0 + b_last.p3);
            e.p4 = q5 * (b_last.p6 * 21.0 + b_last.p5 * 6.0 + b_last.p4);
            e.p5 = q6 * (b_last.p6 * 7.0 + b_last.p5);
            e.p6 = q7 * b_last.p6;

            b.p0 = e.p0 + be0;
            b.p1 = e.p1 + be1;
            b.p2 = e.p2 + be2;
            b.p3 = e.p3 + be3;
            b.p4 = e.p4 + be4;
            b.p5 = e.p5 + be5;
            b.p6 = e.p6 + be6;

        }
    }

}


/// Calculates the optimal timestep for the next integration step.
///
/// # Arguments
///
/// * `particles` - Vector of particles being integrated
/// * `accelerations` - Current accelerations for all particles
/// * `bs` - Current solution coefficients
/// * `last_timestep` - Previous integration timestep
/// * `epsilon` - Desired integration accuracy
///
/// # Returns
///
/// * `f64` - Recommended timestep for the next integration step
///
/// Estimates the optimal timestep based on the local truncation error and 
/// the dynamics of the system. Uses the ratio of successive terms in the
/// series expansion to gauge the convergence rate.
///
/// The timescale is determined by examining the magnitude of acceleration terms 
/// and their derivatives. If the estimated error is too large, the timestep 
/// will be reduced. If the error is well below the tolerance, the timestep may be 
/// increased to improve efficiency.
pub fn calculate_new_timestep(particles: &Vec<SimulationParticle>, initial_particles: &Vec<SimulationParticle>, bs: &Vec<CoefficientSeptet>, last_timestep: &f64, epsilon: &f64) -> f64 {
    let mut min_timescale2 = f64::INFINITY;
    for idx in 0..particles.len() {
        // let particle = &particles[idx];
        let b = &bs[idx];
        let a0 = initial_particles[idx].acceleration.norm_squared();
        let y2 = (initial_particles[idx].acceleration + b.p0 + b.p1 + b.p2 + b.p3 + b.p4 + b.p5 + b.p6).norm_squared();
        let y3 = (b.p0 + 2.0 * b.p1 + 3.0 * b.p2 + 4.0 * b.p3 + 5.0 * b.p4 + 6.0 * b.p5 + 7.0 * b.p6).norm_squared();
        let y4 = (2.0 * b.p1 + 6.0 * b.p2 + 12.0 * b.p3 + 20.0 * b.p4 + 30.0 * b.p5 + 42.0 * b.p6).norm_squared();

        if !a0.is_normal() {
            continue;
        }

        let timescale2 = 2.0 * y2 / (y3 + (y4 * y2).sqrt());
        if (timescale2 < min_timescale2) & timescale2.is_normal() {
            min_timescale2 = timescale2;
        }
    }

    if min_timescale2.is_normal() {
        min_timescale2.sqrt() * last_timestep * (epsilon * 5040.0).powf(1.0 / 7.0)
    } else {
        last_timestep / SAFETY_FACTOR
    }

}


/// REBOUND's original timestep criteria (`Global` and `Individual`): the relative size of the
/// last term of the acceleration series sets the step, `dt_new = (epsilon / error)^(1/7) dt`.
pub fn legacy_new_timestep(
    mode: AdaptiveMode,
    particles: &[SimulationParticle],
    initial_particles: &[SimulationParticle],
    bs: &[CoefficientSeptet],
    last_timestep: f64,
    epsilon: f64,
) -> f64 {
    let mut error = 0.0f64;
    match mode {
        AdaptiveMode::Global => {
            let (mut maxa, mut maxj) = (0.0f64, 0.0f64);
            for (idx, p) in particles.iter().enumerate() {
                let v2 = initial_particles[idx].velocity.norm_squared();
                let x2 = initial_particles[idx].position.norm_squared();
                // Skip slowly varying accelerations.
                if (v2 * last_timestep * last_timestep / x2).abs() < 1e-16 {
                    continue;
                }
                for k in 0..3 {
                    let ak = p.acceleration[k].abs();
                    if ak.is_normal() && ak > maxa {
                        maxa = ak;
                    }
                    let b6k = bs[idx].p6[k].abs();
                    if b6k.is_normal() && b6k > maxj {
                        maxj = b6k;
                    }
                }
            }
            error = maxj / maxa;
        }
        AdaptiveMode::PerParticle => {
            for (idx, p) in particles.iter().enumerate() {
                let (mut maxa, mut maxj) = (0.0f64, 0.0f64);
                for k in 0..3 {
                    maxa = maxa.max(p.acceleration[k].abs());
                    maxj = maxj.max(bs[idx].p6[k].abs());
                }
                let e = maxj / maxa;
                if e.is_normal() && e > error {
                    error = e;
                }
            }
        }
        _ => {
            for (idx, p) in particles.iter().enumerate() {
                for k in 0..3 {
                    let e = (bs[idx].p6[k] / p.acceleration[k]).abs();
                    if e.is_normal() && e > error {
                        error = e;
                    }
                }
            }
        }
    }
    if error.is_normal() {
        (epsilon / error).powf(1.0 / 7.0) * last_timestep
    } else {
        last_timestep / SAFETY_FACTOR
    }
}

/// Result of one integrator step.
pub type StepResult = Result<(), Box<dyn std::error::Error + Send + Sync>>;

/// Recompute accelerations (and, when variational particles exist, the partials and the
/// variational accelerations).
pub fn update_forces(state: &mut SimulationState, forces: &Vec<Box<dyn Force + Send + Sync>>, with_stm: bool) {
    if with_stm {
        update_accelerations_and_stms(state, forces);
        update_variational_accelerations(state);
    } else {
        update_accelerations(state, forces);
    }
}

pub fn update_accelerations(state: &mut SimulationState, forces: &Vec<Box<dyn Force + Send + Sync>>) {
    // first clear the accelerations vector
    for acc in state.particles_1.iter_mut() {
        acc.acceleration = Vector3::zeros();
    }
    for force in forces {
        force.apply_acceleration(state);
    }
}


pub fn update_accelerations_and_stms(state: &mut SimulationState, forces: &Vec<Box<dyn Force + Send + Sync>>) {
    for p in state.particles_1.iter_mut() {
        p.acceleration = Vector3::zeros();
    }
    state.partials.clear();
    state.partials.resize(state.particles_1.len(), Partials::default());
    for force in forces {
        force.apply_acceleration_and_stm(state);
    }
}

/// Accelerations of the variational particles from their parents' acceleration Jacobians:
/// da = (da/dr) dr + (da/dv) dv + (da/dA) dA.
pub fn update_variational_accelerations(state: &mut SimulationState) {
    for p in state.variational_particles_1.iter_mut() {
        let parent = &state.partials[p.parent];
        let j = &parent.stm;
        let g = &parent.nongrav;
        let (r, v, dk) = (p.position, p.velocity, p.nongrav);
        let mut a = [0.0; 3];
        for i in 0..3 {
            let row = &j[3 + i];
            a[i] = row[0] * r.x + row[1] * r.y + row[2] * r.z
                + row[3] * v.x + row[4] * v.y + row[5] * v.z
                + g[i][0] * dk[0] + g[i][1] * dk[1] + g[i][2] * dk[2];
        }
        p.acceleration = Vector3::new(a[0], a[1], a[2]);
    }
}
