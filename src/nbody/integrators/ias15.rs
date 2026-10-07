use crate::nbody::integrators::Integrator;
use crate::nbody::forces::{total_acceleration, Force};
use crate::state::{from_pv, pv, State};

use nalgebra::Vector3;

// Gauss Radau spacings
const H: [f64; 8] = [0.0, 0.056_262_560_536_922_15, 0.180_240_691_736_892_36, 0.352_624_717_113_169_6, 
                     0.547_153_626_330_555_4, 0.734_210_177_215_410_5, 0.885_320_946_839_095_8, 0.977_520_613_561_287_5];
// Other constants
const RR: [f64; 28] = [0.056_262_560_536_922_15, 0.180_240_691_736_892_36, 0.123_978_131_199_970_21, 0.352_624_717_113_169_6, 
                       0.296_362_156_576_247_5, 0.172_384_025_376_277_28, 0.547_153_626_330_555_4, 0.490_891_065_793_633_23, 
                       0.366_912_934_593_663_03, 0.194_528_909_217_385_75, 0.734_210_177_215_410_5, 0.677_947_616_678_488_4, 
                       0.553_969_485_478_518_2, 0.381_585_460_102_240_87, 0.187_056_550_884_855_15, 0.885_320_946_839_095_8, 
                       0.829_058_386_302_173_7, 0.705_080_255_102_203_4, 0.532_696_229_725_926_1, 0.338_167_320_508_540_37, 
                       0.151_110_769_623_685_25, 0.977_520_613_561_287_5, 0.921_258_053_024_365_4, 0.797_279_921_824_395_1, 
                       0.624_895_896_448_117_8, 0.430_366_987_230_732_13, 0.243_310_436_345_876_96, 0.092_199_666_722_191_74];

const C: [f64; 21] = [-0.056_262_560_536_922_15, 0.010_140_802_830_063_63, -0.236_503_252_273_814_52, -0.003_575_897_729_251_617_6, 
                       0.093_537_695_259_462_07, -0.589_127_969_386_984_2, 0.001_956_565_409_947_221, -0.054_755_386_889_068_69, 
                       0.415_881_200_082_306_83, -1.136_281_595_717_539_6, -0.001_436_530_236_370_891_5, 0.042_158_527_721_268_706, 
                       -0.360_099_596_502_056_8, 1.250_150_711_840_691, -1.870_491_772_932_95, 0.001_271_790_309_026_867_8, 
                       -0.038_760_357_915_906_77, 0.360_962_243_452_846, -1.466_884_208_400_427, 2.906_136_259_308_429_4, 
                       -2.755_812_719_772_045_7];

const D: [f64; 21] = [0.056_262_560_536_922_15, 0.003_165_475_718_170_829_3, 0.236_503_252_273_814_52, 0.000_178_097_769_221_743_38, 
                      0.045_792_985_506_027_92, 0.589_127_969_386_984_2, 0.000_010_020_236_522_329_128, 0.008_431_857_153_525_702, 
                      0.253_534_069_054_569_27, 1.136_281_595_717_539_6, 0.000_000_563_764_163_931_820_8, 0.001_529_784_002_500_465_7, 
                      0.097_834_236_532_444_01, 0.875_254_664_684_091_1, 1.870_491_772_932_95, 0.000_000_031_718_815_401_761_364, 
                      0.000_276_293_090_982_647_7, 0.036_028_553_983_736_46, 0.576_733_000_277_078_7, 2.248_588_760_769_16, 
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
    bs: Vec<CoefficientSeptet>,
    /// Intermediate step coefficients
    gs: Vec<CoefficientSeptet>,
    /// Error estimate coefficients
    es: Vec<CoefficientSeptet>,
    /// Previous step coefficients
    bs_last: Vec<CoefficientSeptet>,
    /// Previous error estimate coefficients
    es_last: Vec<CoefficientSeptet>,
    /// The last completed step, for dense output
    dense: Option<DenseOutput>,
    /// Scratch buffer for the accelerations at each substep
    accelerations: Vec<Vector3<f64>>,
}

/// The start of the last completed step. With that step's coefficients (`bs_last`) it gives the
/// trajectory anywhere within the step.
#[derive(Debug, Clone)]
struct DenseOutput {
    /// Time at the start of the step
    t0: f64,
    /// Length of the step
    dt: f64,
    x0: Vec<Vector3<f64>>,
    v0: Vec<Vector3<f64>>,
    a0: Vec<Vector3<f64>>,
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
        IAS15 { timestep, epsilon: 1e-9, last_timestep: 0.0, bs: vec![], gs: vec![], es: vec![], bs_last: vec![], es_last: vec![], dense: None, accelerations: vec![] }
    }

    /// Resets all coefficient vectors to zero for the specified number of particles
    ///
    /// # Arguments
    ///
    /// * `n` - Number of particles in the simulation
    ///
    /// This is called when the number of particles changes or at the start of integration
    /// to ensure proper sizing of all coefficient vectors.
    pub fn reset_coefficients(&mut self, n: usize) {
        self.bs = vec![CoefficientSeptet::zeros(); n];
        self.gs = vec![CoefficientSeptet::zeros(); n];
        self.es = vec![CoefficientSeptet::zeros(); n];
        self.bs_last = vec![CoefficientSeptet::zeros(); n];
        self.es_last = vec![CoefficientSeptet::zeros(); n];
        self.dense = None;
    }
}

impl Integrator for IAS15 {

    /// Advances the system by one step using the IAS15 algorithm.
    fn step(&mut self, states: &mut [State], masses: &[f64], t: &mut f64, forces: &[Box<dyn Force + Send + Sync>]) {
        // Number of particles
        let n = states.len();

        // The start of the step, in the buffers of the previous step's dense output.
        let (mut initial_positions, mut initial_velocities, mut initial_accelerations) = match self.dense.take() {
            Some(d) => (d.x0, d.v0, d.a0),
            None => (Vec::new(), Vec::new(), Vec::new()),
        };
        initial_positions.clear();
        initial_positions.extend(states.iter().map(|p| pv(p).0));
        initial_velocities.clear();
        initial_velocities.extend(states.iter().map(|p| pv(p).1));
        total_acceleration(forces, states, masses, &mut initial_accelerations);

        if (self.bs.len() != n) || (self.gs.len() != n) {
            self.reset_coefficients(n);
        }
      
        for (g, b) in self.gs.iter_mut().zip(self.bs.iter()) {
            g.p0 = b.p6 * D[15] + b.p5 * D[10] + b.p4 * D[6] + b.p3 * D[3] + b.p2 * D[1] + b.p1 * D[0] + b.p0;
            g.p1 = b.p6 * D[16] + b.p5 * D[11] + b.p4 * D[7] + b.p3 * D[4] + b.p2 * D[2] + b.p1;
            g.p2 = b.p6 * D[17] + b.p5 * D[12] + b.p4 * D[8] + b.p3 * D[5] + b.p2;
            g.p3 = b.p6 * D[18] + b.p5 * D[13] + b.p4 * D[9] + b.p3;
            g.p4 = b.p6 * D[19] + b.p5 * D[14] + b.p4;
            g.p5 = b.p6 * D[20] + b.p5;
            g.p6 = b.p6;
        }

        // Multiplying by these is cheaper than dividing by RR for every particle and substep.
        let rri = RR.map(|r| 1.0 / r);
        let mut accelerations = std::mem::take(&mut self.accelerations);

        let mut predictor_corrector_error = 1e300;
        let mut predictor_corrector_error_last = 2.0;
        let mut iterations = 0;

        // This is the predictor-corrector loop, which calculates the coefficients for the next step
        loop {

            if predictor_corrector_error < 1e-16 {
                break;
            }
            if iterations > 2 && predictor_corrector_error_last <= predictor_corrector_error {
                break;
            }
            // Like REBOUND, give up on convergence after 12 iterations and let the timestep control
            // below decide whether to accept the step.
            if iterations >= 12 {
                break;
            }

            predictor_corrector_error_last = predictor_corrector_error;
            predictor_corrector_error = 0.0;
            iterations += 1;

            for substep in 1..8 {
                let hh = H[substep];
                // Factors of the nested position and velocity series at this substep.
                let xf = [7.0 * hh / 9.0, 3.0 * hh / 4.0, 5.0 * hh / 7.0, 2.0 * hh / 3.0, 3.0 * hh / 5.0, hh / 2.0, hh / 3.0, self.timestep * hh / 2.0, self.timestep * hh];
                let vf = [7.0 * hh / 8.0, 6.0 * hh / 7.0, 5.0 * hh / 6.0, 4.0 * hh / 5.0, 3.0 * hh / 4.0, 2.0 * hh / 3.0, hh / 2.0, self.timestep * hh];
                for idx in 0..n {
                    let a0 = initial_accelerations[idx];
                    let v0 = initial_velocities[idx];
                    let b = &self.bs[idx];

                    // Calculate the position
                    let d_position = ((((((((b.p6 * xf[0] + b.p5) * xf[1] + b.p4) * xf[2] + b.p3) * xf[3] + b.p2) * xf[4] + b.p1) * xf[5] + b.p0) * xf[6] + a0) * xf[7] + v0) * xf[8];

                    // Calculate the velocity
                    let d_velocity = (((((((b.p6 * vf[0] + b.p5) * vf[1] + b.p4) * vf[2] + b.p3) * vf[3] + b.p2) * vf[4] + b.p1) * vf[5] + b.p0) * vf[6] + a0) * vf[7];
                    states[idx] = from_pv(&(initial_positions[idx] + d_position), &(v0 + d_velocity));
                }

                total_acceleration(forces, states, masses, &mut accelerations);

                match substep {
                    1 => {
                        for idx in 0..n {
                            let a_old = initial_accelerations[idx];
                            let a_new = accelerations[idx];

                            let temp = self.gs[idx].p0;

                            self.gs[idx].p0 = (a_new - a_old) * rri[0];
                            self.bs[idx].p0 += self.gs[idx].p0 - temp;
                        }
                    },
                    2 => {
                        for idx in 0..n {
                            let a_old = initial_accelerations[idx];
                            let a_new = accelerations[idx];

                            let mut temp = self.gs[idx].p1;
                            self.gs[idx].p1 = ((a_new - a_old) * rri[1] - self.gs[idx].p0) * rri[2];
                            temp = self.gs[idx].p1 - temp;

                            self.bs[idx].p0 += temp * C[0];
                            self.bs[idx].p1 += temp;
                        }
                    },
                    3 => {
                        for idx in 0..n {
                            let a_old = initial_accelerations[idx];
                            let a_new = accelerations[idx];

                            let mut temp = self.gs[idx].p2;
                            self.gs[idx].p2 = (((a_new - a_old) * rri[3] - self.gs[idx].p0) * rri[4] - self.gs[idx].p1) * rri[5];
                            temp = self.gs[idx].p2 - temp;

                            self.bs[idx].p0 += temp * C[1];
                            self.bs[idx].p1 += temp * C[2];
                            self.bs[idx].p2 += temp;
                        }
                    },
                    4 => {
                        for idx in 0..n {
                            let a_old = initial_accelerations[idx];
                            let a_new = accelerations[idx];

                            let mut temp = self.gs[idx].p3;
                            self.gs[idx].p3 = ((((a_new - a_old) * rri[6] - self.gs[idx].p0) * rri[7] - self.gs[idx].p1) * rri[8] - self.gs[idx].p2) * rri[9];
                            temp = self.gs[idx].p3 - temp;

                            self.bs[idx].p0 += temp * C[3];
                            self.bs[idx].p1 += temp * C[4];
                            self.bs[idx].p2 += temp * C[5];
                            self.bs[idx].p3 += temp;
                        }
                    },
                    5 => {
                        for idx in 0..n {
                            let a_old = initial_accelerations[idx];
                            let a_new = accelerations[idx];

                            let mut temp = self.gs[idx].p4;
                            self.gs[idx].p4 = (((((a_new - a_old) * rri[10] - self.gs[idx].p0) * rri[11] - self.gs[idx].p1) * rri[12] - self.gs[idx].p2) * rri[13] - self.gs[idx].p3) * rri[14];
                            temp = self.gs[idx].p4 - temp;

                            self.bs[idx].p0 += temp * C[6];
                            self.bs[idx].p1 += temp * C[7];
                            self.bs[idx].p2 += temp * C[8];
                            self.bs[idx].p3 += temp * C[9];
                            self.bs[idx].p4 += temp;
                        }
                    },
                    6 => {
                        for idx in 0..n {
                            let a_old = initial_accelerations[idx];
                            let a_new = accelerations[idx];

                            let mut temp = self.gs[idx].p5;
                            self.gs[idx].p5 = ((((((a_new - a_old) * rri[15] - self.gs[idx].p0) * rri[16] - self.gs[idx].p1) * rri[17] - self.gs[idx].p2) * rri[18] - self.gs[idx].p3) * rri[19] - self.gs[idx].p4) * rri[20];
                            temp = self.gs[idx].p5 - temp;

                            self.bs[idx].p0 += temp * C[10];
                            self.bs[idx].p1 += temp * C[11];
                            self.bs[idx].p2 += temp * C[12];
                            self.bs[idx].p3 += temp * C[13];
                            self.bs[idx].p4 += temp * C[14];
                            self.bs[idx].p5 += temp;
                        }
                    },
                    7 => {
                        // The error is the largest b6 correction over the largest acceleration.
                        let mut max_acceleration2: f64 = 0.0;
                        let mut max_b6_temp2: f64 = 0.0;
                        for idx in 0..n {
                            let a_old = initial_accelerations[idx];
                            let a_new = accelerations[idx];

                            let mut temp = self.gs[idx].p6;
                            self.gs[idx].p6 = (((((((a_new - a_old) * rri[21] - self.gs[idx].p0) * rri[22] - self.gs[idx].p1) * rri[23] - self.gs[idx].p2) * rri[24] - self.gs[idx].p3) * rri[25] - self.gs[idx].p4) * rri[26] - self.gs[idx].p5) * rri[27];
                            temp = self.gs[idx].p6 - temp;

                            self.bs[idx].p0 += temp * C[15];
                            self.bs[idx].p1 += temp * C[16];
                            self.bs[idx].p2 += temp * C[17];
                            self.bs[idx].p3 += temp * C[18];
                            self.bs[idx].p4 += temp * C[19];
                            self.bs[idx].p5 += temp * C[20];
                            self.bs[idx].p6 += temp;

                            let temp2 = temp.norm_squared();
                            if temp2 > max_b6_temp2 && temp2.is_normal() {
                                max_b6_temp2 = temp2;
                            }
                            let a_new2 = a_new.norm_squared();
                            if a_new2 > max_acceleration2 && a_new2.is_normal() {
                                max_acceleration2 = a_new2;
                            }
                        }
                        let error = (max_b6_temp2 / max_acceleration2).sqrt();
                        if error.is_normal() {
                            predictor_corrector_error = error;
                        }
                    },
                    _ => {}
                }
            }
        }
        self.accelerations = accelerations;

        let old_timestep = self.timestep;
        let mut new_timestep = calculate_new_timestep(&initial_accelerations, &self.bs, &old_timestep, &self.epsilon);
        let timestep_ratio = (new_timestep / old_timestep).abs();
        

        // Step was rejected
        if timestep_ratio < SAFETY_FACTOR {
            // println!("Timestep was rejected. Reducing the timestep to {}", new_timestep);
            self.timestep = new_timestep;

            // reset particles
            for idx in 0..n {
                states[idx] = from_pv(&initial_positions[idx], &initial_velocities[idx]);
            }

            if self.last_timestep != 0.0 {
                let ratio = self.timestep / self.last_timestep;
                predict_next_coefficients(&ratio, &self.es_last, &self.bs_last, &mut self.es, &mut self.bs);
            }

            // Redo the step with the new timestep. The retry advances the particles and the
            // epoch, so this attempt must not.
            self.step(states, masses, t, forces);
            return;
        }

        // The timestep was accepted
        if timestep_ratio > 1.0 {
            // println!("Timestep was accepted. Increasing the timestep to {}", new_timestep);
            if timestep_ratio > 1.0 / SAFETY_FACTOR {
                // println!("The timestep ratio is greater than 1/Safety Factor. Increasing the timestep to {}", old_timestep / SAFETY_FACTOR);
                new_timestep = old_timestep / SAFETY_FACTOR;
            }
        }

        

        let t0 = *t;

        // Update the time
        *t += self.timestep;

        // Update the particles
        for idx in 0..n {
            let b = &self.bs[idx];
            let position = initial_positions[idx] + self.timestep * initial_velocities[idx] + self.timestep.powi(2) * (initial_accelerations[idx] / 2.0 + b.p0 / 6.0 + b.p1 / 12.0 + b.p2 / 20.0 + b.p3 / 30.0 + b.p4 / 42.0 + b.p5 / 56.0 + b.p6 / 72.0);
            let velocity = initial_velocities[idx] + self.timestep * (initial_accelerations[idx] + b.p0 / 2.0 + b.p1 / 3.0 + b.p2 / 4.0 + b.p3 / 5.0 + b.p4 / 6.0 + b.p5 / 7.0 + b.p6 / 8.0);
            states[idx] = from_pv(&position, &velocity);
        }

        

        self.dense = Some(DenseOutput { t0, dt: self.timestep, x0: initial_positions, v0: initial_velocities, a0: initial_accelerations });

        self.last_timestep = self.timestep.clone();
        self.timestep = new_timestep;
        let ratio = self.timestep / self.last_timestep;


        self.es_last.clone_from(&self.es);
        self.bs_last.clone_from(&self.bs);

        predict_next_coefficients(&ratio, &self.es_last, &self.bs_last, &mut self.es, &mut self.bs);        

    }

    fn timestep(&self) -> f64 {
        self.timestep
    }

    fn set_timestep(&mut self, timestep: f64) {
        // The predicted coefficients are scaled to the step they were predicted for; predict
        // them again for the new one, or the predictor-corrector starts from a poor guess
        // (shortening a step to land on an output would otherwise take up to 12 iterations).
        if timestep != self.timestep && self.last_timestep != 0.0 && self.bs_last.len() == self.bs.len() {
            let ratio = timestep / self.last_timestep;
            predict_next_coefficients(&ratio, &self.es_last, &self.bs_last, &mut self.es, &mut self.bs);
        }
        self.timestep = timestep;
    }

    fn has_dense_output(&self) -> bool {
        true
    }

    /// Evaluates the last step's Gauss–Radau polynomial at `t`, as REBOUND and ASSIST do. The
    /// result is as accurate as the step itself.
    fn interpolate(&self, t: f64) -> Option<Vec<State>> {
        let d = self.dense.as_ref()?;
        let h = (t - d.t0) / d.dt;
        // The tolerance absorbs round-off in Julian dates at the ends of the step.
        if !(-1e-8..=1.0 + 1e-8).contains(&h) || self.bs_last.len() != d.x0.len() {
            return None;
        }

        // s[k] and u[k] multiply the position and velocity series terms.
        let mut s = [0.0; 9];
        s[0] = d.dt * h;
        s[1] = s[0] * s[0] / 2.0;
        s[2] = s[1] * h / 3.0;
        s[3] = s[2] * h / 2.0;
        s[4] = 3.0 * s[3] * h / 5.0;
        s[5] = 2.0 * s[4] * h / 3.0;
        s[6] = 5.0 * s[5] * h / 7.0;
        s[7] = 3.0 * s[6] * h / 4.0;
        s[8] = 7.0 * s[7] * h / 9.0;

        let mut u = [0.0; 8];
        u[0] = d.dt * h;
        u[1] = u[0] * h / 2.0;
        u[2] = 2.0 * u[1] * h / 3.0;
        u[3] = 3.0 * u[2] * h / 4.0;
        u[4] = 4.0 * u[3] * h / 5.0;
        u[5] = 5.0 * u[4] * h / 6.0;
        u[6] = 6.0 * u[5] * h / 7.0;
        u[7] = 7.0 * u[6] * h / 8.0;

        let mut states = Vec::with_capacity(d.x0.len());
        for (idx, b) in self.bs_last.iter().enumerate() {
            let (x0, v0, a0) = (d.x0[idx], d.v0[idx], d.a0[idx]);
            let position = x0 + s[8] * b.p6 + s[7] * b.p5 + s[6] * b.p4 + s[5] * b.p3 + s[4] * b.p2 + s[3] * b.p1 + s[2] * b.p0 + s[1] * a0 + s[0] * v0;
            let velocity = v0 + u[7] * b.p6 + u[6] * b.p5 + u[5] * b.p4 + u[4] * b.p3 + u[3] * b.p2 + u[2] * b.p1 + u[1] * b.p0 + u[0] * a0;
            states.push(from_pv(&position, &velocity));
        }
        Some(states)
    }
}


/// A 7-element coefficient set used by the IAS15 integrator to represent series expansions
/// of position, velocity and acceleration for each body. Each component (p0 through p6) 
/// represents a term in the series approximation (third through ninth order derivatives of position).
#[derive(Clone, Debug)]
pub struct CoefficientSeptet {
    pub p0: Vector3<f64>,
    pub p1: Vector3<f64>,
    pub p2: Vector3<f64>,
    pub p3: Vector3<f64>,
    pub p4: Vector3<f64>,
    pub p5: Vector3<f64>,
    pub p6: Vector3<f64>,
}

impl CoefficientSeptet {
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

    if rat.abs() > 20.0 {
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

pub fn calculate_new_timestep(accelerations: &[Vector3<f64>], bs: &[CoefficientSeptet], last_timestep: &f64, epsilon: &f64) -> f64 {
    let mut min_timescale2 = f64::INFINITY;
    for idx in 0..accelerations.len() {
        // let particle = &particles[idx];
        let b = &bs[idx];
        let a0 = accelerations[idx].norm_squared();
        let y2 = (accelerations[idx] + b.p0 + b.p1 + b.p2 + b.p3 + b.p4 + b.p5 + b.p6).norm_squared();
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