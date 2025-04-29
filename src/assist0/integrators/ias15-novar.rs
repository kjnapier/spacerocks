use crate::time::Time;
use crate::spice::{SpiceBody, SpiceKernel};

use crate::assist::{SimulationParticle, SimulationState};
use crate::assist::integrators::Integrator;
use crate::assist::forces::Force;

use std::collections::HashMap;
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
        IAS15 { timestep, epsilon: 1e-9, last_timestep: 0.0, bs: vec![], gs: vec![], es: vec![], bs_last: vec![], es_last: vec![] }
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
    }

}

impl Integrator for IAS15 {
    fn step(&mut self, state: &mut SimulationState, forces: &Vec<Box<dyn Force + Send + Sync>>, kernel: &SpiceKernel) {
        // Label the outer integration loop to avoid recursion.
        'integration_loop: loop {
            let n = state.particles.len();

            // Check if the number of particles has changed.
            // if state.spice_particles[0].epoch != state.particles_1[0].epoch {
            //     let state_vectors = kernel.get_barycentric_states(&state.spice_bodies, state.particles_1[0].epoch).unwrap();
            //     for (i, state_vector) in state_vectors.iter().enumerate() {
            //         let (x, y, z, vx, vy, vz) = state_vector;
            //         let position = Vector3::new(*x, *y, *z);
            //         let velocity = Vector3::new(*vx, *vy, *vz);
            //         state.spice_particles[i].position = position;
            //         state.spice_particles[i].velocity = velocity;
            //         state.spice_particles[i].epoch = state.particles_1[0].epoch;
            //     }
            // }
            
            let state_vectors = kernel.get_barycentric_states(&state.spice_bodies, state.particles_1[0].epoch).unwrap();
            for (i, state_vector) in state_vectors.iter().enumerate() {
                let (x, y, z, vx, vy, vz) = state_vector;
                let position = Vector3::new(*x, *y, *z);
                let velocity = Vector3::new(*vx, *vy, *vz);
                state.spice_particles[i].position = position;
                state.spice_particles[i].velocity = velocity;
            }

            update_accelerations(state, forces);

            // Save initial conditions.
            let initial_positions: Vec<Vector3<f64>> = state.particles_1.iter().map(|p| p.position).collect();
            let initial_velocities: Vec<Vector3<f64>> = state.particles_1.iter().map(|p| p.velocity).collect();
            let initial_accelerations: Vec<Vector3<f64>> = state.particles_1.iter().map(|p| p.acceleration).collect();
            let initial_epoch = state.particles_1[0].epoch;

            // Ensure coefficient vectors are sized correctly.
            if self.bs.len() != n || self.gs.len() != n {
                self.reset_coefficients(n);
            }

            // Update intermediate coefficients.
            for (g, b) in self.gs.iter_mut().zip(&self.bs) {
                g.p0 = b.p6 * D[15] + b.p5 * D[10] + b.p4 * D[6] + b.p3 * D[3] + b.p2 * D[1] + b.p1 * D[0] + b.p0;
                g.p1 = b.p6 * D[16] + b.p5 * D[11] + b.p4 * D[7] + b.p3 * D[4] + b.p2 * D[2] + b.p1;
                g.p2 = b.p6 * D[17] + b.p5 * D[12] + b.p4 * D[8] + b.p3 * D[5] + b.p2;
                g.p3 = b.p6 * D[18] + b.p5 * D[13] + b.p4 * D[9] + b.p3;
                g.p4 = b.p6 * D[19] + b.p5 * D[14] + b.p4;
                g.p5 = b.p6 * D[20] + b.p5;
                g.p6 = b.p6;
            }

            let mut predictor_corrector_error = 1e300;
            let mut predictor_corrector_error_last = 2.0;
            let mut iterations = 0;

          
            let mut substep_map: Vec<Vec<(Vector3<f64>, Vector3<f64>)>> = vec![Vec::new(); 8];
            for substep in 1..8 {
                let hh = H[substep];
                let state_vectors = kernel.get_barycentric_states(&state.spice_bodies, initial_epoch + self.timestep * hh).unwrap();
                for (i, state_vector) in state_vectors.iter().enumerate() {
                    let (x, y, z, vx, vy, vz) = state_vector;
                    let position = Vector3::new(*x, *y, *z);
                    let velocity = Vector3::new(*vx, *vy, *vz);
                    substep_map[substep].push((position, velocity));
                }
            }           

            // Predictor-corrector iteration loop.
            loop {

                // if predictor_corrector_error < 1e-16 {
                //     break;
                // }
                if predictor_corrector_error < 2e-16 {
                    break;
                }

                if (iterations > 2 && predictor_corrector_error_last <= predictor_corrector_error) {
                    break;
                }

                if iterations >= 10 {
                    println!("Convergence failed after 10 iterations. Reducing timestep.");
                    self.timestep /= 2.0;
                    continue 'integration_loop;
                }
                predictor_corrector_error_last = predictor_corrector_error;
                predictor_corrector_error = 0.0;
                iterations += 1;

                // Loop over substeps 1..7.
                for substep in 1..8 {

                    let hh = H[substep];
                    for idx in 0..n {
                        let a0 = initial_accelerations[idx];
                        let v0 = initial_velocities[idx];
                        let b = &self.bs[idx];

                        // Compute position increment.
                        let d_position = ((((((((b.p6 * 7.0 * hh / 9.0 + b.p5) * 3.0 * hh / 4.0 + b.p4)
                            * 5.0 * hh / 7.0 + b.p3) * 2.0 * hh / 3.0 + b.p2)
                            * 3.0 * hh / 5.0 + b.p1) * hh / 2.0 + b.p0)
                            * hh / 3.0 + a0) * self.timestep * hh / 2.0 + v0)
                            * self.timestep * hh;
                        state.particles_1[idx].position = initial_positions[idx] + d_position;

                        // Compute velocity increment.
                        let d_velocity = (((((((b.p6 * 7.0 * hh / 8.0 + b.p5) * 6.0 * hh / 7.0 + b.p4)
                            * 5.0 * hh / 6.0 + b.p3) * 4.0 * hh / 5.0 + b.p2)
                            * 3.0 * hh / 4.0 + b.p1) * 2.0 * hh / 3.0 + b.p0)
                            * hh / 2.0 + a0) * self.timestep * hh;
                        state.particles_1[idx].velocity = initial_velocities[idx] + d_velocity;

                        state.particles_1[idx].epoch = initial_epoch + self.timestep * hh;

                    }

                    // // Update the state of the simulation with the new positions and velocities.
                    for (idx, particle) in state.spice_particles.iter_mut().enumerate() {
                        let sm = &substep_map[substep];
                        let (position, velocity) = sm[idx];
                        particle.position = position;
                        particle.velocity = velocity;
                        particle.epoch = state.particles_1[0].epoch;
                    }

                    update_accelerations(state, forces);

                    // Update coefficients based on substep.
                    match substep {
                        1 => {
                            for idx in 0..n {
                                let a_old = initial_accelerations[idx];
                                let a_new = state.particles_1[idx].acceleration;
                                let temp = self.gs[idx].p0;
                                self.gs[idx].p0 = (a_new - a_old) / RR[0];
                                self.bs[idx].p0 += self.gs[idx].p0 - temp;
                            }
                        },
                        2 => {
                            for idx in 0..n {
                                let a_old = initial_accelerations[idx];
                                let a_new = state.particles_1[idx].acceleration;
                                let mut temp = self.gs[idx].p1;
                                self.gs[idx].p1 = ((a_new - a_old) / RR[1] - self.gs[idx].p0) / RR[2];
                                temp = self.gs[idx].p1 - temp;
                                self.bs[idx].p0 += temp * C[0];
                                self.bs[idx].p1 += temp;
                            }
                        },
                        3 => {
                            for idx in 0..n {
                                let a_old = initial_accelerations[idx];
                                let a_new = state.particles_1[idx].acceleration;
                                let mut temp = self.gs[idx].p2;
                                self.gs[idx].p2 = (((a_new - a_old) / RR[3] - self.gs[idx].p0) / RR[4] - self.gs[idx].p1) / RR[5];
                                temp = self.gs[idx].p2 - temp;
                                self.bs[idx].p0 += temp * C[1];
                                self.bs[idx].p1 += temp * C[2];
                                self.bs[idx].p2 += temp;
                            }
                        },
                        4 => {
                            for idx in 0..n {
                                let a_old = initial_accelerations[idx];
                                let a_new = state.particles_1[idx].acceleration;
                                let mut temp = self.gs[idx].p3;
                                self.gs[idx].p3 = ((((a_new - a_old) / RR[6] - self.gs[idx].p0) / RR[7] - self.gs[idx].p1) / RR[8] - self.gs[idx].p2) / RR[9];
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
                                let a_new = state.particles_1[idx].acceleration;
                                let mut temp = self.gs[idx].p4;
                                self.gs[idx].p4 = (((((a_new - a_old) / RR[10] - self.gs[idx].p0) / RR[11] - self.gs[idx].p1) / RR[12] - self.gs[idx].p2) / RR[13] - self.gs[idx].p3) / RR[14];
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
                                let a_new = state.particles_1[idx].acceleration;
                                let mut temp = self.gs[idx].p5;
                                self.gs[idx].p5 = ((((((a_new - a_old) / RR[15] - self.gs[idx].p0) / RR[16] - self.gs[idx].p1) / RR[17] - self.gs[idx].p2) / RR[18] - self.gs[idx].p3) / RR[19] - self.gs[idx].p4) / RR[20];
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
                            let mut max_b6_temp = 0.0;
                            let mut max_acceleration = 0.0;
                            for idx in 0..n {
                                let a_old = initial_accelerations[idx];
                                let a_new = state.particles_1[idx].acceleration;
                                let mut temp = self.gs[idx].p6;
                                self.gs[idx].p6 = (((((((a_new - a_old) / RR[21] - self.gs[idx].p0) / RR[22] - self.gs[idx].p1) / RR[23] - self.gs[idx].p2) / RR[24] - self.gs[idx].p3) / RR[25] - self.gs[idx].p4) / RR[26] - self.gs[idx].p5) / RR[27];
                                temp = self.gs[idx].p6 - temp;
                                self.bs[idx].p0 += temp * C[15];
                                self.bs[idx].p1 += temp * C[16];
                                self.bs[idx].p2 += temp * C[17];
                                self.bs[idx].p3 += temp * C[18];
                                self.bs[idx].p4 += temp * C[19];
                                self.bs[idx].p5 += temp * C[20];
                                self.bs[idx].p6 += temp;

                                if true {
                                    let temp_norm = temp.norm();
                                    if temp_norm > max_b6_temp && temp_norm.is_normal() {
                                        max_b6_temp = temp_norm;
                                    }
                                    let a_new_norm = a_new.norm();
                                    if a_new_norm > max_acceleration && a_new_norm.is_normal() {
                                        max_acceleration = a_new.norm();
                                    }
                                    let error = max_b6_temp / max_acceleration;
                                    if (error.is_normal()) & (error > predictor_corrector_error) {
                                        predictor_corrector_error = error;
                                    }
                                } else {
                                    predictor_corrector_error = temp.norm() / a_new.norm();
                                }
                            }
                        },
                        _ => {}
                    } // end match
                } // end substep loop
            } // end predictor-corrector loop

            // Compute the new timestep.
            let old_timestep = self.timestep;
            let mut new_timestep = calculate_new_timestep(&state.particles_1, &initial_accelerations, &self.bs, &old_timestep, &self.epsilon);
            let timestep_ratio = (new_timestep / old_timestep).abs();

            // If the new timestep is too small, reject the step and try again.
            if timestep_ratio < SAFETY_FACTOR {
                self.timestep = new_timestep;
                for idx in 0..n {
                    state.particles_1[idx].position = initial_positions[idx];
                    state.particles_1[idx].velocity = initial_velocities[idx];
                    state.particles_1[idx].acceleration = initial_accelerations[idx];
                    state.particles_1[idx].epoch = initial_epoch;
                }
                if self.last_timestep != 0.0 {
                    predict_next_coefficients(&timestep_ratio, &self.es_last, &self.bs_last, &mut self.es, &mut self.bs);
                }
                continue 'integration_loop;
            }

            update_accelerations(state, forces);

            // Accept the step: update epoch and particles.
            // *(&mut state.epoch) += self.timestep;
            let new_epoch = initial_epoch + self.timestep;
            for idx in 0..n {
                let b = &self.bs[idx];
                state.particles_1[idx].epoch = new_epoch;
                state.particles_1[idx].position = initial_positions[idx]
                    + self.timestep * initial_velocities[idx]
                    + self.timestep.powi(2)
                        * (initial_accelerations[idx] / 2.0
                            + b.p0 / 6.0
                            + b.p1 / 12.0
                            + b.p2 / 20.0
                            + b.p3 / 30.0
                            + b.p4 / 42.0
                            + b.p5 / 56.0
                            + b.p6 / 72.0);
                state.particles_1[idx].velocity = initial_velocities[idx]
                    + self.timestep
                        * (initial_accelerations[idx]
                            + b.p0 / 2.0
                            + b.p1 / 3.0
                            + b.p2 / 4.0
                            + b.p3 / 5.0
                            + b.p4 / 6.0
                            + b.p5 / 7.0
                            + b.p6 / 8.0);
                // state.particles_1[idx].acceleration = accelerations[idx];

                state.particles_0[idx].position = initial_positions[idx];
                state.particles_0[idx].velocity = initial_velocities[idx];
                state.particles_0[idx].acceleration = initial_accelerations[idx];
                state.particles_0[idx].epoch = initial_epoch;
            }
            self.last_timestep = self.timestep.clone();
            self.timestep = new_timestep;
            let ratio = self.timestep / self.last_timestep;            

            self.es_last = self.es.clone();
            self.bs_last = self.bs.clone();

            predict_next_coefficients(&ratio, &self.es_last, &self.bs_last, &mut self.es, &mut self.bs);
            break 'integration_loop;
        }
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
pub fn calculate_new_timestep(particles: &Vec<SimulationParticle>, accelerations: &Vec<Vector3<f64>>, bs: &Vec<CoefficientSeptet>, last_timestep: &f64, epsilon: &f64) -> f64 {
    let mut min_timescale2 = f64::INFINITY;
    for idx in 0..particles.len() {
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


pub fn update_accelerations(state: &mut SimulationState, forces: &Vec<Box<dyn Force + Send + Sync>>) {
    // first clear the accelerations vector
    for acc in state.particles_1.iter_mut() {
        acc.acceleration = Vector3::zeros();
    }
    for force in forces {
        force.apply_acceleration(state);
    }
}