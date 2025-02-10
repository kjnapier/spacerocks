use crate::nbody::integrators::Integrator;
use crate::nbody::forces::Force;
use crate::SpaceRock;
use crate::time::Time;

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
pub struct IAS15 {
    pub timestep: f64,
    pub epsilon: f64,
    pub last_timestep: f64,
    bs: Vec<CoefficientSeptet>,
    gs: Vec<CoefficientSeptet>,
    es: Vec<CoefficientSeptet>,
    bs_last: Vec<CoefficientSeptet>,
    es_last: Vec<CoefficientSeptet>,
}

impl IAS15 {
    pub fn new(timestep: f64) -> IAS15 {
        IAS15 { timestep, epsilon: 1e-9, last_timestep: 0.0, bs: vec![], gs: vec![], es: vec![], bs_last: vec![], es_last: vec![] }
    }

    pub fn reset_coefficients(&mut self, n: usize) {
        self.bs = vec![CoefficientSeptet::zeros(); n];
        self.gs = vec![CoefficientSeptet::zeros(); n];
        self.es = vec![CoefficientSeptet::zeros(); n];
        self.bs_last = vec![CoefficientSeptet::zeros(); n];
        self.es_last = vec![CoefficientSeptet::zeros(); n];
    }
}

impl Integrator for IAS15 {

    fn step(&mut self, particles: &mut Vec<SpaceRock>, epoch: &mut Time, forces: &Vec<Box<dyn Force + Send + Sync>>) {
        // for now I'll only integrate the particles, just to keep things simple

        let mut accelerations: Vec<Vector3<f64>> = vec![Vector3::zeros(); particles.len()];
        for force in forces {
            let acc = force.calculate_acceleration(particles);
            for (idx, a) in acc.iter().enumerate() {
                accelerations[idx] += a;
            }
        }

        // We don't want to clone the original SpaceRock objects because some of the contents are heap allocated, making the clone operation expensive.
        let initial_positions: Vec<Vector3<f64>> = particles.iter().map(|p| p.position).collect();
        let initial_velocities: Vec<Vector3<f64>> = particles.iter().map(|p| p.velocity).collect();
        let initial_accelerations: Vec<Vector3<f64>> = accelerations.clone();
        let initial_epoch = particles[0].epoch.clone();

        // Number of particles
        let n = particles.len();

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
            if iterations >= 10 {
                println!("At least 10 predictor corrector loops in IAS15 did not converge. This is typically an indication of the timestep being too large.");
                self.timestep /= 2.0;
                println!("Reducing the timestep to {}", self.timestep);
                self.step(particles, epoch, forces);
            }

            predictor_corrector_error_last = predictor_corrector_error;
            predictor_corrector_error = 0.0;
            iterations += 1;

            for substep in 1..8 {
                for idx in 0..n {
                    let a0 = initial_accelerations[idx];
                    let v0 = initial_velocities[idx];

                    let b = &self.bs[idx];
                    // let g = &self.gs[idx];
                    let hh = H[substep];

                    // Calculate the position
                    let d_position = ((((((((b.p6 * 7.0 * hh / 9.0 + b.p5) * 3.0 * hh / 4.0 + b.p4) * 5.0 * hh / 7.0 + b.p3) * 2.0 * hh / 3.0 + b.p2) * 3.0 * hh / 5.0 + b.p1) * hh / 2.0 + b.p0) * hh / 3.0 + a0) * self.timestep * hh / 2.0 + v0) * self.timestep * hh;
                    particles[idx].position = initial_positions[idx] + d_position;

                    // Calculate the velocity
                    let d_velocity = (((((((b.p6 * 7.0 * hh / 8.0 + b.p5) * 6.0 * hh / 7.0 + b.p4) * 5.0 * hh / 6.0 + b.p3) * 4.0 * hh / 5.0 + b.p2) * 3.0 * hh / 4.0 + b.p1) * 2.0 * hh / 3.0 + b.p0) * hh / 2.0 + a0) * self.timestep * hh;
                    particles[idx].velocity = initial_velocities[idx] + d_velocity;

                    particles[idx].epoch += self.timestep * hh;
                }

                let mut accelerations: Vec<Vector3<f64>> = vec![Vector3::zeros(); particles.len()];
                for force in forces {
                    let acc = force.calculate_acceleration(particles);
                    for (idx, a) in acc.iter().enumerate() {
                        accelerations[idx] += a;
                    }
                }

                match substep {
                    1 => {
                        for idx in 0..n {
                            let a_old = initial_accelerations[idx];
                            let a_new = accelerations[idx];

                            let temp = self.gs[idx].p0;

                            self.gs[idx].p0 = (a_new - a_old) / RR[0];
                            self.bs[idx].p0 += self.gs[idx].p0 - temp;
                        }
                    },
                    2 => {
                        for idx in 0..n {
                            let a_old = initial_accelerations[idx];
                            let a_new = accelerations[idx];

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
                            let a_new = accelerations[idx];

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
                            let a_new = accelerations[idx];

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
                            let a_new = accelerations[idx];

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
                            let a_new = accelerations[idx];

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

                        let mut max_acceleration = 0.0;
                        let mut max_b6_temp = 0.0;
                        for idx in 0..n {
                            let a_old = initial_accelerations[idx];
                            let a_new = accelerations[idx];

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
                }
            }
        }

        let old_timestep = self.timestep;
        let mut new_timestep = calculate_new_timestep(particles, &initial_accelerations, &self.bs, &old_timestep, &self.epsilon);
        let timestep_ratio = (new_timestep / old_timestep).abs();
        

        // Step was rejected
        if timestep_ratio < SAFETY_FACTOR {
            println!("Timestep was rejected. Reducing the timestep to {}", new_timestep);
            self.timestep = new_timestep;

            // reset particles
            for idx in 0..n {
                let perturber = &mut particles[idx];
                perturber.position = initial_positions[idx];
                perturber.velocity = initial_velocities[idx];
                accelerations[idx] = initial_accelerations[idx];
                // perturber.epoch.epoch = epoch.epoch;
                perturber.epoch = initial_epoch.clone();
            }

            if self.last_timestep != 0.0 {
                // let ratio = self.timestep / self.last_timestep;
                // predict_next_coefficients(&timestep_ratio, &mut self.es, &mut self.bs);
                predict_next_coefficients(&timestep_ratio, &self.es_last, &self.bs_last, &mut self.es, &mut self.bs);
            }

            // recursively call step with the new timestep
            self.step(particles, epoch, forces);
            // return;
        }

        // The timestep was accepted
        if timestep_ratio > 1.0 {
            // println!("Timestep was accepted. Increasing the timestep to {}", new_timestep);
            if timestep_ratio > 1.0 / SAFETY_FACTOR {
                // println!("The timestep ratio is greater than 1/Safety Factor. Increasing the timestep to {}", old_timestep / SAFETY_FACTOR);
                new_timestep = old_timestep / SAFETY_FACTOR;
            }
        }

        

        // Update the epoch
        *epoch += self.timestep; //self.timestep;

        // Update the particles
        for idx in 0..n {
            let b = &self.bs[idx];
            // let g = &self.gs[idx];
            let particle = &mut particles[idx];

            particle.epoch = epoch.clone();
            particle.position = initial_positions[idx] + self.timestep * initial_velocities[idx] + self.timestep.powi(2) * (initial_accelerations[idx] / 2.0 + b.p0 / 6.0 + b.p1 / 12.0 + b.p2 / 20.0 + b.p3 / 30.0 + b.p4 / 42.0 + b.p5 / 56.0 + b.p6 / 72.0);
            particle.velocity = initial_velocities[idx] + self.timestep * (initial_accelerations[idx] + b.p0 / 2.0 + b.p1 / 3.0 + b.p2 / 4.0 + b.p3 / 5.0 + b.p4 / 6.0 + b.p5 / 7.0 + b.p6 / 8.0);
        }

        

        self.last_timestep = self.timestep.clone();
        self.timestep = new_timestep;
        let ratio = self.timestep / self.last_timestep;


        self.es_last = self.es.clone();
        self.bs_last = self.bs.clone();

        predict_next_coefficients(&ratio, &self.es_last, &self.bs_last, &mut self.es, &mut self.bs);        

    }

    fn timestep(&self) -> f64 {
        self.timestep
    }

    fn set_timestep(&mut self, timestep: f64) {
        self.timestep = timestep;
    }
}


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


pub fn calculate_new_timestep(particles: &Vec<SpaceRock>, accelerations: &Vec<Vector3<f64>>, bs: &Vec<CoefficientSeptet>, last_timestep: &f64, epsilon: &f64) -> f64 {
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