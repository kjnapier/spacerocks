use crate::nbody::integrators::Integrator;
use crate::nbody::forces::Force;
use crate::SpaceRock;
use crate::time::Time;

use rayon::prelude::*;

use nalgebra::Vector3;

// Copilot filled these in. Check that they're correct later
// Gauss Radau spacings
const h: [f64; 8] = [0.0, 0.0562625605369221464656521910318, 0.180240691736892364987579942780, 0.352624717113169637373907769648, 
                     0.547153626330555383001448554766, 0.734210177215410531523210605558, 0.885320946839095768090359771030, 0.977520613561287501891174488626];
// Other constants
const rr: [f64; 28] = [0.0562625605369221464656522, 0.1802406917368923649875799, 0.1239781311999702185219278, 0.3526247171131696373739078, 
                       0.2963621565762474909082556, 0.1723840253762772723863278, 0.5471536263305553830014486, 0.4908910657936332365357964, 
                       0.3669129345936630180138686, 0.1945289092173857456275408, 0.7342101772154105315232106, 0.6779476166784883850575584, 
                       0.5539694854785181665356307, 0.3815854601022408941493028, 0.1870565508848551485217621, 0.8853209468390957680903598, 
                       0.8290583863021736216247076, 0.7050802551022034031027798, 0.5326962297259261307164520, 0.3381673205085403850889112, 
                       0.1511107696236852365671492, 0.9775206135612875018911745, 0.9212580530243653554255223, 0.7972799218243951369035945, 
                       0.6248958964481178645172667, 0.4303669872307321188897259, 0.2433104363458769703679639, 0.0921996667221917338008147];

const c: [f64; 21] = [-0.0562625605369221464656522, 0.0101408028300636299864818, -0.2365032522738145114532321, -0.0035758977292516175949345, 
                       0.0935376952594620658957485, -0.5891279693869841488271399, 0.0019565654099472210769006, -0.0547553868890686864408084, 
                       0.4158812000823068616886219, -1.1362815957175395318285885, -0.0014365302363708915424460, 0.0421585277212687077072973, 
                       -0.3600995965020568122897665, 1.2501507118406910258505441, -1.8704917729329500633517991, 0.0012717903090268677492943, 
                       -0.0387603579159067703699046, 0.3609622434528459832253398, -1.4668842084004269643701553, 2.9061362593084293014237913, 
                       -2.7558127197720458314421588];

const d: [f64; 21] = [0.0562625605369221464656522, 0.0031654757181708292499905, 0.2365032522738145114532321, 0.0001780977692217433881125, 
                      0.0457929855060279188954539, 0.5891279693869841488271399, 0.0000100202365223291272096, 0.0084318571535257015445000, 
                      0.2535340690545692665214616, 1.1362815957175395318285885, 0.0000005637641639318207610, 0.0015297840025004658189490, 
                      0.0978342365324440053653648, 0.8752546646840910912297246, 1.8704917729329500633517991, 0.0000000317188154017613665, 
                      0.0002762930909826476593130, 0.0360285539837364596003871, 0.5767330002770787313544596, 2.2485887607691597933926895, 
                      2.7558127197720458314421588];

const SAFETY_FACTOR: f64 = 0.1;

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
        for force in forces {
            force.apply(particles);
        }

        // We don't want to clone the original SpaceRock objects because some of the contents are heap allocated, making the clone operation expensive.
        let initial_positions: Vec<Vector3<f64>> = particles.iter().map(|p| p.position).collect();
        let initial_velocities: Vec<Vector3<f64>> = particles.iter().map(|p| p.velocity).collect();
        let initial_accelerations: Vec<Vector3<f64>> = particles.iter().map(|p| p.acceleration).collect();

        // Number of particles
        let n = particles.len();

        if (self.bs.len() != n) || (self.gs.len() != n) {
            self.reset_coefficients(n);
        }
      
        for (g, b) in self.gs.iter_mut().zip(self.bs.iter()) {
            g.p0 = b.p6 * d[15] + b.p5 * d[10] + b.p4 * d[6] + b.p3 * d[3] + b.p2 * d[1] + b.p1 * d[0] + b.p0;
            g.p1 = b.p6 * d[16] + b.p5 * d[11] + b.p4 * d[7] + b.p3 * d[4] + b.p2 * d[2] + b.p1;
            g.p2 = b.p6 * d[17] + b.p5 * d[12] + b.p4 * d[8] + b.p3 * d[5] + b.p2;
            g.p3 = b.p6 * d[18] + b.p5 * d[13] + b.p4 * d[9] + b.p3;
            g.p4 = b.p6 * d[19] + b.p5 * d[14] + b.p4;
            g.p5 = b.p6 * d[20] + b.p5;
            g.p6 = b.p6;
        }

        let mut predictor_corrector_error = 1e300;
        let mut predictor_corrector_error_last = 2.0;
        let mut iterations = 0;

        // This is the predictor-corrector loop, which calculates the coefficients for the next step
        while true {

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
                    let g = &self.gs[idx];
                    let hh = h[substep];

                    // Calculate the position
                    let d_position = ((((((((b.p6 * 7.0 * hh / 9.0 + b.p5) * 3.0 * hh / 4.0 + b.p4) * 5.0 * hh / 7.0 + b.p3) * 2.0 * hh / 3.0 + b.p2) * 3.0 * hh / 5.0 + b.p1) * hh / 2.0 + b.p0) * hh / 3.0 + a0) * self.timestep * hh / 2.0 + v0) * self.timestep * hh;
                    particles[idx].position = initial_positions[idx] + d_position;

                    // Calculate the velocity
                    let d_velocity = (((((((b.p6 * 7.0 * hh / 8.0 + b.p5) * 6.0 * hh / 7.0 + b.p4) * 5.0 * hh / 6.0 + b.p3) * 4.0 * hh / 5.0 + b.p2) * 3.0 * hh / 4.0 + b.p1) * 2.0 * hh / 3.0 + b.p0) * hh / 2.0 + a0) * self.timestep * hh;
                    particles[idx].velocity = initial_velocities[idx] + d_velocity;

                    particles[idx].epoch += self.timestep * hh;
                }


                for perturber in particles.iter_mut() {
                    perturber.acceleration = Vector3::zeros();
                }
                for force in forces {
                    force.apply(particles);
                }

                match substep {
                    1 => {
                        for idx in 0..n {
                            let a_old = initial_accelerations[idx];
                            let a_new = particles[idx].acceleration;

                            let temp = self.gs[idx].p0.clone();

                            self.gs[idx].p0 = (a_new - a_old) / rr[0];
                            self.bs[idx].p0 += self.gs[idx].p0 - temp;
                        }
                    },
                    2 => {
                        for idx in 0..n {
                            let a_old = initial_accelerations[idx];
                            let a_new = particles[idx].acceleration;

                            let mut temp = self.gs[idx].p1.clone();
                            self.gs[idx].p1 = ((a_new - a_old) / rr[1] - self.gs[idx].p0) / rr[2];
                            temp = self.gs[idx].p1 - temp;

                            self.bs[idx].p0 += temp * c[0];
                            self.bs[idx].p1 += temp;
                        }
                    },
                    3 => {
                        for idx in 0..n {
                            let a_old = initial_accelerations[idx];
                            let a_new = particles[idx].acceleration;

                            let mut temp = self.gs[idx].p2.clone();
                            self.gs[idx].p2 = (((a_new - a_old) / rr[3] - self.gs[idx].p0) / rr[4] - self.gs[idx].p1) / rr[5];
                            temp = self.gs[idx].p2 - temp;

                            self.bs[idx].p0 += temp * c[1];
                            self.bs[idx].p1 += temp * c[2];
                            self.bs[idx].p2 += temp;
                        }
                    },
                    4 => {
                        for idx in 0..n {
                            let a_old = initial_accelerations[idx];
                            let a_new = particles[idx].acceleration;

                            let mut temp = self.gs[idx].p3.clone();
                            self.gs[idx].p3 = ((((a_new - a_old) / rr[6] - self.gs[idx].p0) / rr[7] - self.gs[idx].p1) / rr[8] - self.gs[idx].p2) / rr[9];
                            temp = self.gs[idx].p3 - temp;

                            self.bs[idx].p0 += temp * c[3];
                            self.bs[idx].p1 += temp * c[4];
                            self.bs[idx].p2 += temp * c[5];
                            self.bs[idx].p3 += temp;
                        }
                    },
                    5 => {
                        for idx in 0..n {
                            let a_old = initial_accelerations[idx];
                            let a_new = particles[idx].acceleration;

                            let mut temp = self.gs[idx].p4.clone();
                            self.gs[idx].p4 = (((((a_new - a_old) / rr[10] - self.gs[idx].p0) / rr[11] - self.gs[idx].p1) / rr[12] - self.gs[idx].p2) / rr[13] - self.gs[idx].p3) / rr[14];
                            temp = self.gs[idx].p4 - temp;

                            self.bs[idx].p0 += temp * c[6];
                            self.bs[idx].p1 += temp * c[7];
                            self.bs[idx].p2 += temp * c[8];
                            self.bs[idx].p3 += temp * c[9];
                            self.bs[idx].p4 += temp;
                        }
                    },
                    6 => {
                        for idx in 0..n {
                            let a_old = initial_accelerations[idx];
                            let a_new = particles[idx].acceleration;

                            let mut temp = self.gs[idx].p5.clone();
                            self.gs[idx].p5 = ((((((a_new - a_old) / rr[15] - self.gs[idx].p0) / rr[16] - self.gs[idx].p1) / rr[17] - self.gs[idx].p2) / rr[18] - self.gs[idx].p3) / rr[19] - self.gs[idx].p4) / rr[20];
                            temp = self.gs[idx].p5 - temp;

                            self.bs[idx].p0 += temp * c[10];
                            self.bs[idx].p1 += temp * c[11];
                            self.bs[idx].p2 += temp * c[12];
                            self.bs[idx].p3 += temp * c[13];
                            self.bs[idx].p4 += temp * c[14];
                            self.bs[idx].p5 += temp;
                        }
                    },
                    7 => {

                        let mut max_acceleration = 0.0;
                        let mut max_b6_temp = 0.0;
                        for idx in 0..n {
                            let a_old = initial_accelerations[idx];
                            let a_new = particles[idx].acceleration;

                            let mut temp = self.gs[idx].p6.clone();
                            self.gs[idx].p6 = (((((((a_new - a_old) / rr[21] - self.gs[idx].p0) / rr[22] - self.gs[idx].p1) / rr[23] - self.gs[idx].p2) / rr[24] - self.gs[idx].p3) / rr[25] - self.gs[idx].p4) / rr[26] - self.gs[idx].p5) / rr[27];
                            temp = self.gs[idx].p6 - temp;

                            self.bs[idx].p0 += temp * c[15];
                            self.bs[idx].p1 += temp * c[16];
                            self.bs[idx].p2 += temp * c[17];
                            self.bs[idx].p3 += temp * c[18];
                            self.bs[idx].p4 += temp * c[19];
                            self.bs[idx].p5 += temp * c[20];
                            self.bs[idx].p6 += temp;

                            if true {
                                let temp_norm = temp.norm();
                                if temp_norm > max_b6_temp {
                                    if temp_norm.is_normal() {
                                        max_b6_temp = temp_norm;
                                    }
                                }
                                let a_new_norm = a_new.norm();
                                if a_new_norm > max_acceleration {
                                    if a_new_norm.is_normal() {
                                        max_acceleration = a_new.norm();
                                    }
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
        let mut new_timestep = calculate_new_timestep(particles, &self.bs, &old_timestep, &self.epsilon);
        let timestep_ratio = (new_timestep / old_timestep).abs();

        if timestep_ratio < SAFETY_FACTOR {
            self.timestep = new_timestep;

            // reset particles
            for idx in 0..n {
                let mut perturber = &mut particles[idx];
                perturber.position = initial_positions[idx];
                perturber.velocity = initial_velocities[idx];
                perturber.acceleration = initial_accelerations[idx];
                perturber.epoch.epoch = epoch.epoch;
            }

            if self.last_timestep != 0.0 {
                let ratio = self.timestep / self.last_timestep;
                // predict_next_coefficients(&timestep_ratio, &mut self.es, &mut self.bs);
                predict_next_coefficients(&timestep_ratio, &self.es_last, &self.bs_last, &mut self.es, &mut self.bs);
            }

            // recursively call step with the new timestep
            self.step(particles, epoch, forces);
        }

        // The timestep was accepted
        if timestep_ratio > 1.0 / SAFETY_FACTOR {
            new_timestep = old_timestep / SAFETY_FACTOR;
        }

        // Update the epoch
        *epoch += self.timestep;

        // Update the perturbers
        for idx in 0..n {
            let b = &self.bs[idx];
            let g = &self.gs[idx];
            let mut particle = &mut particles[idx];

            particle.epoch.epoch = epoch.epoch;
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
    fn new(p0: Vector3<f64>, p1: Vector3<f64>, p2: Vector3<f64>, p3: Vector3<f64>, p4: Vector3<f64>, p5: Vector3<f64>, p6: Vector3<f64>) -> CoefficientSeptet {
        CoefficientSeptet { p0, p1, p2, p3, p4, p5, p6 }
    }

    fn zeros() -> CoefficientSeptet {
        CoefficientSeptet { p0: Vector3::zeros(), p1: Vector3::zeros(), p2: Vector3::zeros(), p3: Vector3::zeros(), p4: Vector3::zeros(), p5: Vector3::zeros(), p6: Vector3::zeros() }
    }
}


fn predict_next_coefficients(ratio: &f64, es_last: &Vec<CoefficientSeptet>, bs_last: &Vec<CoefficientSeptet>, es: &mut Vec<CoefficientSeptet>, bs: &mut Vec<CoefficientSeptet>) {

    let rat = ratio.clone();

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


pub fn calculate_new_timestep(particles: &Vec<SpaceRock>, bs: &Vec<CoefficientSeptet>, last_timestep: &f64, epsilon: &f64) -> f64 {
    let mut min_timescale2 = f64::INFINITY;
    for idx in 0..particles.len() {
        let particle = &particles[idx];
        let b = &bs[idx];
        let a0 = particle.acceleration.norm_squared();
        let y2 = (particle.acceleration + b.p0 + b.p1 + b.p2 + b.p3 + b.p4 + b.p5 + b.p6).norm_squared();
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
        return min_timescale2.sqrt() * last_timestep * (epsilon * 5040.0).powf(1.0 / 7.0);
    } else {
        return last_timestep / SAFETY_FACTOR;
    }

}