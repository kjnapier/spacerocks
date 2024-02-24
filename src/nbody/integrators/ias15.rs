use crate::nbody::integrators::Integrator;
use crate::nbody::forces::Force;
use crate::SpaceRock;

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



pub struct IAS15 {
    pub timestep: f64,
    pub tolerance: f64,
    pub last_step: f64,
    pub bs: Vec<CoefficientSeptet>,
    pub gs: Vec<CoefficientSeptet>,
    adaptive_timestep: bool,
}

impl IAS15 {
    pub fn new(timestep: f64, tolerance: f64) -> IAS15 {
        IAS15 { timestep, tolerance, last_step: 0.0, bs: vec![], gs: vec![], adaptive_timestep: true }
    }

    pub fn reset_coefficients(&mut self, n: usize) {
        self.bs = vec![CoefficientSeptet::zeros(); n];
        self.gs = vec![CoefficientSeptet::zeros(); n];
    }

}

impl Integrator for IAS15 {

    fn step(&mut self, particles: &mut Vec<SpaceRock>, perturbers: &mut Vec<SpaceRock>, forces: &Vec<Box<dyn Force + Send + Sync>>) {
        // for now I'll only integrate the perturbers, just to keep things simple
        for force in forces {
            force.apply(particles, perturbers);
        }


        // Clone the perturbers so we can reset to the initial state if we need to.
        let mut perturbers_clone = perturbers.clone();

        // Number of perturbers
        let n = perturbers.len();

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

        let t_beginning = perturbers[0].epoch.clone();
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
            if iterations >= 12 {
                println!("At least 10 predictor corrector loops in IAS15 did not converge. This is typically an indication of the timestep being too large.");
                break;
            }

            predictor_corrector_error_last = predictor_corrector_error;
            predictor_corrector_error = 0.0;
            iterations += 1;

            for substep in 1..8 {
                for idx in 0..n {
                    let a0 = perturbers[idx].acceleration;
                    let v0 = perturbers[idx].velocity;
                    let b = &self.bs[idx];
                    let g = &self.gs[idx];
                    let hh = h[substep];

                    // Calculate the position
                    let d_position = ((((((((b.p6 * 7.0 * hh / 9.0 + b.p5) * 3.0 * hh / 4.0 + b.p4) * 5.0 * hh / 7.0 + b.p3) * 2.0 * hh / 3.0 + b.p2) * 3.0 * hh / 5.0 + b.p1) * hh / 2.0 + b.p0) * hh / 3.0 + a0) * self.timestep * hh / 2.0 + v0) * self.timestep * hh;
                    perturbers_clone[idx].position = perturbers[idx].position + d_position;

                    // Calculate the velocity
                    // let d_velocity = (((((((b.p6 * 7.0 * hh / 8.0 + b.p5) * 6.0 * hh / 7.0 + b.p4) * 5.0 * hh / 6.0 + b.p3) * 4.0 * hh / 5.0 + b.p2) * 3.0 * hh / 4.0 + b.p1) * 2.0 * hh / 3.0 + b.p0) * hh / 2.0 + a0) * self.timestep * hh;
                    // perturbers_clone[idx].velocity = perturbers[idx].velocity + d_velocity;

                    perturbers_clone[idx].epoch = perturbers[idx].epoch.clone() + self.timestep * hh;
                }

                // Calculate the accelerations
                // zero out the accelerations
                for perturber in perturbers_clone.iter_mut() {
                    perturber.acceleration = Vector3::zeros();
                }
                for particle in particles.iter_mut() {
                    particle.acceleration = Vector3::zeros();
                }
                for force in forces {
                    force.apply(particles, &mut perturbers_clone);
                }

                match substep {
                    1 => {
                        for idx in 0..n {
                            let a_old = perturbers[idx].acceleration;
                            let a_new = perturbers_clone[idx].acceleration;

                            let temp = self.gs[idx].p0.clone();

                            self.gs[idx].p0 = (a_new - a_old) / rr[0];
                            self.bs[idx].p0 += self.gs[idx].p0 - temp;
                        }
                    },
                    2 => {
                        for idx in 0..n {
                            let a_old = perturbers[idx].acceleration;
                            let a_new = perturbers_clone[idx].acceleration;

                            let mut temp = self.gs[idx].p1.clone();
                            self.gs[idx].p1 = ((a_new - a_old) / rr[1] - self.gs[idx].p0) / rr[2];
                            temp = self.gs[idx].p1 - temp;

                            self.bs[idx].p0 += temp * c[0];
                            self.bs[idx].p1 += temp;
                        }
                    },
                    3 => {
                        for idx in 0..n {
                            let a_old = perturbers[idx].acceleration;
                            let a_new = perturbers_clone[idx].acceleration;

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
                            let a_old = perturbers[idx].acceleration;
                            let a_new = perturbers_clone[idx].acceleration;

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
                            let a_old = perturbers[idx].acceleration;
                            let a_new = perturbers_clone[idx].acceleration;

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
                            let a_old = perturbers[idx].acceleration;
                            let a_new = perturbers_clone[idx].acceleration;

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
                            let a_old = perturbers[idx].acceleration;
                            let a_new = perturbers_clone[idx].acceleration;

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

                            if self.adaptive_timestep {
                                if temp.norm() > max_b6_temp {
                                    max_b6_temp = temp.norm();
                                }
                                if a_new.norm() > max_acceleration {
                                    max_acceleration = a_new.norm();
                                }
                                let error = max_b6_temp / max_acceleration;
                                if error > predictor_corrector_error {
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

        // Estimate the error
        
        // Update the perturbers
        for idx in 0..n {
            let b = &self.bs[idx];
            let g = &self.gs[idx];
            let mut perturber = &mut perturbers[idx];

            perturber.position += self.timestep * perturber.velocity + self.timestep.powi(2) * (perturber.acceleration / 2.0 + b.p0 / 6.0 + b.p1 / 12.0 + b.p2 / 20.0 + b.p3 / 30.0 + b.p4 / 42.0 + b.p5 / 56.0 + b.p6 / 72.0);
            perturber.velocity += self.timestep * (perturber.acceleration + b.p0 / 2.0 + b.p1 / 3.0 + b.p2 / 4.0 + b.p3 / 5.0 + b.p4 / 6.0 + b.p5 / 7.0 + b.p6 / 8.0);
            perturber.epoch += self.timestep;
        }

        println!("Predictor-corrector loop iterations: {}", iterations);

    }

    fn timestep(&self) -> f64 {
        self.timestep
    }

    fn set_timestep(&mut self, timestep: f64) {
        self.timestep = timestep;
    }
}

#[derive(Clone)]
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