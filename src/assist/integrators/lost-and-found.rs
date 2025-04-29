
// impl Integrator for IAS15 {
//     fn step(&mut self, state: &mut SimulationState, forces: &Vec<Box<dyn Force + Send + Sync>>, kernel: &SpiceKernel) {
//         // Label the outer integration loop to avoid recursion.
//         'integration_loop: loop {
//             let n = state.particles.len();

//             // Check if the number of particles has changed.
//             // if state.spice_particles[0].epoch != state.particles_1[0].epoch {
//             //     let state_vectors = kernel.get_barycentric_states(&state.spice_bodies, state.particles_1[0].epoch).unwrap();
//             //     for (i, state_vector) in state_vectors.iter().enumerate() {
//             //         let (x, y, z, vx, vy, vz) = state_vector;
//             //         let position = Vector3::new(*x, *y, *z);
//             //         let velocity = Vector3::new(*vx, *vy, *vz);
//             //         state.spice_particles[i].position = position;
//             //         state.spice_particles[i].velocity = velocity;
//             //         state.spice_particles[i].epoch = state.particles_1[0].epoch;
//             //     }
//             // }
            
//             let state_vectors = kernel.get_barycentric_states(&state.spice_bodies, state.particles_1[0].epoch).unwrap();
//             for (i, state_vector) in state_vectors.iter().enumerate() {
//                 let (x, y, z, vx, vy, vz) = state_vector;
//                 let position = Vector3::new(*x, *y, *z);
//                 let velocity = Vector3::new(*vx, *vy, *vz);
//                 state.spice_particles[i].position = position;
//                 state.spice_particles[i].velocity = velocity;
//             }

//             // Allocate and fill temporary acceleration vector once.
//             let mut accelerations = vec![Vector3::zeros(); n];
//             for force in forces {
//                 let acc = force.calculate_acceleration(state);
//                 for (i, a) in acc.iter().enumerate() {
//                     accelerations[i] += *a;
//                 }
//             }

//             // Save initial conditions.
//             let initial_positions: Vec<Vector3<f64>> = state.particles_1.iter().map(|p| p.position).collect();
//             let initial_velocities: Vec<Vector3<f64>> = state.particles_1.iter().map(|p| p.velocity).collect();
//             let initial_accelerations: Vec<Vector3<f64>> = accelerations.clone();
//             let initial_epoch = state.particles_1[0].epoch;

//             // Ensure coefficient vectors are sized correctly.
//             if self.bs.len() != n || self.gs.len() != n {
//                 self.reset_coefficients(n);
//             }

//             // Update intermediate coefficients.
//             for (g, b) in self.gs.iter_mut().zip(&self.bs) {
//                 g.p0 = b.p6 * D[15] + b.p5 * D[10] + b.p4 * D[6] + b.p3 * D[3] + b.p2 * D[1] + b.p1 * D[0] + b.p0;
//                 g.p1 = b.p6 * D[16] + b.p5 * D[11] + b.p4 * D[7] + b.p3 * D[4] + b.p2 * D[2] + b.p1;
//                 g.p2 = b.p6 * D[17] + b.p5 * D[12] + b.p4 * D[8] + b.p3 * D[5] + b.p2;
//                 g.p3 = b.p6 * D[18] + b.p5 * D[13] + b.p4 * D[9] + b.p3;
//                 g.p4 = b.p6 * D[19] + b.p5 * D[14] + b.p4;
//                 g.p5 = b.p6 * D[20] + b.p5;
//                 g.p6 = b.p6;
//             }

//             let mut predictor_corrector_error = 1e300;
//             let mut predictor_corrector_error_last = 2.0;
//             let mut iterations = 0;

          
//             let mut substep_map: Vec<Vec<(Vector3<f64>, Vector3<f64>)>> = vec![Vec::new(); 8];
//             for substep in 1..8 {
//                 let hh = H[substep];
//                 let state_vectors = kernel.get_barycentric_states(&state.spice_bodies, initial_epoch + self.timestep * hh).unwrap();
//                 for (i, state_vector) in state_vectors.iter().enumerate() {
//                     let (x, y, z, vx, vy, vz) = state_vector;
//                     let position = Vector3::new(*x, *y, *z);
//                     let velocity = Vector3::new(*vx, *vy, *vz);
//                     substep_map[substep].push((position, velocity));
//                 }
//             }           

//             // Predictor-corrector iteration loop.
//             loop {

//                 // if predictor_corrector_error < 1e-16 {
//                 //     break;
//                 // }
//                 if predictor_corrector_error < 2e-16 {
//                     break;
//                 }

//                 if (iterations > 2 && predictor_corrector_error_last <= predictor_corrector_error) {
//                     break;
//                 }

//                 if iterations >= 10 {
//                     println!("Convergence failed after 10 iterations. Reducing timestep.");
//                     self.timestep /= 2.0;
//                     continue 'integration_loop;
//                 }
//                 predictor_corrector_error_last = predictor_corrector_error;
//                 predictor_corrector_error = 0.0;
//                 iterations += 1;

//                 // Loop over substeps 1..7.
//                 for substep in 1..8 {

//                     let hh = H[substep];
//                     for idx in 0..n {
//                         let a0 = initial_accelerations[idx];
//                         let v0 = initial_velocities[idx];
//                         let b = &self.bs[idx];

//                         // Compute position increment.
//                         let d_position = ((((((((b.p6 * 7.0 * hh / 9.0 + b.p5) * 3.0 * hh / 4.0 + b.p4)
//                             * 5.0 * hh / 7.0 + b.p3) * 2.0 * hh / 3.0 + b.p2)
//                             * 3.0 * hh / 5.0 + b.p1) * hh / 2.0 + b.p0)
//                             * hh / 3.0 + a0) * self.timestep * hh / 2.0 + v0)
//                             * self.timestep * hh;
//                         state.particles_1[idx].position = initial_positions[idx] + d_position;

//                         // Compute velocity increment.
//                         let d_velocity = (((((((b.p6 * 7.0 * hh / 8.0 + b.p5) * 6.0 * hh / 7.0 + b.p4)
//                             * 5.0 * hh / 6.0 + b.p3) * 4.0 * hh / 5.0 + b.p2)
//                             * 3.0 * hh / 4.0 + b.p1) * 2.0 * hh / 3.0 + b.p0)
//                             * hh / 2.0 + a0) * self.timestep * hh;
//                         state.particles_1[idx].velocity = initial_velocities[idx] + d_velocity;

//                         state.particles_1[idx].epoch = initial_epoch + self.timestep * hh;

//                     }

//                     // // Update the state of the simulation with the new positions and velocities.
//                     for (idx, particle) in state.spice_particles.iter_mut().enumerate() {
//                         let sm = &substep_map[substep];
//                         let (position, velocity) = sm[idx];
//                         particle.position = position;
//                         particle.velocity = velocity;
//                         particle.epoch = state.particles_1[0].epoch;
//                     }

//                     // Reset the temporary accelerations vector to reuse it.
//                     for acc in accelerations.iter_mut() {
//                         *acc = Vector3::zeros();
//                     }
//                     for force in forces {
//                         let acc = force.calculate_acceleration(state);
//                         for (i, a) in acc.iter().enumerate() {
//                             accelerations[i] += *a;
//                         }
//                     }

//                     // Update coefficients based on substep.
//                     match substep {
//                         1 => {
//                             for idx in 0..n {
//                                 let a_old = initial_accelerations[idx];
//                                 let a_new = accelerations[idx];
//                                 let temp = self.gs[idx].p0;
//                                 self.gs[idx].p0 = (a_new - a_old) / RR[0];
//                                 self.bs[idx].p0 += self.gs[idx].p0 - temp;
//                             }
//                         },
//                         2 => {
//                             for idx in 0..n {
//                                 let a_old = initial_accelerations[idx];
//                                 let a_new = accelerations[idx];
//                                 let mut temp = self.gs[idx].p1;
//                                 self.gs[idx].p1 = ((a_new - a_old) / RR[1] - self.gs[idx].p0) / RR[2];
//                                 temp = self.gs[idx].p1 - temp;
//                                 self.bs[idx].p0 += temp * C[0];
//                                 self.bs[idx].p1 += temp;
//                             }
//                         },
//                         3 => {
//                             for idx in 0..n {
//                                 let a_old = initial_accelerations[idx];
//                                 let a_new = accelerations[idx];
//                                 let mut temp = self.gs[idx].p2;
//                                 self.gs[idx].p2 = (((a_new - a_old) / RR[3] - self.gs[idx].p0) / RR[4] - self.gs[idx].p1) / RR[5];
//                                 temp = self.gs[idx].p2 - temp;
//                                 self.bs[idx].p0 += temp * C[1];
//                                 self.bs[idx].p1 += temp * C[2];
//                                 self.bs[idx].p2 += temp;
//                             }
//                         },
//                         4 => {
//                             for idx in 0..n {
//                                 let a_old = initial_accelerations[idx];
//                                 let a_new = accelerations[idx];
//                                 let mut temp = self.gs[idx].p3;
//                                 self.gs[idx].p3 = ((((a_new - a_old) / RR[6] - self.gs[idx].p0) / RR[7] - self.gs[idx].p1) / RR[8] - self.gs[idx].p2) / RR[9];
//                                 temp = self.gs[idx].p3 - temp;
//                                 self.bs[idx].p0 += temp * C[3];
//                                 self.bs[idx].p1 += temp * C[4];
//                                 self.bs[idx].p2 += temp * C[5];
//                                 self.bs[idx].p3 += temp;
//                             }
//                         },
//                         5 => {
//                             for idx in 0..n {
//                                 let a_old = initial_accelerations[idx];
//                                 let a_new = accelerations[idx];
//                                 let mut temp = self.gs[idx].p4;
//                                 self.gs[idx].p4 = (((((a_new - a_old) / RR[10] - self.gs[idx].p0) / RR[11] - self.gs[idx].p1) / RR[12] - self.gs[idx].p2) / RR[13] - self.gs[idx].p3) / RR[14];
//                                 temp = self.gs[idx].p4 - temp;
//                                 self.bs[idx].p0 += temp * C[6];
//                                 self.bs[idx].p1 += temp * C[7];
//                                 self.bs[idx].p2 += temp * C[8];
//                                 self.bs[idx].p3 += temp * C[9];
//                                 self.bs[idx].p4 += temp;
//                             }
//                         },
//                         6 => {
//                             for idx in 0..n {
//                                 let a_old = initial_accelerations[idx];
//                                 let a_new = accelerations[idx];
//                                 let mut temp = self.gs[idx].p5;
//                                 self.gs[idx].p5 = ((((((a_new - a_old) / RR[15] - self.gs[idx].p0) / RR[16] - self.gs[idx].p1) / RR[17] - self.gs[idx].p2) / RR[18] - self.gs[idx].p3) / RR[19] - self.gs[idx].p4) / RR[20];
//                                 temp = self.gs[idx].p5 - temp;
//                                 self.bs[idx].p0 += temp * C[10];
//                                 self.bs[idx].p1 += temp * C[11];
//                                 self.bs[idx].p2 += temp * C[12];
//                                 self.bs[idx].p3 += temp * C[13];
//                                 self.bs[idx].p4 += temp * C[14];
//                                 self.bs[idx].p5 += temp;
//                             }
//                         },
//                         7 => {
//                             let mut max_b6_temp = 0.0;
//                             let mut max_acceleration = 0.0;
//                             for idx in 0..n {
//                                 let a_old = initial_accelerations[idx];
//                                 let a_new = accelerations[idx];
//                                 let mut temp = self.gs[idx].p6;
//                                 self.gs[idx].p6 = (((((((a_new - a_old) / RR[21] - self.gs[idx].p0) / RR[22] - self.gs[idx].p1) / RR[23] - self.gs[idx].p2) / RR[24] - self.gs[idx].p3) / RR[25] - self.gs[idx].p4) / RR[26] - self.gs[idx].p5) / RR[27];
//                                 temp = self.gs[idx].p6 - temp;
//                                 self.bs[idx].p0 += temp * C[15];
//                                 self.bs[idx].p1 += temp * C[16];
//                                 self.bs[idx].p2 += temp * C[17];
//                                 self.bs[idx].p3 += temp * C[18];
//                                 self.bs[idx].p4 += temp * C[19];
//                                 self.bs[idx].p5 += temp * C[20];
//                                 self.bs[idx].p6 += temp;

//                                 if true {
//                                     let temp_norm = temp.norm();
//                                     if temp_norm > max_b6_temp && temp_norm.is_normal() {
//                                         max_b6_temp = temp_norm;
//                                     }
//                                     let a_new_norm = a_new.norm();
//                                     if a_new_norm > max_acceleration && a_new_norm.is_normal() {
//                                         max_acceleration = a_new.norm();
//                                     }
//                                     let error = max_b6_temp / max_acceleration;
//                                     if (error.is_normal()) & (error > predictor_corrector_error) {
//                                         predictor_corrector_error = error;
//                                     }
//                                 } else {
//                                     predictor_corrector_error = temp.norm() / a_new.norm();
//                                 }
//                             }
//                         },
//                         _ => {}
//                     } // end match
//                 } // end substep loop
//             } // end predictor-corrector loop

//             // Compute the new timestep.
//             let old_timestep = self.timestep;
//             let mut new_timestep = calculate_new_timestep(&state.particles_1, &initial_accelerations, &self.bs, &old_timestep, &self.epsilon);
//             let timestep_ratio = (new_timestep / old_timestep).abs();

//             // If the new timestep is too small, reject the step and try again.
//             if timestep_ratio < SAFETY_FACTOR {
//                 self.timestep = new_timestep;
//                 for idx in 0..n {
//                     state.particles_1[idx].position = initial_positions[idx];
//                     state.particles_1[idx].velocity = initial_velocities[idx];
//                     state.particles_1[idx].acceleration = initial_accelerations[idx];
//                     // accelerations[idx] = initial_accelerations[idx];
//                     state.particles_1[idx].epoch = initial_epoch;
//                 }
//                 if self.last_timestep != 0.0 {
//                     predict_next_coefficients(&timestep_ratio, &self.es_last, &self.bs_last, &mut self.es, &mut self.bs);
//                 }
//                 continue 'integration_loop;
//             }

//             // calculate the acceleration at the new timestep
//             // clear the accelerations vector
//             for acc in accelerations.iter_mut() {
//                 *acc = Vector3::zeros();
//             }
//             for force in forces {
//                 let acc = force.calculate_acceleration(state);
//                 for (i, a) in acc.iter().enumerate() {
//                     accelerations[i] += *a;
//                 }
//             }

//             // Accept the step: update epoch and particles.
//             // *(&mut state.epoch) += self.timestep;
//             let new_epoch = initial_epoch + self.timestep;
//             for idx in 0..n {
//                 let b = &self.bs[idx];
//                 state.particles_1[idx].epoch = new_epoch;
//                 state.particles_1[idx].position = initial_positions[idx]
//                     + self.timestep * initial_velocities[idx]
//                     + self.timestep.powi(2)
//                         * (initial_accelerations[idx] / 2.0
//                             + b.p0 / 6.0
//                             + b.p1 / 12.0
//                             + b.p2 / 20.0
//                             + b.p3 / 30.0
//                             + b.p4 / 42.0
//                             + b.p5 / 56.0
//                             + b.p6 / 72.0);
//                 state.particles_1[idx].velocity = initial_velocities[idx]
//                     + self.timestep
//                         * (initial_accelerations[idx]
//                             + b.p0 / 2.0
//                             + b.p1 / 3.0
//                             + b.p2 / 4.0
//                             + b.p3 / 5.0
//                             + b.p4 / 6.0
//                             + b.p5 / 7.0
//                             + b.p6 / 8.0);
//                 state.particles_1[idx].acceleration = accelerations[idx];

//                 state.particles_0[idx].position = initial_positions[idx];
//                 state.particles_0[idx].velocity = initial_velocities[idx];
//                 state.particles_0[idx].acceleration = initial_accelerations[idx];
//                 state.particles_0[idx].epoch = initial_epoch;
//             }
//             self.last_timestep = self.timestep.clone();
//             self.timestep = new_timestep;
//             let ratio = self.timestep / self.last_timestep;            

//             self.es_last = self.es.clone();
//             self.bs_last = self.bs.clone();

//             predict_next_coefficients(&ratio, &self.es_last, &self.bs_last, &mut self.es, &mut self.bs);
//             break 'integration_loop;
//         }
//     }

//     fn timestep(&self) -> f64 {
//         self.timestep
//     }

//     fn set_timestep(&mut self, timestep: f64) {
//         self.timestep = timestep;
//     }

//     fn last_timestep(&self) -> f64 {
//         self.last_timestep
//     }

//     fn bs_last(&self) -> &Vec<CoefficientSeptet> {
//         &self.bs_last
//     }
// }