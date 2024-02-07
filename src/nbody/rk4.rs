
// // implementing the Runge-Kutta 4th order method for the n-body problem
// use crate::nbody::acceleration::{calculate_and_set_test_particle_accelerations, calculate_and_set_perturber_accelerations};
// use crate::spacerock::SpaceRock;


// pub fn rk4_step(perturbers: &mut Vec<SpaceRock>, test_particles: &mut Vec<SpaceRock>, timestep: f64) {
//     let mut K1x_test = Vec::new();
//     let mut K2x_test = Vec::new();
//     let mut K3x_test = Vec::new();
//     let mut K4x_test = Vec::new();

//     let mut K1v_test = Vec::new();
//     let mut K2v_test = Vec::new();
//     let mut K3v_test = Vec::new();
//     let mut K4v_test = Vec::new();

//     let mut K1x_perturber = Vec::new();
//     let mut K2x_perturber = Vec::new();
//     let mut K3x_perturber = Vec::new();
//     let mut K4x_perturber = Vec::new();
    
//     let mut K1v_perturber = Vec::new();
//     let mut K2v_perturber = Vec::new();
//     let mut K3v_perturber = Vec::new();
//     let mut K4v_perturber = Vec::new();

//     // get the initial positions of the test particles
//     let mut x0_test = Vec::new();
//     for particle in &*test_particles {
//         x0_test.push(particle.position);
//     }

//     // get the initial positions of the perturbers
//     let mut x0_perturber = Vec::new();
//     for particle in &*perturbers {
//         x0_perturber.push(particle.position);
//     }

//     // calculate the K1 values
//     calculate_and_set_perturber_accelerations(perturbers);
//     calculate_and_set_test_particle_accelerations(perturbers, test_particles);
//     for particle in &mut *test_particles {
//         K1x_test.push(particle.velocity);
//         K1v_test.push(particle.acceleration);
//     }
//     for particle in &mut *perturbers {
//         K1x_perturber.push(particle.velocity);
//         K1v_perturber.push(particle.acceleration);
//     }

//     // calculate the K2 values
//     for (i, particle) in test_particles.iter_mut().enumerate() {
//         particle.position += 0.5 * timestep * K1x_test[i];
//         particle.velocity += 0.5 * timestep * K1v_test[i];
//     }
//     for (i, perturber) in perturbers.iter_mut().enumerate() {
//         perturber.position += 0.5 * timestep * K1x_perturber[i];
//         perturber.velocity += 0.5 * timestep * K1v_perturber[i];
//     }
//     calculate_and_set_perturber_accelerations(perturbers);
//     calculate_and_set_test_particle_accelerations(perturbers, test_particles);
//     for particle in &mut *test_particles {
//         K2x_test.push(particle.velocity);
//         K2v_test.push(particle.acceleration);
//     }
//     for particle in &mut *perturbers {
//         K2x_perturber.push(particle.velocity);
//         K2v_perturber.push(particle.acceleration);
//     }

//     // calculate the K3 values
//     for (i, particle) in test_particles.iter_mut().enumerate() {
//         particle.position += 0.5 * timestep * K2x_test[i];
//         particle.velocity += 0.5 * timestep * K2v_test[i];
//     }
//     for (i, perturber) in perturbers.iter_mut().enumerate() {
//         perturber.position += 0.5 * timestep * K2x_perturber[i];
//         perturber.velocity += 0.5 * timestep * K2v_perturber[i];
//     }
//     calculate_and_set_perturber_accelerations(perturbers);
//     calculate_and_set_test_particle_accelerations(perturbers, test_particles);

//     for particle in &mut *test_particles {
//         K3x_test.push(particle.velocity);
//         K3v_test.push(particle.acceleration);
//     }
//     for particle in &mut *perturbers {
//         K3x_perturber.push(particle.velocity);
//         K3v_perturber.push(particle.acceleration);
//     }

//     // calculate the K4 values
//     for (i, particle) in test_particles.iter_mut().enumerate() {
//         particle.position += timestep * K3x_test[i];
//         particle.velocity += timestep * K3v_test[i];
//     }
//     for (i, perturber) in perturbers.iter_mut().enumerate() {
//         perturber.position += timestep * K3x_perturber[i];
//         perturber.velocity += timestep * K3v_perturber[i];
//     }
//     calculate_and_set_perturber_accelerations(perturbers);
//     calculate_and_set_test_particle_accelerations(perturbers, test_particles);
//     for particle in &mut *test_particles {
//         K4x_test.push(particle.velocity);
//         K4v_test.push(particle.acceleration);
//     }
//     for particle in &mut *perturbers {
//         K4x_perturber.push(particle.velocity);
//         K4v_perturber.push(particle.acceleration);
//     }

//     // // calculate the new positions and velocities
//     // for (i, particle) in test_particles.iter_mut().enumerate() {
//     //     particle.position += timestep * (K1x_test[i] + 2.0 * K2x_test[i] + 2.0 * K3x_test[i] + K4x_test[i]) / 6.0;
//     //     particle.velocity += timestep * (K1v_test[i] + 2.0 * K2v_test[i] + 2.0 * K3v_test[i] + K4v_test[i]) / 6.0;
//     // }
//     // for (i, perturber) in perturbers.iter_mut().enumerate() {
//     //     perturber.position += timestep * (K1x_perturber[i] + 2.0 * K2x_perturber[i] + 2.0 * K3x_perturber[i] + K4x_perturber[i]) / 6.0;
//     //     perturber.velocity += timestep * (K1v_perturber[i] + 2.0 * K2v_perturber[i] + 2.0 * K3v_perturber[i] + K4v_perturber[i]) / 6.0;
//     // }

//     // calculate the new positions and velocities
//     for (i, particle) in test_particles.iter_mut().enumerate() {
//         particle.position = x0_test[i] + timestep * (K1x_test[i] + 2.0 * K2x_test[i] + 2.0 * K3x_test[i] + K4x_test[i]) / 6.0;
//         particle.velocity = K1x_test[i] + timestep * (K1v_test[i] + 2.0 * K2v_test[i] + 2.0 * K3v_test[i] + K4v_test[i]) / 6.0;
//     }

//     for (i, perturber) in perturbers.iter_mut().enumerate() {
//         perturber.position = x0_perturber[i] + timestep * (K1x_perturber[i] + 2.0 * K2x_perturber[i] + 2.0 * K3x_perturber[i] + K4x_perturber[i]) / 6.0;
//         perturber.velocity = K1x_perturber[i] + timestep * (K1v_perturber[i] + 2.0 * K2v_perturber[i] + 2.0 * K3v_perturber[i] + K4v_perturber[i]) / 6.0;
//     }
   
// }