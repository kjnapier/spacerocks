// //implement the Runge-Kutta-Fehlberg 4(5) method

// use crate::nbody::acceleration::{calculate_and_set_test_particle_accelerations, calculate_and_set_perturber_accelerations};
// use crate::spacerock::SpaceRock;
// use nalgebra::Vector3;

// pub fn rkf45_step(perturbers: &mut Vec<SpaceRock>, test_particles: &mut Vec<SpaceRock>, timestep: f64) {
//     let mut K1 = Vec::new();
//     let mut K2 = Vec::new();
//     let mut K3 = Vec::new();
//     let mut K4 = Vec::new();
//     let mut K5 = Vec::new();
//     let mut K6 = Vec::new();
//     let mut perturber_accelerations: Vec<Vector3<f64>> = Vec::new();
//     let mut test_particle_accelerations: Vec<Vector3<f64>> = Vec::new();

//     calculate_and_set_perturber_accelerations(perturbers);
//     calculate_and_set_test_particle_accelerations(perturbers, test_particles);
//     for particle in &mut *test_particles {
//         K1.push(particle.acceleration);    
//     }
//     for particle in &mut *perturbers {
//         perturber_accelerations.push(particle.acceleration);
//     }

//     // K2
//     for (i, particle) in test_particles.iter_mut().enumerate() {
//         particle.position += 0.25 * timestep * particle.velocity;
//         particle.velocity += 0.25 * timestep * K1[i];
//     }
//     for (i, perturber) in perturbers.iter_mut().enumerate() {
//         perturber.position += 0.25 * timestep * perturber.velocity;
//         perturber.velocity += 0.25 * timestep * perturber_accelerations[i];
//     }
//     calculate_and_set_perturber_accelerations(perturbers);
//     calculate_and_set_test_particle_accelerations(perturbers, test_particles);
//     for particle in &mut *test_particles {
//         K2.push(particle.acceleration);    
//     }
//     for particle in &mut *perturbers {
//         perturber_accelerations.push(particle.acceleration);
//     }

//     // K3
//     for (i, particle) in test_particles.iter_mut().enumerate() {
//         particle.position += 3.0/32.0 * timestep * particle.velocity;
//         particle.velocity += 3.0/32.0 * timestep * K2[i];
//     }
//     for (i, perturber) in perturbers.iter_mut().enumerate() {
//         perturber.position += 3.0/32.0 * timestep * perturber.velocity;
//         perturber.velocity += 3.0/32.0 * timestep * perturber_accelerations[i];
//     }
//     calculate_and_set_perturber_accelerations(perturbers);
//     calculate_and_set_test_particle_accelerations(perturbers, test_particles);
//     for particle in &mut *test_particles {
//         K3.push(particle.acceleration);    
//     }
//     for particle in &mut *perturbers {
//         perturber_accelerations.push(particle.acceleration);
//     }

//     //K4
//     for (i, particle) in test_particles.iter_mut().enumerate() {
//         particle.position += 1932.0/2197.0 * timestep * particle.velocity;
//         particle.velocity += 1932.0/2197.0 * timestep * K3[i];
//     }
//     for (i, perturber) in perturbers.iter_mut().enumerate() {
//         perturber.position += 1932.0/2197.0 * timestep * perturber.velocity;
//         perturber.velocity += 1932.0/2197.0 * timestep * perturber_accelerations[i];
//     }
//     calculate_and_set_perturber_accelerations(perturbers);
//     calculate_and_set_test_particle_accelerations(perturbers, test_particles);

//     for particle in &mut *test_particles {
//         K4.push(particle.acceleration);    
//     }
//     for particle in &mut *perturbers {
//         perturber_accelerations.push(particle.acceleration);
//     }

//     //K5
//     for (i, particle) in test_particles.iter_mut().enumerate() {
//         particle.position += 439.0/216.0 * timestep * particle.velocity;
//         particle.velocity += 439.0/216.0 * timestep * K4[i];
//     }
//     for (i, perturber) in perturbers.iter_mut().enumerate() {
//         perturber.position += 439.0/216.0 * timestep * perturber.velocity;
//         perturber.velocity += 439.0/216.0 * timestep * perturber_accelerations[i];
//     }
//     calculate_and_set_perturber_accelerations(perturbers);
//     calculate_and_set_test_particle_accelerations(perturbers, test_particles);
//     for particle in &mut *test_particles {
//         K5.push(particle.acceleration);    
//     }
//     for particle in &mut *perturbers {
//         perturber_accelerations.push(particle.acceleration);
//     }

//     //K6
//     for (i, particle) in test_particles.iter_mut().enumerate() {
//         particle.position += -8.0 * timestep * particle.velocity;
//         particle.velocity += -8.0 * timestep * K5[i];
//     }
//     for (i, perturber) in perturbers.iter_mut().enumerate() {
//         perturber.position += -8.0 * timestep * perturber.velocity;
//         perturber.velocity += -8.0 * timestep * perturber_accelerations[i];
//     }
//     calculate_and_set_perturber_accelerations(perturbers);
//     calculate_and_set_test_particle_accelerations(perturbers, test_particles);
//     for particle in &mut *test_particles {
//         K6.push(particle.acceleration);    
//     }
//     for particle in &mut *perturbers {
//         perturber_accelerations.push(particle.acceleration);
//     }

//     //update positions and velocities
//     for (i, particle) in test_particles.iter_mut().enumerate() {
//         particle.position += 25.0/216.0 * timestep * particle.velocity + 1408.0/2565.0 * timestep * K3[i] + 2197.0/4104.0 * timestep * K4[i] - 1.0/5.0 * timestep * K5[i];
//         particle.velocity += 25.0/216.0 * timestep * K1[i] + 1408.0/2565.0 * timestep * K3[i] + 2197.0/4104.0 * timestep * K4[i] - 1.0/5.0 * timestep * K6[i];
//     }

//     for (i, perturber) in perturbers.iter_mut().enumerate() {
//         perturber.position += 25.0/216.0 * timestep * perturber.velocity + 1408.0/2565.0 * timestep * perturber_accelerations[i] + 2197.0/4104.0 * timestep * perturber_accelerations[i] - 1.0/5.0 * timestep * perturber_accelerations[i];
//         perturber.velocity += 25.0/216.0 * timestep * perturber_accelerations[i] + 1408.0/2565.0 * timestep * perturber_accelerations[i] + 2197.0/4104.0 * timestep * perturber_accelerations[i] - 1.0/5.0 * timestep * perturber_accelerations[i];
//     }
// }

