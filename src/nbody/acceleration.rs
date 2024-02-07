use crate::constants::G;
use crate::nbody::simparticle::SimParticle;
use nalgebra::Vector3;

fn calculate_acceleration_test_particle(perturbers: &mut Vec<SimParticle>, test_particle: &mut SimParticle) {
    test_particle.acceleration = Vector3::new(0.0, 0.0, 0.0);
    for perturber in &mut *perturbers {
        let r_vec = test_particle.position - perturber.position;
        let r = r_vec.norm();
        let r_hat = r_vec / r;
        let a = -G * perturber.mass * r_hat / r.powi(2);
        test_particle.acceleration += a;
    }
}

pub fn calculate_and_set_test_particle_accelerations(perturbers: &mut Vec<SimParticle>, test_particles: &mut Vec<SimParticle>) {
    for test_particle in test_particles.iter_mut() {
        calculate_acceleration_test_particle(perturbers, test_particle);
    }
}

pub fn calculate_and_set_perturber_accelerations(perturbers: &mut Vec<SimParticle>) {
    for perturber in perturbers.iter_mut() {
        perturber.acceleration = Vector3::new(0.0, 0.0, 0.0);
    }
    let N = perturbers.len();
    for idx in 0..N {
        for jdx in (idx + 1)..N {
            let r_vec = perturbers[idx].position - perturbers[jdx].position;
            let r = r_vec.norm();
            let xi = -G * r_vec / (r * r * r);
            let da_i = xi * perturbers[jdx].mass;
            let da_j = -xi * perturbers[idx].mass;

            perturbers[idx].acceleration += da_i;
            perturbers[jdx].acceleration += da_j;
        }
    }
}