use crate::constants::G;
use crate::spacerock::SpaceRock;
use crate::nbody::simparticle::SimParticle;
use nalgebra::Vector3;



pub fn euler_cromer_step(perturbers: &mut Vec<SimParticle>, test_particles: &mut Vec<SimParticle>, timestep: f64) -> Result<(), Box<dyn std::error::Error>> {
    
    for particle in test_particles {
        for perturber in &mut *perturbers {
            let r_vec = particle.position - perturber.position;
            let r = r_vec.norm();
            let r_hat = r_vec / r;
            let a = -G * perturber.mass * r_hat / r.powi(2);
            particle.acceleration += a;
        }
        particle.position += 0.5 * timestep * particle.velocity;
        particle.velocity += timestep * particle.acceleration;
        particle.position += 0.5 * timestep * particle.velocity;
        
        particle.epoch += timestep;
        particle.acceleration = Vector3::new(0.0, 0.0, 0.0);
    }

    let N = perturbers.len();
    for idx in 0..N {
        for jdx in (idx + 1)..N {
            let r_vec = perturbers[idx].position - perturbers[jdx].position;
            let r = r_vec.norm();
            let xi = -G * r_vec / r.powi(3);
            let da_i = xi * perturbers[jdx].mass;
            let da_j = -xi * perturbers[idx].mass;

            perturbers[idx].acceleration += da_i;
            perturbers[jdx].acceleration += da_j;
        }
    }

    for mut perturber in perturbers.iter_mut() {
        perturber.position += 0.5 * timestep * perturber.velocity;
        perturber.velocity += timestep * perturber.acceleration;
        perturber.position += 0.5 * timestep * perturber.velocity;

        perturber.epoch += timestep;
        perturber.acceleration = Vector3::new(0.0, 0.0, 0.0);
    }

    Ok(())

}