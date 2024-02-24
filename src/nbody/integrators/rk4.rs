
use crate::SpaceRock;
use crate::Simulation;
use crate::nbody::integrators::Integrator;
use crate::nbody::forces::Force;

use rayon::prelude::*;

#[derive(PartialEq, Debug, Clone, Copy)]
pub struct RK4 {
    pub timestep: f64,
}

impl RK4 {
    pub fn new(timestep: f64) -> RK4 {
        RK4 { timestep }
    }
}

impl Integrator for RK4 {

    fn step(&mut self, particles: &mut Vec<SpaceRock>, perturbers: &mut Vec<SpaceRock>, forces: &Vec<Box<dyn Force + Send + Sync>>) {

        let mut k1 = particles.clone();
        let mut l1 = perturbers.clone();
        for force in forces {
            force.apply(&mut k1, &mut l1);
        }

        let mut k2 = particles.clone();
        let mut l2 = perturbers.clone();
        for idx in 0..particles.len() {
            k2[idx].position += 0.5 * self.timestep * k1[idx].velocity;
            k2[idx].velocity += 0.5 * self.timestep * k1[idx].acceleration;
        }
        for idx in 0..perturbers.len() {
            l2[idx].position += 0.5 * self.timestep * l1[idx].velocity;
            l2[idx].velocity += 0.5 * self.timestep * l1[idx].acceleration;
        }
        for force in forces {
            force.apply(&mut k2, &mut l2);
        }

        let mut k3 = particles.clone();
        let mut l3 = perturbers.clone();
        for idx in 0..particles.len() {
            k3[idx].position += 0.5 * self.timestep * k2[idx].velocity;
            k3[idx].velocity += 0.5 * self.timestep * k2[idx].acceleration;
        }
        for idx in 0..perturbers.len() {
            l3[idx].position += 0.5 * self.timestep * l2[idx].velocity;
            l3[idx].velocity += 0.5 * self.timestep * l2[idx].acceleration;
        }
        for force in forces {
            force.apply(&mut k3, &mut l3);
        }

        let mut k4 = particles.clone();
        let mut l4 = perturbers.clone();
        for idx in 0..particles.len() {
            k4[idx].position += self.timestep * k3[idx].velocity;
            k4[idx].velocity += self.timestep * k3[idx].acceleration;
        }
        for idx in 0..perturbers.len() {
            l4[idx].position += self.timestep * l3[idx].velocity;
            l4[idx].velocity += self.timestep * l3[idx].acceleration;
        }
        for force in forces {
            force.apply(&mut k4, &mut l4);
        }

        // update positions and velocities
        for idx in 0..particles.len() {
            particles[idx].position += self.timestep / 6.0 * (k1[idx].velocity + 2.0 * k2[idx].velocity + 2.0 * k3[idx].velocity + k4[idx].velocity);
            particles[idx].velocity += self.timestep / 6.0 * (k1[idx].acceleration + 2.0 * k2[idx].acceleration + 2.0 * k3[idx].acceleration + k4[idx].acceleration);
            particles[idx].epoch += self.timestep;
        }
        for idx in 0..perturbers.len() {
            perturbers[idx].position += self.timestep / 6.0 * (l1[idx].velocity + 2.0 * l2[idx].velocity + 2.0 * l3[idx].velocity + l4[idx].velocity);
            perturbers[idx].velocity += self.timestep / 6.0 * (l1[idx].acceleration + 2.0 * l2[idx].acceleration + 2.0 * l3[idx].acceleration + l4[idx].acceleration);
            perturbers[idx].epoch += self.timestep;
        }


    }

    fn timestep(&self) -> f64 {
        self.timestep
    }

    fn set_timestep(&mut self, timestep: f64) {
        self.timestep = timestep;
    }
}