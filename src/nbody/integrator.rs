use crate::nbody::leapfrog::Leapfrog;
use crate::nbody::freezer::Freezer;
use crate::nbody::rk4::RK4;
use crate::SpaceRock;

#[derive(PartialEq, Debug, Clone)]  
pub enum Integrator {
    Leapfrog(Leapfrog),
    Freezer(Freezer),
    RK4(RK4),
}

impl Integrator {
    pub fn step(&self, perturbers: &mut Vec<SpaceRock>, particles: &mut Vec<SpaceRock>) {
        match self {
            Integrator::Leapfrog(integrator) => integrator.step(perturbers, particles),
            Integrator::Freezer(integrator) => integrator.step(perturbers, particles),
            Integrator::RK4(integrator) => integrator.step(perturbers, particles),
        }
    }

    pub fn threaded_step(&self, perturbers: &mut Vec<SpaceRock>, particles: &mut Vec<SpaceRock>) {
        match self {
            Integrator::Leapfrog(integrator) => integrator.threaded_step(perturbers, particles),
            Integrator::Freezer(integrator) => integrator.step(perturbers, particles),
            Integrator::RK4(integrator) => integrator.step(perturbers, particles),
        }
    }

    pub fn timestep(&self) -> f64 {
        match self {
            Integrator::Leapfrog(integrator) => integrator.timestep,
            Integrator::Freezer(integrator) => integrator.timestep,
            Integrator::RK4(integrator) => integrator.timestep,
        }
    }

    pub fn set_timestep(&mut self, timestep: f64) {
        match self {
            Integrator::Leapfrog(integrator) => integrator.timestep = timestep,
            Integrator::Freezer(integrator) => integrator.timestep = timestep,
            Integrator::RK4(integrator) => integrator.timestep = timestep,
        }
    }
}

