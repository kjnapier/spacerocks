use crate::SpaceRock;
use crate::Simulation;
use crate::nbody::forces::Force;

use core::fmt::Formatter;
use core::fmt;
use core::fmt::Debug;

pub trait Integrator: Send + Sync {
    fn step(&mut self, particles: &mut Vec<SpaceRock>, perturbers: &mut Vec<SpaceRock>, forces: &Vec<Box<dyn Force + Send + Sync>>);
    fn timestep(&self) -> f64;
    fn set_timestep(&mut self, timestep: f64);
}