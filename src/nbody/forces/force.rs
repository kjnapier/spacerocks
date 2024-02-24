use crate::spacerock::SpaceRock;
use nalgebra::Vector3;

pub trait Force: Send + Sync {
    fn apply(&self, particles: &mut Vec<SpaceRock>, perturbers: &mut Vec<SpaceRock>);
    // fn calculate_accelerations(&self, particles: &Vec<SpaceRock>, perturbers: &Vec<SpaceRock>) -> [Vec<Vector3<f64>>; 2];
}