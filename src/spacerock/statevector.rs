use nalgebra::Vector3;

#[derive(Debug, Clone, PartialEq)]
pub struct StateVector {
    pub position: Vector3<f64>,
    pub velocity: Vector3<f64>,
}

impl StateVector {
    pub fn new(x: f64, y: f64, z: f64, vx: f64, vy: f64, vz: f64) -> Self {
        StateVector {
            position: Vector3::new(x, y, z),
            velocity: Vector3::new(vx, vy, vz),
        }
    }

    pub fn from_vectors(position: Vector3<f64>, velocity: Vector3<f64>) -> Self {
        StateVector {
            position,
            velocity,
        }
    }
}