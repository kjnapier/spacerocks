use crate::constants::MU_BARY;
use crate::StateVector;
use crate::KeplerOrbit;
use nalgebra::Vector3;

const EMIN: f64 = 1e-10;
const IMIN: f64 = 1e-10;

pub fn calc_kep_from_xyz(position: Vector3<f64>, velocity: Vector3<f64>, mu: f64) -> KeplerOrbit {

    // energy
    // h
    // hz
    // 


    let r = position.norm();
    let vsq = velocity.dot(&velocity);
    let hvec = position.cross(&velocity);
    let hsq = hvec.dot(&hvec);
    let rdot = position.dot(&velocity) / r;

    let specific_energy = vsq / 2.0 - mu / r;
 
    let cos_inc = hvec.z / hvec.norm();
    

}