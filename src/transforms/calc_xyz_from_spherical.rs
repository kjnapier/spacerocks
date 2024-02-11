// use crate::statevector::StateVector;
// use crate::sphericalstate::SphericalState;
// use crate::constants::MU_BARY;

// use crate::transforms::calc_E_from_M;
// use crate::transforms::calc_f_from_E;

// use nalgebra::Vector3;

// pub fn ahat(phi: f64, theta: f64) -> Vector3<f64> {
//     return Vector3::new(-phi.sin(), phi.cos(), 0.0);
// }

// pub fn dhat(phi: f64, theta: f64) -> Vector3<f64> {
//     return Vector3::new(-theta.sin() * phi.cos(), -theta.sin() * phi.sin(), theta.cos());
// }

// pub fn calc_xyz_from_spherical(spherical: SphericalState) -> StateVector {

//     let xhat = spherical.phi.cos() * spherical.theta.cos();
//     let yhat = spherical.phi.sin() * spherical.theta.cos();
//     let zhat = spherical.theta.sin();

//     let rhat = Vector3::new(xhat, yhat, zhat);

//     let a = ahat(spherical.phi, spherical.theta);
//     let d = dhat(spherical.phi, spherical.theta);

//     let psi = (spherical.inc.cos() / spherical.theta.cos()).acos();
//     let ohat = psi.cos() * a + psi.sin() * d;

//     let position = spherical.r * rhat;
//     let velocity = spherical.vr * rhat + spherical.vo * ohat;


//     return StateVector { position: position, velocity: velocity };

// }