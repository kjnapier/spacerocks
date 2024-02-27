// use crate::observing::Detection;
// use crate::constants::{MU_BARY, SPEED_OF_LIGHT};
// use crate::StateVector;
// use crate::SpaceRock;

// use nalgebra::Matrix3;
// use nalgebra::matrix;

// pub fn gauss_fit(dets: &Vec<&Detection>, min_distance: f64) -> Option<Vec<SpaceRock>> {
//     let mut detections = dets.clone();
//     detections.sort_by(|b, a| a.epoch.epoch.partial_cmp(&b.epoch.epoch).unwrap());
//     let first = &detections[0];
//     let last = &detections[detections.len() - 1];
//     let middle = &detections[detections.len() / 2];
//     let triplet = vec![first.clone(), middle.clone(), last.clone()];

//     return gauss(&triplet, min_distance);
// }

// pub fn gauss(triplet: &Vec<&Detection>, min_distance: f64) -> Option<Vec<SpaceRock>> {

//     let R1 = triplet[0].observer.position;
//     let R2 = triplet[1].observer.position;
//     let R3 = triplet[2].observer.position;

//     let rho1 = triplet[0].pointing();
//     let rho2 = triplet[1].pointing();
//     let rho3 = triplet[2].pointing();
    
//     let t1 = triplet[0].epoch.epoch;
//     let t2 = triplet[1].epoch.epoch;
//     let t3 = triplet[2].epoch.epoch;

//     let tau1 = t1 - t2;
//     let tau3 = t3 - t2;
//     let tau = t3 - t1;

//     let p1 = rho2.cross(&rho3);
//     let p2 = rho1.cross(&rho3);
//     let p3 = rho1.cross(&rho2);

//     let D0 = rho1.dot(&p1);

//     let D: Matrix3<f64> = Matrix3::new(
//         R1.dot(&p1), R1.dot(&p2), R1.dot(&p3),
//         R2.dot(&p1), R2.dot(&p2), R2.dot(&p3),
//         R3.dot(&p1), R3.dot(&p2), R3.dot(&p3),
//     );

//     // get the item in the first row, second column
//     let A = (1.0/D0) * (-D[(0,1)] * (tau3/tau) + D[(1,1)] + D[(2,1)] * (tau1/tau));
//     let B = (1.0/(6.0 * D0)) * (D[(0,1)] * (tau3.powi(2) - tau.powi(2)) * (tau3/tau) + D[(2,1)] * (tau.powi(2) - tau1.powi(2)) * (tau1/tau));
//     let E = R2.dot(&rho2);

//     let R2sq = R2.dot(&R2);

//     let a = -(A.powi(2) + 2.0 * A * E + R2sq);
//     let b = -2.0 * MU_BARY * B * (A + E);
//     let c = -MU_BARY.powi(2) * B.powi(2);

//     let mat = matrix![0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
//                       0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0;
//                       0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0;
//                       0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0;
//                       0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0;
//                       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0;
//                       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0;
//                       -c,  0.0, 0.0, -b,  0.0, 0.0, -a,  0.0];

//     let mat = match mat.try_schur(0.00001, 10000) {
//         Some(mat) => mat,
//         None => return None,
//     };

//     let complex_roots = mat.complex_eigenvalues();
//     let roots: Vec<f64> = complex_roots.iter().filter(|x| x.im == 0.0 && x.re > min_distance).map(|x| x.re).collect();
//     if roots.len() == 0 {
//         return None;
//     }


//     // create an empty vector to hold the orbits
//     let mut res: Vec<SpaceRock> = Vec::new();
//     for root in &roots {

//         let a1 = (1.0/D0) * ((6.0 * (D[(2,0)] * (tau1/tau3) + D[(1,0)] * (tau/tau3)) * root.powi(3) + MU_BARY * D[(2,0)] * (tau.powi(2) - tau1.powi(2)) * (tau1/tau3)) / (6.0 * root.powi(3) + MU_BARY * (tau.powi(2) - tau3.powi(2))) - D[(0,0)]);
//         let a2 = A + (MU_BARY * B) / root.powi(3);
//         let a3 = (1.0/D0) * ((6.0 * (D[(0,2)] * (tau3/tau1) - D[(1,2)] * (tau/tau1)) * root.powi(3) + MU_BARY * D[(0,2)] * (tau.powi(2) - tau3.powi(2)) * (tau3/tau1)) / (6.0 * root.powi(3) + MU_BARY * (tau.powi(2) - tau1.powi(2))) - D[(2,2)]);

//         let r1 = R1 + a1 * rho1;
//         let r2 = R2 + a2 * rho2;
//         let r3 = R3 + a3 * rho3;

//         let f1 = 1.0 - 0.5 * (MU_BARY/root.powi(3)) * tau1.powi(2);
//         let f3 = 1.0 - 0.5 * (MU_BARY/root.powi(3)) * tau3.powi(2);
//         let g1 = tau1 - (1.0/6.0) * (MU_BARY / root.powi(3)) * tau1.powi(3);
//         let g3 = tau3 - (1.0/6.0) * (MU_BARY / root.powi(3)) * tau3.powi(3);

//         let v2 = (-f3 * r1 + f1 * r3) / (f1 * g3 - f3 * g1);

//         let state = StateVector::new(r2.x, r2.y, r2.z, v2.x, v2.y, v2.z);
//         let ltt = state.position.norm() / SPEED_OF_LIGHT;
//         let mut corrected_t = triplet[1].epoch.clone();
//         corrected_t.epoch -= ltt;
//         let rock = SpaceRock::from_state("rock", state, corrected_t, "J2000", "SSB");
//         res.push(rock);
//     }

//     return Some(res);

// }