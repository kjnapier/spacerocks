use crate::{Observation, SpaceRock};
use crate::data::{SPEED_OF_LIGHT, MU_BARY};

use nalgebra::Matrix3;
use nalgebra::matrix;


/// Implements Gauss' method for initial orbit determination from three observations.
/// This method uses three on-sky measurements to determine possible orbital solutions.
///
/// # Arguments
/// * `o1`, `o2`, `o3` - Three observations, ordered by time
/// * `min_distance` - Minimum distance (in AU) to consider valid solutions
///
/// # Returns
/// * `Some(Vec<SpaceRock>)` - Vector of possible orbit solutions
/// * `None` - If no valid solutions found
///
/// The method:
/// 1. Orders observations chronologically
/// 2. Solves 8th degree polynomial for orbital distance
/// 3. Uses roots to generate candidate orbits
/// 4. Applies light-time correction to final solutions
pub fn gauss(o1: &Observation, o2: &Observation, o3: &Observation, min_distance: f64) -> Option<Vec<SpaceRock>> {

    // get the order of the epochs
    let mut triplet = [o1, o2, o3];
    triplet.sort_by(|a, b| a.epoch.epoch.partial_cmp(&b.epoch.epoch).unwrap());

    let o1 = triplet[0].observer.position;
    let o2 = triplet[1].observer.position;
    let o3 = triplet[2].observer.position;

    let rho1 = triplet[0].pointing();
    let rho2 = triplet[1].pointing();
    let rho3 = triplet[2].pointing();
    
    let t1 = triplet[0].epoch.utc().jd();
    let t2 = triplet[1].epoch.utc().jd();
    let t3 = triplet[2].epoch.utc().jd();

    let tau1 = t1 - t2;
    let tau3 = t3 - t2;
    let tau = t3 - t1;

    let p1 = rho2.cross(&rho3);
    let p2 = rho1.cross(&rho3);
    let p3 = rho1.cross(&rho2);

    let d0 = rho1.dot(&p1);

    let d: Matrix3<f64> = Matrix3::new(
        o1.dot(&p1), o1.dot(&p2), o1.dot(&p3),
        o2.dot(&p1), o2.dot(&p2), o2.dot(&p3),
        o3.dot(&p1), o3.dot(&p2), o3.dot(&p3),
    );

    // get the item in the first row, second column
    let a = (1.0/d0) * (-d[(0,1)] * (tau3/tau) + d[(1,1)] + d[(2,1)] * (tau1/tau));
    let b = (1.0/(6.0 * d0)) * (d[(0,1)] * (tau3.powi(2) - tau.powi(2)) * (tau3/tau) + d[(2,1)] * (tau.powi(2) - tau1.powi(2)) * (tau1/tau));
    let e = o2.dot(&rho2);

    let o2sq = o2.dot(&o2);

    let aa = -(a.powi(2) + 2.0 * a * e + o2sq);
    let bb = -2.0 * MU_BARY * b * (a + e);
    let cc = -MU_BARY.powi(2) * b.powi(2);

    let mat = matrix![0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
                      0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0;
                      0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0;
                      0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0;
                      0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0;
                      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0;
                      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0;
                      -cc, 0.0, 0.0, -bb, 0.0, 0.0, -aa, 0.0];

    let mat = match mat.try_schur(0.00001, 10000) {
        Some(mat) => mat,
        None => return None,
    };

    let complex_roots = mat.complex_eigenvalues();
    let roots: Vec<f64> = complex_roots.iter().filter(|x| x.im == 0.0 && x.re > min_distance).map(|x| x.re).collect();
    if roots.len() == 0 {
        return None;
    }


    // create an empty vector to hold the orbits
    let mut res: Vec<SpaceRock> = Vec::new();
    for root in &roots {

        let a1 = (1.0/d0) * ((6.0 * (d[(2,0)] * (tau1/tau3) + d[(1,0)] * (tau/tau3)) * root.powi(3) + MU_BARY * d[(2,0)] * (tau.powi(2) - tau1.powi(2)) * (tau1/tau3)) / (6.0 * root.powi(3) + MU_BARY * (tau.powi(2) - tau3.powi(2))) - d[(0,0)]);
        let a2 = a + (MU_BARY * b) / root.powi(3);
        let a3 = (1.0/d0) * ((6.0 * (d[(0,2)] * (tau3/tau1) - d[(1,2)] * (tau/tau1)) * root.powi(3) + MU_BARY * d[(0,2)] * (tau.powi(2) - tau3.powi(2)) * (tau3/tau1)) / (6.0 * root.powi(3) + MU_BARY * (tau.powi(2) - tau1.powi(2))) - d[(2,2)]);

        let r1 = o1 + a1 * rho1;
        let r2 = o2 + a2 * rho2;
        let r3 = o3 + a3 * rho3;

        let f1 = 1.0 - 0.5 * (MU_BARY/root.powi(3)) * tau1.powi(2);
        let f3 = 1.0 - 0.5 * (MU_BARY/root.powi(3)) * tau3.powi(2);
        let g1 = tau1 - (1.0/6.0) * (MU_BARY / root.powi(3)) * tau1.powi(3);
        let g3 = tau3 - (1.0/6.0) * (MU_BARY / root.powi(3)) * tau3.powi(3);

        let v2 = (-f3 * r1 + f1 * r3) / (f1 * g3 - f3 * g1);

        let x = r2.x;
        let y = r2.y;
        let z = r2.z;
        let vx = v2.x;
        let vy = v2.y;
        let vz = v2.z;

        let ltt = r2.norm() / SPEED_OF_LIGHT;
        let mut corrected_t = triplet[1].epoch.clone();
        corrected_t.epoch -= ltt;
        let rock = SpaceRock::from_xyz("rock", x, y, z, vx, vy, vz, corrected_t, "J2000", "SSb").expect("Failed to create SpaceRock from XYZ");
        res.push(rock);
    }

    return Some(res);

}