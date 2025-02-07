use crate::time::Time;
use crate::spacerock::SpaceRock;
use crate::nbody::Simulation;
use crate::observing::Observation;
use crate::observing::observation::ObservationType;

use nalgebra::{DMatrix, DVector};

use std::time::Instant;


pub fn residuals(detections: &Vec<&Observation>, theta: &[f64; 7], mut sim: Simulation) -> Result<DVector<f64>, Box<dyn std::error::Error>> {


    let mut trial = SpaceRock::from_xyz("rock", theta[0], theta[1], theta[2], theta[3], theta[4], theta[5], Time::new(theta[6], "tdb", "jd")?, "J2000", "SSB")?;
    sim.integrate(&trial.epoch);
    sim.add(trial);

    let mut residuals: DVector<f64> = DVector::zeros(detections.len());
    // let mut rock = trial.clone();
    
    for (idx, detection) in detections.iter().enumerate() {

        let observed_parameters = match detection.observation_type {
            ObservationType::Astrometry { ra, dec } => DVector::from_vec(vec![ra, dec]),
            ObservationType::Streak { ra, dec, ra_rate, dec_rate } => DVector::from_vec(vec![ra, dec, ra_rate, dec_rate]),
            _ => return Err("Observation type not supported".into()),
        };

        // somehow get the rock to the epoch of the detection
        // rock.analytic_propagate(&detection.epoch);
        sim.integrate(&detection.epoch);
        let mut rock = sim.get_particle("rock")?.clone();

        // calculate the model observations
        let astro = rock.observe(&detection.observer)?;
        

        let model_parameters = match detection.observation_type {
            ObservationType::Astrometry { ra, dec } => DVector::from_vec(vec![astro.ra(), astro.dec()]),
            ObservationType::Streak { ra, dec, ra_rate, dec_rate } => DVector::from_vec(vec![astro.ra(), astro.dec(), astro.ra_rate().expect("Must have ra rate"), astro.dec_rate().expect("Must have ra rate")]),
            _ => return Err("Observation type not supported".into()),
        };

        let d = observed_parameters - model_parameters;
        //println!("Residual: {:?}", d);

        let m = &d.transpose() * &detection.inverse_covariance.clone().unwrap() * &d;
        //println!("Chi squared: {}", m[0]);

        residuals[idx] = m[0].sqrt();

    }
    Ok(residuals)
}

// pub fn residuals_and_derivatives(detections: &Vec<&Observation>, theta: &[f64; 7], sim: Simulation) -> (DVector<f64>, DMatrix<f64>) {
//     let central_residuals = residuals(detections, theta, sim.clone()).unwrap();

//     let mut jac = DMatrix::zeros(detections.len(), 6);
//     let mut theta_plus = theta.clone();
    
//     // let eps = (f64::EPSILON).sqrt();
//     let eps = 1.0e-8;

//     for i in 0..6 {
//         theta_plus[i] += eps;
//         let res_plus = residuals(detections, &theta_plus, sim.clone()).unwrap();
//         let deriv = (res_plus - &central_residuals) / (eps);
//         for j in 0..detections.len() {
//             jac[(j, i)] = deriv[j]
//         }
//         theta_plus[i] -= eps;
//     }

//     (central_residuals, jac)
// }

pub fn residuals_and_derivatives(detections: &Vec<&Observation>, theta: &[f64; 7], sim: Simulation) -> (DVector<f64>, DMatrix<f64>) {
    let central_residuals = residuals(detections, theta, sim.clone()).unwrap();

    let mut jac = DMatrix::zeros(detections.len(), 6);
    let mut theta_plus = theta.clone();
    let mut theta_minus = theta.clone();
    
    // let eps = (f64::EPSILON).sqrt();
    let eps = 1.0e-8;

    for i in 0..6 {
        theta_plus[i] += eps;
        theta_minus[i] -= eps;
        let res_plus = residuals(detections, &theta_plus, sim.clone()).unwrap();
        let res_minus = residuals(detections, &theta_minus, sim.clone()).unwrap();
        let deriv = (res_plus - res_minus) / (2.0 * eps);
        for j in 0..detections.len() {
            jac[(j, i)] = deriv[j]
        }
        theta_plus[i] -= eps;
        theta_minus[i] += eps;
    }

    (central_residuals, jac)
}


pub fn cost(detections: &Vec<&Observation>, theta: &[f64; 7], sim: Simulation) -> f64 {
    let epoch = Time::new(theta[6], "tdb", "jd").unwrap();
    let trial = SpaceRock::from_xyz("rock", theta[0], theta[1], theta[2], theta[3], theta[4], theta[5], epoch, "J2000", "SSB");
    let res = residuals(detections, theta, sim.clone()).unwrap();
    //let csq = res.dot(&res);
    //csq
    // sum of res
    let mut s = 0.0;
    for i in 0..res.len() {
        s += res[i];
    }
    s
}



#[derive(Debug, Clone)]
pub struct FitResult {
    pub chisq: f64,
    pub rock: SpaceRock,
    pub niter: usize,
    pub dof: f64,
    // pub residuals: Vec<f64>,
    pub covariance: DMatrix<f64>,
}


pub fn fit_orbit_lm(detections: &Vec<&Observation>, initial_guess: &[f64; 7], sim: Simulation) -> Result<Option<FitResult>, Box<dyn std::error::Error>> {  

    // do timing
    let start = Instant::now();
    
    let csq_tol = 1e-3;
    let grad_tol = 1e-3;
    let theta_tol = 1e-10;
    let rho_accept = -0.01;

    // let csq_tol = 1e-3;
    // let grad_tol = 1e-15;
    // let theta_tol = 1e-15;
    // let rho_accept = 0.2;

    let maxiter = 100;
    let mut niter = 0;

    let dof = detections.len() as f64 * 2.0 - 6.0;
    let mut theta = initial_guess.clone();
    
    let mut csq = cost(detections, &theta, sim.clone());
    let mut new_csq: f64 = csq;
    
    let mut csq_change: f64 = 10000.0;
    
    let mut lambda: f64 = 0.01;

    while csq_change.abs() > csq_tol {
    // while true {

        println!("Iteration: {}, chisq: {}, lambda: {}, ndof: {}", niter, csq, lambda, dof);

        // let res = residuals(detections, &theta);
        // let mut j = jacobian(detections, &theta);

        let (res, mut j) = residuals_and_derivatives(detections, &theta, sim.clone());

        let mut a = &j.transpose() * &j;
        for idx in 0..6 {
            a[(idx, idx)] += lambda;
        }
        
        let a_inv = match a.try_inverse() {
            Some(a_inv) => a_inv,
            None => return Err("Matrix inversion failed".into())
        };

        //check if the gradient has converged
        let grad = &j.transpose() * &res;
        let mut max_grad = 0.0;
        for idx in 0..6 {
            if grad[idx].abs() > max_grad {
                max_grad = grad[idx].abs();
            }
        }
        if max_grad < grad_tol {
            println!("Gradient converged");
            break;
        }
        let h = a_inv * &grad;

        // if no parameter is being changed by more than theta_tol, break
        let mut max_h = 0.0;
        for idx in 0..6 {
            let perturbation = (h[idx] / theta[idx]).abs();
            if perturbation > max_h {
                max_h = perturbation;
            }
        }
        if max_h < theta_tol {
            println!("Parameters converged");
            break;
        }

        for idx in 0..6 {
            theta[idx] -= h[idx];
        }

        new_csq = cost(detections, &theta, sim.clone());
        csq_change = new_csq - csq;

        // let rho = (csq - new_csq) / (&h.transpose() * (lambda * &h + &grad)).norm();
        // let actual_reduction = csq - new_csq; 

        // // The predicted reduction in standard LM:
        // let predicted_reduction = - (grad.transpose() * &h + 0.5 * (h.transpose() * &j.transpose() * &j * &h))[0];

        //  let rho = actual_reduction / predicted_reduction; 

        let rho = (csq - new_csq) / (h.clone().transpose() * (lambda * h.clone() + grad.clone())).norm();
        if rho > rho_accept {
            // accept the step and shrink lambda
            lambda *= 0.1;
            csq = new_csq;
        } else {
            // reject the step, increase lambda, and reset the parameters
            lambda *= 10.0;
            for idx in 0..6 {
                theta[idx] += h[idx];
            }
        }

        niter += 1;
        if niter >= maxiter {
            println!("Max iterations reached");
            break;
        }
        
    }

    let rock = SpaceRock::from_xyz("rock", theta[0], theta[1], theta[2], theta[3], theta[4], theta[5], Time::new(theta[6], "tdb", "jd")?, "J2000", "SSB")?;

    

    let (res, mut j) = residuals_and_derivatives(detections, &theta, sim.clone());
    let a = &j.transpose() * &j;
    let cov = a.pseudo_inverse(1e-10)?;

    // let rd_resid = radec_residuals(detections, &theta, epoch);

    println!("Final chisq: {}", new_csq);
    println!("Final chisq/dof: {}", new_csq / dof);
    println!("Final parameters: {:?}", theta);
    let duration = start.elapsed();
    println!("Time elapsed in LM fit is: {:?}", duration);

    return Ok(Some(FitResult {
        chisq: new_csq * dof,
        dof: dof,
        // residuals: residuals(detections, &theta).iter().map(|r| *r).collect::<Vec<_>>(),
        rock: rock,
        niter: niter,
        covariance: cov
    }));


}


// pub fn fit_orbit_lm(detections: &Vec<&Observation>, initial_guess: &[f64; 7], sim: Simulation) -> Result<Option<FitResult>, Box<dyn std::error::Error>> {  

//     // do timing
//     let start = Instant::now();
    
//     let csq_tol = 1e-1;
//     let grad_tol = 1e-3;
//     let theta_tol = 1e-4;
//     let rho_accept = 0.1;
//     let tol = 1.0e-10;

//     let maxiter = 10;
//     let dof = detections.len() as f64 * 2.0 - 6.0;
//     let mut theta = initial_guess.clone();
    
//     let mut chisq = cost(detections, &theta, sim.clone()) / dof;
//     println!("Initial chisq: {}", chisq);
//     let mut new_csq: f64 = chisq;
    
//     let mut csq_change: f64 = 1000.0;
    
//     let mut lambda: f64 = 0.00001;

//     for iter in 0..maxiter {

//         let (r, J) = residuals_and_derivatives(detections, &theta, sim.clone());

//         // Compute the gradient g = J^T * r.
//         let g = J.transpose() * &r;

//         // Convergence check: if the gradient is small, we are done.
//         if g.norm() < tol {
//             println!("Converged (gradient norm < tol) at iteration {}.", iter);
//             break;
//         }

//         // Hessian approximation H = J^T * J.
//         let H = J.transpose() * &J;
//         // Form the damped Hessian: H_damped = H + lambda * I.
//         let H_damped = &H + lambda * DMatrix::identity(H.nrows(), H.ncols());

//         // Solve for update: h = - (H_damped)^{-1} * g.
//         let h = match H_damped.lu().solve(&(-&g)) {
//             Some(sol) => sol,
//             None => {
//                 println!("Failed to solve linear system at iteration {}.", iter);
//                 break;
//             }
//         };

//         println!("Update: {:?}", h);

//         // Candidate new parameters.
//         // let theta_new = &theta + &h;
//         let mut theta_new = theta.clone();
//         for idx in 0..6 {
//             theta_new[idx] += h[idx];
//         }
//         let r_new = residuals(detections, &theta_new, sim.clone())?;
//         let chisq_new = r_new.norm_squared();
//         println!("New chisq: {}", chisq_new);

//         // Compute the actual reduction.
//         let actual_reduction = chisq - chisq_new;
//         // Compute the predicted reduction from the quadratic model:
//         // predicted = - [g^T h + 0.5 * h^T H h]
//         let predicted_reduction = - (g.dot(&h) + 0.5 * &h.dot(&(H * &h)));

//         // Compute the ratio of actual to predicted reduction.
//         let rho = actual_reduction / predicted_reduction;

//         // Decide whether to accept the step.
//         if rho > 0.0 {
//             // Accept the update.
//             theta = theta_new;
//             chisq = chisq_new;
//             // Optionally reduce lambda.
//             lambda *= 0.1;
//         } else {
//             // Reject the update and increase lambda.
//             lambda *= 10.0;
//         }

//         println!(
//             "Iteration {}: chisq = {}, lambda = {}, rho = {}",
//             iter, chisq, lambda, rho
//         );

//         // Convergence check: if the update h is very small.
//         let mut theta_norm_sq = 0.0;
//         for idx in 0..6 {
//             theta_norm_sq += theta[idx] * theta[idx];
//         }
//         let theta_norm = theta_norm_sq.sqrt();

//         if h.norm() < tol * (theta_norm + tol) {
//             println!("Converged (parameter update small) at iteration {}.", iter);
//             break;
//         }
        
//     }

//     let rock = SpaceRock::from_xyz("rock", theta[0], theta[1], theta[2], theta[3], theta[4], theta[5], Time::new(theta[6], "tdb", "jd")?, "J2000", "SSB")?; 

//     let (res, mut j) = residuals_and_derivatives(detections, &theta, sim.clone());
//     let a = &j.transpose() * &j;
//     let cov = a.pseudo_inverse(1e-10)?;

//     // let rd_resid = radec_residuals(detections, &theta, epoch);

//     println!("Final chisq: {}", new_csq * dof);
//     println!("Final parameters: {:?}", theta);
//     let duration = start.elapsed();
//     println!("Time elapsed in LM fit is: {:?}", duration);

//     return Ok(Some(FitResult {
//         chisq: new_csq * dof,
//         dof: dof,
//         // residuals: residuals(detections, &theta).iter().map(|r| *r).collect::<Vec<_>>(),
//         rock: rock,
//         niter: 0,
//         covariance: cov
//     }));


// }