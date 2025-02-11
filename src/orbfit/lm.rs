/// Configuration parameters for the Levenberg-Marquardt optimization algorithm
pub struct LevenbergMarquardt {
     /// Maximum number of iterations before stopping
    pub max_iter: usize,
    /// Convergence threshold for the parameter changes
    pub param_tol: f64,
    /// Convergence threshold for the gradient
    pub grad_tol: f64,
    /// Initial damping parameter
    pub lambda: f64,
    /// Minimum acceptable reduction ratio for step acceptance
    pub rho_accept: f64,
}

// implement default values for the LevenbergMarquardt struct
impl LevenbergMarquardt {
    /// Creates a new LevenbergMarquardt optimizer with specified parameters
    ///
    /// # Arguments
    /// * `max_iter` - Maximum iterations before stopping
    /// * `param_tol` - Parameter change tolerance for convergence
    /// * `grad_tol` - Gradient norm tolerance for convergence
    /// * `lambda` - Initial damping parameter value
    /// * `rho_accept` - Minimum reduction ratio to accept step
    pub fn new(max_iter: usize, param_tol: f64, grad_tol: f64, lambda: f64, rho_accept: f64) -> Self {
        Self { max_iter, param_tol, grad_tol, lambda, rho_accept }
    }

    /// Creates a LevenbergMarquardt optimizer with default parameters
    pub fn default() -> Self {
        Self {
            max_iter: 1_000,
            param_tol: 1e-6,
            grad_tol: 1e-6,
            lambda: 1e-3,
            rho_accept: 0.0,
        }
    }

}


impl Minimizer for LevenbergMarquardt {
    /// Minimize the cost function using the Levenberg-Marquardt algorithm
    ///
    /// # Arguments
    /// * `model` - Model to fit to the data
    /// * `initial_guess` - Initial guess for the parameters
    /// * `data` - Data to fit the model to
    ///
    /// # Returns
    /// Vector of optimized parameters
    
    fn minimize(&self, model: &impl Model, initial_guess: Vec<f64>, data: Vec<f64>) -> Vec<f64> {

    let dof = detections.len() as f64 * 2.0 - 6.0;
    let mut theta = initial_guess.clone();
    
    // let mut csq = cost(detections, &theta, sim.clone());
    let (mut res, mut j) = residuals_and_derivatives(detections, &theta, sim.clone());
    let mut grad = &j.transpose() * &res;
    let mut a = &j.transpose() * &j;

    let mut csq = res.sum();
    let mut new_csq: f64 = csq.clone();

    let mut lambda: f64 = 0.001;
    let eye = DMatrix::identity(6, 6);

    while grad.norm() > grad_tol {

        println!("Iteration: {}, chisq: {}, lambda: {}, ndof: {}", niter, csq, lambda, dof);

        let h = match (&a + lambda * &eye).lu().solve(&(-&grad)) {
            Some(h) => h,
            None => return Err("Matrix inversion failed".into())
        };
        if h.norm() < theta_tol {
            println!("Parameters converged");
            break;
        }

        let mut theta_new = theta.clone();
        for idx in 0..6 {
            theta_new[idx] += h[idx];
        }

        let (res2, j2) = residuals_and_derivatives(detections, &theta_new, sim.clone());
        new_csq = res2.sum();

        let rho = (csq - new_csq) / (&h.transpose() * (lambda * &h - &grad)).norm();
        if rho > rho_accept {
            res = res2;
            theta = theta_new;
            j = j2;
            grad = &j.transpose() * &res;
            a = &j.transpose() * &j;
            lambda *= 0.1;
            csq = new_csq;
        } else {
            lambda *= 10.0;
        }

        niter += 1;
        if niter >= maxiter {
            println!("Max iterations reached");
            break;
        }
        
    }
        
    }
}