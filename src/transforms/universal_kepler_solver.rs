use crate::transforms::stumpff::{stumpff_c, stumpff_s};

/// Evaluates the universal Kepler's function f(χ) = 0.
/// 
/// This function is used in the Newton-Raphson iteration to solve Kepler's equation 
/// in universal form. 
/// 
/// # Arguments
/// * `chi` - Universal anomaly value χ
/// * `r0` - Initial radius [AU]
/// * `vr0` - Initial radial velocity [AU/day]
/// * `alpha` - Reciprocal of semi-major axis (-2E/μ) [1/AU]
/// * `mu` - Gravitational parameter [AU³/day²]
/// * `dt` - Time interval [days]
/// 
/// # Returns
/// Value of the universal Kepler function at χ
fn f(chi: f64, r0: f64, vr0: f64, alpha: f64, mu: f64, dt: f64) -> f64 {
    let z = alpha * chi.powi(2);
    let first_term = r0 * vr0 / mu.sqrt() * chi.powi(2) * stumpff_c(z);
    let second_term = (1.0 - alpha * r0) * chi.powi(3) * stumpff_s(z);
    let third_term = r0 * chi;
    let fourth_term = dt * mu.sqrt();
    first_term + second_term + third_term - fourth_term
}

/// Evaluates the derivative of the universal Kepler's function df(χ)/dχ.
/// 
/// This function calculates the derivative needed for Newton-Raphson iteration:
/// 
/// # Arguments
/// * `chi` - Universal anomaly value χ
/// * `r0` - Initial radius [AU]
/// * `vr0` - Initial radial velocity [AU/day]
/// * `alpha` - Reciprocal of semi-major axis (-2E/μ) [1/AU]
/// * `mu` - Gravitational parameter [AU³/day²]
/// 
/// # Returns
/// Value of the derivative at χ
fn df_dchi(chi: f64, r0: f64, vr0: f64, alpha: f64, mu: f64) -> f64 {
    let z = alpha * chi.powi(2);
    let first_term = r0 * vr0 / mu.sqrt() * chi * (1.0 - z * stumpff_s(z));
    let second_term = (1.0 - alpha * r0) * chi.powi(2) * stumpff_c(z);
    let third_term = r0;
    first_term + second_term + third_term
}

fn d2f_dchi2(chi: f64, r0: f64, vr0: f64, alpha: f64, mu: f64) -> f64 {
    let z = alpha * chi.powi(2);
    let dz = 2.0 * alpha * chi;
    let dS = (stumpff_c(z) - 3.0 * stumpff_s(z)) / (2.0 * z);
    let dC = (1.0 - z*stumpff_s(z) - 2.0*stumpff_c(z)) / (2.0 * z);

    let first_term = r0 * vr0 / mu.sqrt() * (1.0- 2.0*alpha*chi*stumpff_s(z) - alpha*chi.powi(2)*dS*dz);
    let second_term = 2.0*(1.0 - alpha*r0) * chi * stumpff_c(z);
    let third_term = (1.0 - alpha*r0) * chi.powi(2) * dC * dz;
    first_term + second_term + third_term
}

/// Solves the universal Kepler equation using Newton-Raphson iteration.
/// 
/// This function implements a numerical solver for the universal form of Kepler's equation,
/// which works for all orbit types (elliptical, parabolic, and hyperbolic).
/// It uses Newton-Raphson iteration to find the universal anomaly χ that satisfies
/// the time equation.
/// 
/// # Arguments
/// * `r0` - Initial radius [AU]
/// * `vr0` - Initial radial velocity [AU/day]
/// * `alpha` - Reciprocal of semi-major axis (-2E/μ) [1/AU]
/// * `mu` - Gravitational parameter [AU³/day²]
/// * `dt` - Time interval [days]
/// * `tol` - Convergence tolerance
/// * `max_iter` - Maximum number of iterations
/// 
/// # Returns
/// * `Ok(f64)` - Universal anomaly χ if solution converges
/// * `Err(Box<dyn Error>)` - Error if solution fails to converge
/// 
/// # Example
/// ```
/// let chi = solve_for_universal_anomaly(1.0, 0.1, -1.0, 1.0, 1.0, 1e-12, 1000)?;
/// ```

// We can bring this back at some point, but for now we will use Laguerre's method
// pub fn solve_for_universal_anomaly(r0: f64, vr0: f64, alpha: f64, mu: f64, dt: f64, tol: f64, max_iter: usize) -> Result<f64, Box<dyn std::error::Error>> {
//     let mut chi = mu.sqrt() * alpha.abs() * dt;
//     let mut iter = 0;
//     let mut error = f(chi, r0, vr0, alpha, mu, dt).abs();

//     while error > tol {
//         if iter > max_iter {
//             println!("\nSolver failed with:");
//             println!("r0: {}, vr0: {}, alpha: {}, mu: {}, dt: {}", r0, vr0, alpha, mu, dt);
//             println!("Initial chi: {}", mu.sqrt() * alpha.abs() * dt);
//             println!("Final chi: {}, Final error: {}", chi, error);
//             println!("z value: {}", alpha * chi.powi(2));
//             println!("Last f value: {}", f(chi, r0, vr0, alpha, mu, dt));
//             println!("Last df value: {}", df_dchi(chi, r0, vr0, alpha, mu));
//             return Err("Universal Kepler solver did not converge. Pretty bad.".into());
//         }
//         let f_val = f(chi, r0, vr0, alpha, mu, dt);
//         let df_val = df_dchi(chi, r0, vr0, alpha, mu);
//         let delta_chi = f_val / df_val;
//         chi -= delta_chi;
//         error = f(chi, r0, vr0, alpha, mu, dt).abs();
//         iter += 1;
//     }
//     Ok(chi)
// }

// Implement solver using Laguerre's method
pub fn solve_for_universal_anomaly(r0: f64, vr0: f64, alpha: f64, mu: f64, dt: f64, tol: f64, max_iter: usize) -> Result<f64, Box<dyn std::error::Error>> {
    let mut chi = mu.sqrt() * alpha.abs() * dt;
    let mut iter = 0;
    let mut error = f(chi, r0, vr0, alpha, mu, dt).abs();
    let n = 3.0; // Polynomial degree

    while error > tol {
        if iter > max_iter {
            println!("\nSolver failed with:");
            println!("r0: {}, vr0: {}, alpha: {}, mu: {}, dt: {}", r0, vr0, alpha, mu, dt);
            println!("Initial chi: {}", mu.sqrt() * alpha.abs() * dt);
            println!("Final chi: {}, Final error: {}", chi, error);
            println!("z value: {}", alpha * chi.powi(2));
            println!("Last f value: {}", f(chi, r0, vr0, alpha, mu, dt));
            println!("Last df value: {}", df_dchi(chi, r0, vr0, alpha, mu));
            return Err("Universal Kepler solver did not converge. Pretty bad.".into());
        }
        
        let G = df_dchi(chi, r0, vr0, alpha, mu) / error; 
        let H = G.powi(2) - d2f_dchi2(chi, r0, vr0, alpha, mu) / error;

        let denom_plus = G + ((n-1.0) * (n*H - G.powi(2))).sqrt();
        let denom_minus = G - ((n-1.0) * (n*H - G.powi(2))).sqrt();
        let denom = if denom_plus.abs() > denom_minus.abs() { denom_plus } else { denom_minus };

        let a = n / denom;
        chi = chi - a;
    }
   
    Ok(chi)
}