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

// / Solves the universal Kepler equation using Newton-Raphson iteration.
// / 
// / This function implements a numerical solver for the universal form of Kepler's equation,
// / which works for all orbit types (elliptical, parabolic, and hyperbolic).
// / It uses Newton-Raphson iteration to find the universal anomaly χ that satisfies
// / the time equation.
// / 
// / # Arguments
// / * `r0` - Initial radius [AU]
// / * `vr0` - Initial radial velocity [AU/day]
// / * `alpha` - Reciprocal of semi-major axis (-2E/μ) [1/AU]
// / * `mu` - Gravitational parameter [AU³/day²]
// / * `dt` - Time interval [days]
// / * `tol` - Convergence tolerance
// / * `max_iter` - Maximum number of iterations
// / 
// / # Returns
// / * `Ok(f64)` - Universal anomaly χ if solution converges
// / * `Err(Box<dyn Error>)` - Error if solution fails to converge
// / 
// / # Example
// / ```
// / let chi = solve_for_universal_anomaly(1.0, 0.1, -1.0, 1.0, 1.0, 1e-12, 1000)?;
// / ```

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

// // Implement solver using Laguerre's method
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
        error = f(chi, r0, vr0, alpha, mu, dt).abs();
        iter += 1;
    }
   
    Ok(chi)
}

// use std::error::Error;
// use std::io;

// /// Stumpff functions c_k(x) for k = 0,1,2,3.
// /// Returns (c0, c1, c2, c3).
// pub fn stumpff(mut x: f64) -> (f64, f64, f64, f64) {
//     let mut n: i32 = 0;
//     let xm = 0.1_f64;

//     while x.abs() > xm {
//         n += 1;
//         x /= 4.0;
//     }

//     // Match the exact Horner polynomials in the Python
//     let mut d2 = (1.0
//         - x * (1.0
//         - x * (1.0
//         - x * (1.0
//         - x * (1.0
//         - x * (1.0 - x / 182.0) / 132.0) / 90.0) / 56.0) / 30.0) / 12.0) / 2.0;

//     let mut d3 = (1.0
//         - x * (1.0
//         - x * (1.0
//         - x * (1.0
//         - x * (1.0
//         - x * (1.0 - x / 210.0) / 156.0) / 110.0) / 72.0) / 42.0) / 20.0) / 6.0;

//     let mut d1 = 1.0 - x * d3;
//     let mut d0 = 1.0 - x * d2;

//     while n > 0 {
//         n -= 1;
//         d3 = (d2 + d0 * d3) / 4.0;
//         d2 = (d1 * d1) / 2.0;
//         d1 = d0 * d1;
//         d0 = 2.0 * d0 * d0 - 1.0;
//     }

//     (d0, d1, d2, d3)
// }

// /// Root function for the universal Kepler equation.
// /// Returns (f, fp, fpp, fppp).
// pub fn root_function(s: f64, mu: f64, alpha: f64, r0: f64, r0dot: f64, t: f64)
//     -> (f64, f64, f64, f64)
// {
//     let (c0, c1, c2, c3) = stumpff(alpha * s * s);
//     let zeta = mu - alpha * r0;

//     let f   = r0 * s * c1 + r0 * r0dot * s * s * c2 + mu * s * s * s * c3 - t;
//     let fp  = r0 * c0 + r0 * r0dot * s * c1 + mu * s * s * c2; // ≡ r(s)
//     let fpp = zeta * s * c1 + r0 * r0dot * c0;
//     let fppp= zeta * c0 - r0 * r0dot * alpha * s * c1;

//     (f, fp, fpp, fppp)
// }

// /// Faithful port of the Python halley_safe:
// /// Returns (converged, root_or_nan, fp_at_return).
// pub fn halley_safe(
//     x1: f64,
//     x2: f64,
//     mu: f64,
//     alpha: f64,
//     r0: f64,
//     r0dot: f64,
//     t: f64,
//     xacc: f64,
//     maxit: usize,
// ) -> (bool, f64, f64) {
//     // Use these values later (match Python slicing [0:3])
//     let (mut fl, fpl, fppl, _) = root_function(x1, mu, alpha, r0, r0dot, t);
//     let (mut fh, fph, fpph, _) = root_function(x2, mu, alpha, r0, r0dot, t);

//     // Verify bracket
//     if (fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0) {
//         return (false, f64::NAN, fl);
//     }
//     if fl == 0.0 {
//         return (true, x1, fpl);
//     }
//     if fh == 0.0 {
//         return (true, x2, fph);
//     }

//     // Orient so that f(xl) < 0 < f(xh)
//     let (mut xl, mut xh) = if fl < 0.0 { (x1, x2) } else { (x2, x1) };

//     // Choose endpoint with smaller |f|
//     let (mut rts, mut f, mut fp, mut fpp) = if fl.abs() < fh.abs() {
//         (xl, fl, fpl, fppl)
//     } else {
//         (xh, fh, fph, fpph)
//     };

//     // Initialize midpoint, steps
//     rts = 0.5 * (x1 + x2);
//     let mut dxold = (x2 - x1).abs();
//     let mut dx = dxold;

//     let (f0, fp0, fpp0, _) = root_function(rts, mu, alpha, r0, r0dot, t);
//     f = f0; fp = fp0; fpp = fpp0;

//     for _ in 0..maxit {
//         // Safeguard condition: fall back to bisection
//         if (((rts - xh) * fp - f) * ((rts - xl) * fp - f) > 0.0)
//             || ( (2.0 * f).abs() > (dxold * fp).abs() )
//         {
//             dxold = dx;
//             dx = 0.5 * (xh - xl);
//             rts = xl + dx;

//             let rel = if rts != 0.0 { (dx / rts).abs() } else { dx.abs() };
//             if rel < xacc {
//                 return (true, rts, fp);
//             }
//         } else {
//             dxold = dx;
//             // Python code computes Newton step then overwrites with Halley step
//             // dx = f / fp
//             dx = (2.0 * f * fp) / (2.0 * fp * fp - f * fpp);  // Halley
//             rts -= dx;

//             let rel = if rts != 0.0 { (dx / rts).abs() } else { dx.abs() };
//             if rel < xacc {
//                 return (true, rts, fp);
//             }
//         }

//         if {
//             let rel = if rts != 0.0 { (dx / rts).abs() } else { dx.abs() };
//             rel < xacc
//         } {
//             return (true, rts, fp);
//         }

//         let (ff, fp_new, fpp_new, _) = root_function(rts, mu, alpha, r0, r0dot, t);
//         f = ff; fp = fp_new; fpp = fpp_new;

//         // Maintain the bracket
//         if f < 0.0 {
//             xl = rts;
//             fl = f;
//         } else {
//             xh = rts;
//             fh = f;
//         }
//     }

//     (false, f64::NAN, fp)
// }

// /// Public API: solve for the universal anomaly `s` with the requested signature.
// /// Sets a bracket (x1=0, expanding x2 along the sign of dt) and calls `halley_safe`.
// pub fn solve_for_universal_anomaly(
//     r0: f64,
//     vr0: f64,
//     alpha: f64,
//     mu: f64,
//     dt: f64,
//     tol: f64,
//     max_iter: usize,
// ) -> Result<f64, Box<dyn Error>> {
//     if !r0.is_finite() || !vr0.is_finite() || !alpha.is_finite() || !mu.is_finite() || !dt.is_finite() {
//         return Err(Box::new(io::Error::new(io::ErrorKind::InvalidInput, "non-finite input")));
//     }
//     if mu <= 0.0 {
//         return Err(Box::new(io::Error::new(io::ErrorKind::InvalidInput, "mu must be > 0")));
//     }
//     if tol <= 0.0 {
//         return Err(Box::new(io::Error::new(io::ErrorKind::InvalidInput, "tol must be > 0")));
//     }
//     if dt == 0.0 {
//         return Ok(0.0);
//     }

//     // Bracket: f(0) = -dt. March x2 outward until sign flips.
//     let x1 = 0.0;
//     let mut x2 = if dt.abs() > 0.0 { dt } else { dt.signum() }; // start with ~dt
//     if x2 == 0.0 { x2 = if dt > 0.0 { 1.0 } else { -1.0 }; }

//     let f0 = -dt; // root_function(0, ...).f
//     let mut f2 = root_function(x2, mu, alpha, r0, vr0, dt).0;

//     let mut k = 0usize;
//     while f0.signum() == f2.signum() && k < 200 {
//         x2 *= 2.0;
//         f2 = root_function(x2, mu, alpha, r0, vr0, dt).0;
//         k += 1;
//     }
//     if f0.signum() == f2.signum() {
//         return Err(Box::new(io::Error::new(
//             io::ErrorKind::Other,
//             "failed to bracket the universal anomaly root",
//         )));
//     }

//     let (ok, s, _fp) = halley_safe(x1, x2, mu, alpha, r0, vr0, dt, tol, max_iter);
//     if ok && s.is_finite() {
//         Ok(s)
//     } else {
//         Err(Box::new(io::Error::new(
//             io::ErrorKind::Other,
//             "Halley solver did not converge",
//         )))
//     }
// }
