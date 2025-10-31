// use crate::transforms::stumpff::{stumpff_c, stumpff_s};

// /// Evaluates the universal Kepler's function f(χ) = 0.
// /// 
// /// This function is used in the Newton-Raphson iteration to solve Kepler's equation 
// /// in universal form. 
// /// 
// /// # Arguments
// /// * `chi` - Universal anomaly value χ
// /// * `r0` - Initial radius [AU]
// /// * `vr0` - Initial radial velocity [AU/day]
// /// * `alpha` - Reciprocal of semi-major axis (-2E/μ) [1/AU]
// /// * `mu` - Gravitational parameter [AU³/day²]
// /// * `dt` - Time interval [days]
// /// 
// /// # Returns
// /// Value of the universal Kepler function at χ
// fn f(chi: f64, r0: f64, vr0: f64, alpha: f64, mu: f64, dt: f64) -> f64 {
//     let z = alpha * chi.powi(2);
//     let first_term = r0 * vr0 / mu.sqrt() * chi.powi(2) * stumpff_c(z);
//     let second_term = (1.0 - alpha * r0) * chi.powi(3) * stumpff_s(z);
//     let third_term = r0 * chi;
//     let fourth_term = dt * mu.sqrt();
//     first_term + second_term + third_term - fourth_term
// }

// /// Evaluates the derivative of the universal Kepler's function df(χ)/dχ.
// /// 
// /// This function calculates the derivative needed for Newton-Raphson iteration:
// /// 
// /// # Arguments
// /// * `chi` - Universal anomaly value χ
// /// * `r0` - Initial radius [AU]
// /// * `vr0` - Initial radial velocity [AU/day]
// /// * `alpha` - Reciprocal of semi-major axis (-2E/μ) [1/AU]
// /// * `mu` - Gravitational parameter [AU³/day²]
// /// 
// /// # Returns
// /// Value of the derivative at χ
// fn df_dchi(chi: f64, r0: f64, vr0: f64, alpha: f64, mu: f64) -> f64 {
//     let z = alpha * chi.powi(2);
//     let first_term = r0 * vr0 / mu.sqrt() * chi * (1.0 - z * stumpff_s(z));
//     let second_term = (1.0 - alpha * r0) * chi.powi(2) * stumpff_c(z);
//     let third_term = r0;
//     first_term + second_term + third_term
// }

// fn d2f_dchi2(chi: f64, r0: f64, vr0: f64, alpha: f64, mu: f64) -> f64 {
//     let z = alpha * chi.powi(2);
//     let dz = 2.0 * alpha * chi;
//     let dS = (stumpff_c(z) - 3.0 * stumpff_s(z)) / (2.0 * z);
//     let dC = (1.0 - z*stumpff_s(z) - 2.0*stumpff_c(z)) / (2.0 * z);

//     let first_term = r0 * vr0 / mu.sqrt() * (1.0- 2.0*alpha*chi*stumpff_s(z) - alpha*chi.powi(2)*dS*dz);
//     let second_term = 2.0*(1.0 - alpha*r0) * chi * stumpff_c(z);
//     let third_term = (1.0 - alpha*r0) * chi.powi(2) * dC * dz;
//     first_term + second_term + third_term
// }

// // / Solves the universal Kepler equation using Newton-Raphson iteration.
// // / 
// // / This function implements a numerical solver for the universal form of Kepler's equation,
// // / which works for all orbit types (elliptical, parabolic, and hyperbolic).
// // / It uses Newton-Raphson iteration to find the universal anomaly χ that satisfies
// // / the time equation.
// // / 
// // / # Arguments
// // / * `r0` - Initial radius [AU]
// // / * `vr0` - Initial radial velocity [AU/day]
// // / * `alpha` - Reciprocal of semi-major axis (-2E/μ) [1/AU]
// // / * `mu` - Gravitational parameter [AU³/day²]
// // / * `dt` - Time interval [days]
// // / * `tol` - Convergence tolerance
// // / * `max_iter` - Maximum number of iterations
// // / 
// // / # Returns
// // / * `Ok(f64)` - Universal anomaly χ if solution converges
// // / * `Err(Box<dyn Error>)` - Error if solution fails to converge
// // / 
// // / # Example
// // / ```
// // / let chi = solve_for_universal_anomaly(1.0, 0.1, -1.0, 1.0, 1.0, 1e-12, 1000)?;
// // / ```

// // We can bring this back at some point, but for now we will use Laguerre's method
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
// // pub fn solve_for_universal_anomaly(r0: f64, vr0: f64, alpha: f64, mu: f64, dt: f64, tol: f64, max_iter: usize) -> Result<f64, Box<dyn std::error::Error>> {
// //     let mut chi = mu.sqrt() * alpha.abs() * dt;
// //     let mut iter = 0;
// //     let mut error = f(chi, r0, vr0, alpha, mu, dt).abs();
// //     let n = 3.0; // Polynomial degree

// //     while error > tol {
// //         if iter > max_iter {
// //             println!("\nSolver failed with:");
// //             println!("r0: {}, vr0: {}, alpha: {}, mu: {}, dt: {}", r0, vr0, alpha, mu, dt);
// //             println!("Initial chi: {}", mu.sqrt() * alpha.abs() * dt);
// //             println!("Final chi: {}, Final error: {}", chi, error);
// //             println!("z value: {}", alpha * chi.powi(2));
// //             println!("Last f value: {}", f(chi, r0, vr0, alpha, mu, dt));
// //             println!("Last df value: {}", df_dchi(chi, r0, vr0, alpha, mu));
// //             return Err("Universal Kepler solver did not converge. Pretty bad.".into());
// //         }
        
// //         let G = df_dchi(chi, r0, vr0, alpha, mu) / error; 
// //         let H = G.powi(2) - d2f_dchi2(chi, r0, vr0, alpha, mu) / error;

// //         let denom_plus = G + ((n-1.0) * (n*H - G.powi(2))).sqrt();
// //         let denom_minus = G - ((n-1.0) * (n*H - G.powi(2))).sqrt();
// //         let denom = if denom_plus.abs() > denom_minus.abs() { denom_plus } else { denom_minus };

// //         let a = n / denom;
// //         chi = chi - a;
// //         error = f(chi, r0, vr0, alpha, mu, dt).abs();
// //         iter += 1;
// //     }
   
// //     Ok(chi)
// // }


use std::error::Error;
use std::io;

/// Computes the Stumpff functions c_k(x) for k = 0, 1, 2, 3.
/// Returns (c0, c1, c2, c3).
pub fn stumpff(mut x: f64) -> (f64, f64, f64, f64) {
    let mut n: i32 = 0;
    let xm = 0.1_f64;

    while x.abs() > xm {
        n += 1;
        x *= 0.25;
    }

    let mut c2 = (1.0
        - x * (1.0
        - x * (1.0
        - x * (1.0
        - x * (1.0
        - x * (1.0 - x / 182.0) / 132.0) / 90.0) / 56.0) / 30.0) / 12.0) / 2.0;

    let mut c3 = (1.0
        - x * (1.0
        - x * (1.0
        - x * (1.0
        - x * (1.0
        - x * (1.0 - x / 210.0) / 156.0) / 110.0) / 72.0) / 42.0) / 20.0) / 6.0;

    let mut c1 = 1.0 - x * c3;
    let mut c0 = 1.0 - x * c2;

    // Undo argument scaling
    while n > 0 {
        n -= 1;
        c3 = (c2 + c0 * c3) * 0.25;
        c2 = 0.5 * (c1 * c1);
        c1 = c0 * c1;
        c0 = 2.0 * c0 * c0 - 1.0;
    }

    (c0, c1, c2, c3)
}

/// Root function for the universal Kepler equation.
/// Returns (f, f', f'', f''').
pub fn root_function(s: f64, mu: f64, alpha: f64, r0: f64, vr0: f64, dt: f64) -> (f64, f64, f64, f64) {
    let (c0, c1, c2, c3) = stumpff(alpha * s * s);
    let zeta = mu - alpha * r0;

    let f    = r0 * s * c1 + r0 * vr0 * s * s * c2 + mu * s * s * s * c3 - dt;
    let fp   = r0 * c0 + r0 * vr0 * s * c1 + mu * s * s * c2;
    let fpp  = zeta * s * c1 + r0 * vr0 * c0;
    let fppp = zeta * c0 - r0 * vr0 * alpha * s * c1;

    (f, fp, fpp, fppp)
}

/// Solve the universal Kepler equation for the universal anomaly `s`.
pub fn solve_for_universal_anomaly(
    r0: f64,
    vr0: f64,
    alpha: f64,
    mu: f64,
    dt: f64,
    tol: f64,
    max_iter: usize,
) -> Result<f64, Box<dyn Error>> {
    if !r0.is_finite() || !vr0.is_finite() || !alpha.is_finite() || !mu.is_finite() || !dt.is_finite() {
        return Err(Box::new(io::Error::new(io::ErrorKind::InvalidInput, "non-finite input")));
    }
    if mu <= 0.0 {
        return Err(Box::new(io::Error::new(io::ErrorKind::InvalidInput, "mu must be > 0")));
    }
    if tol <= 0.0 {
        return Err(Box::new(io::Error::new(io::ErrorKind::InvalidInput, "tol must be > 0")));
    }
    if dt == 0.0 {
        return Ok(0.0);
    }

    // Bracket the root
    let mut s_lo = 0.0;
    let mut s_hi = dt;
    if s_hi == 0.0 {
        s_hi = dt.signum();
    }

    let f_lo = -dt;
    let mut f_hi = root_function(s_hi, mu, alpha, r0, vr0, dt).0;

    let mut expand_iter = 0usize;
    const MAX_EXPAND: usize = 200;
    while f_lo.signum() == f_hi.signum() && expand_iter < MAX_EXPAND {
        s_hi *= 2.0;
        f_hi = root_function(s_hi, mu, alpha, r0, vr0, dt).0;
        expand_iter += 1;
    }

    if f_lo.signum() == f_hi.signum() {
        return Err(Box::new(io::Error::new(io::ErrorKind::Other, "failed to bracket root")));
    }

    let (mut xl, mut xh, mut fl, mut fh) = if f_lo < 0.0 {
        (s_lo, s_hi, f_lo, f_hi)
    } else {
        (s_hi, s_lo, f_hi, f_lo)
    };

    // Halley-with-bisection loop
    let mut rts = 0.5 * (xl + xh);
    let mut dxold = (xh - xl).abs();
    let mut dx = dxold;

    let (mut f, mut fp, mut fpp, _) = root_function(rts, mu, alpha, r0, vr0, dt);

    for _ in 0..max_iter {
        let use_bisect = ((rts - xh) * fp - f) * ((rts - xl) * fp - f) > 0.0
            || (2.0 * f).abs() > (dxold * fp).abs();

        if use_bisect {
            dxold = dx;
            dx = 0.5 * (xh - xl);
            rts = xl + dx;
            let rel = if rts != 0.0 { (dx / rts).abs() } else { dx.abs() };
            if rel < tol {
                return Ok(rts);
            }
        } else {
            dxold = dx;
            let denom = 2.0 * fp * fp - f * fpp;
            if denom == 0.0 || !denom.is_finite() {
                dx = 0.5 * (xh - xl);
                rts = xl + dx;
            } else {
                dx = (2.0 * f * fp) / denom;
                rts -= dx;
            }
            let rel = if rts != 0.0 { (dx / rts).abs() } else { dx.abs() };
            if rel < tol {
                return Ok(rts);
            }
        }

        let (ff, fp_new, fpp_new, _) = root_function(rts, mu, alpha, r0, vr0, dt);
        f = ff;
        fp = fp_new;
        fpp = fpp_new;

        if f < 0.0 {
            xl = rts;
            fl = f;
        } else {
            xh = rts;
            fh = f;
        }
    }

    Err(Box::new(io::Error::new(
        io::ErrorKind::Other,
        "solver did not converge within max_iter",
    )))
}
