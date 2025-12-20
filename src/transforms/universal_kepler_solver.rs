use crate::transforms::stumpff::{stumpff_c, stumpff_s};

fn f(chi: f64, r0: f64, vr0: f64, alpha: f64, mu: f64, dt: f64) -> f64 {
    let z = alpha * chi * chi;
    let sqrt_mu = mu.sqrt();
    let c = stumpff_c(z);
    let s = stumpff_s(z);

    (r0 * vr0 / sqrt_mu) * chi * chi * c
        + (1.0 - alpha * r0) * chi * chi * chi * s
        + r0 * chi
        - sqrt_mu * dt
}

fn df_dchi(chi: f64, r0: f64, vr0: f64, alpha: f64, mu: f64) -> f64 {
    let z = alpha * chi * chi;
    let sqrt_mu = mu.sqrt();
    let c = stumpff_c(z);
    let s = stumpff_s(z);

    (r0 * vr0 / sqrt_mu) * chi * (1.0 - z * s)
        + (1.0 - alpha * r0) * chi * chi * c
        + r0
}

/// Safe derivatives of Stumpff S and C wrt z.
/// These are ONLY needed if you want Halley (second derivative).
/// They must handle z ~ 0 safely.
fn dS_dz(z: f64) -> f64 {
    // Series near 0:
    // S(z) = 1/6 - z/120 + z^2/5040 - ...
    // dS/dz = -1/120 + z/2520 - z^2/100800 + ...
    let az = z.abs();
    if az < 1e-8 {
        return -1.0 / 120.0 + z / 2520.0 - (z * z) / 100800.0;
    }
    // identity: dS/dz = (C(z) - 3S(z)) / (2z)
    (stumpff_c(z) - 3.0 * stumpff_s(z)) / (2.0 * z)
}

fn dC_dz(z: f64) -> f64 {
    // Series near 0:
    // C(z) = 1/2 - z/24 + z^2/720 - ...
    // dC/dz = -1/24 + z/360 - z^2/13440 + ...
    let az = z.abs();
    if az < 1e-8 {
        return -1.0 / 24.0 + z / 360.0 - (z * z) / 13440.0;
    }
    // identity: dC/dz = (1 - zS(z) - 2C(z)) / (2z)
    (1.0 - z * stumpff_s(z) - 2.0 * stumpff_c(z)) / (2.0 * z)
}

fn d2f_dchi2(chi: f64, r0: f64, vr0: f64, alpha: f64, mu: f64) -> f64 {
    // Differentiate df safely:
    // df = (r0*vr0/sqrt(mu))*chi*(1 - zS) + (1-alpha*r0)*chi^2*C + r0
    // with z = alpha*chi^2, dz/dchi = 2*alpha*chi
    let sqrt_mu = mu.sqrt();
    let z = alpha * chi * chi;
    let dz = 2.0 * alpha * chi;

    let s = stumpff_s(z);
    let c = stumpff_c(z);
    let ds = dS_dz(z);
    let dc = dC_dz(z);

    let a = r0 * vr0 / sqrt_mu;
    let b = 1.0 - alpha * r0;

    // d/dchi [ chi*(1 - zS) ] = (1 - zS) + chi * d/dchi(1 - zS)
    // d/dchi(1 - zS) = -(dz*S + z*dS*dz)
    let term1 = a * ((1.0 - z * s) + chi * (-(dz * s + z * ds * dz)));

    // d/dchi [ chi^2 * C(z) ] = 2chi*C + chi^2 * dC/dz * dz
    let term2 = b * (2.0 * chi * c + chi * chi * dc * dz);

    term1 + term2
}

/// Robust universal anomaly solve: bracket + safeguarded Halley/Newton + bisection fallback.
/// - Works for elliptical / parabolic-near / hyperbolic, as long as stumpff_* are stable.
pub fn solve_for_universal_anomaly(
    r0: f64,
    vr0: f64,
    alpha: f64,
    mu: f64,
    dt: f64,
    tol: f64,
    max_iter: usize,
) -> Result<f64, Box<dyn std::error::Error>> {
    if !r0.is_finite() || !vr0.is_finite() || !alpha.is_finite() || !mu.is_finite() || !dt.is_finite() {
        return Err("Non-finite input to universal anomaly solver".into());
    }
    if mu <= 0.0 {
        return Err("mu must be > 0".into());
    }
    if tol <= 0.0 {
        return Err("tol must be > 0".into());
    }

    // Handle dt = 0 exactly
    if dt == 0.0 {
        return Ok(0.0);
    }

    let sqrt_mu = mu.sqrt();

    // --- Initial guess (much safer) ---
    // Near-parabolic alpha~0 => chi ~ sqrt(mu)*dt/r0
    // Otherwise chi ~ sqrt(mu)*|alpha|*dt is often okay but can be too aggressive;
    // so blend toward the parabolic guess.
    let chi_par = sqrt_mu * dt.abs() / r0.max(1e-15);
    let chi_lin = sqrt_mu * dt.abs() * alpha.abs();
    let mut chi0 = if alpha.abs() < 1e-6 {
        chi_par
    } else {
        // soft blend: prevents chi0 from being absurdly tiny or huge
        (0.5 * chi_par + 0.5 * chi_lin).max(1e-12)
    };

    // Give chi the sign of dt (common convention; f is odd-ish in chi for many cases)
    if dt < 0.0 {
        chi0 = -chi0;
    }

    // --- Bracket the root ---
    // For dt>0: f(0) = -sqrt(mu)*dt < 0. We want f(chi_hi) > 0.
    // For dt<0: f(0) > 0. We want f(chi_hi) < 0, with chi_hi negative.
    let mut chi_lo = 0.0;
    let mut f_lo = f(chi_lo, r0, vr0, alpha, mu, dt);

    let mut chi_hi = chi0;
    let mut f_hi = f(chi_hi, r0, vr0, alpha, mu, dt);

    // Expand until sign change or we give up.
    // Expansion factor 2 is usually plenty.
    let mut expand = 0;
    while f_lo.signum() == f_hi.signum() && expand < 60 {
        chi_hi *= 2.0;
        f_hi = f(chi_hi, r0, vr0, alpha, mu, dt);
        expand += 1;
    }
    if f_lo.signum() == f_hi.signum() {
        // If we cannot bracket, something is numerically wrong (often Stumpff instability),
        // or dt is absurdly large for the units.
        return Err("Failed to bracket universal Kepler root (check stumpff stability / units)".into());
    }

    // Ensure ordering: chi_lo < chi_hi
    if chi_lo > chi_hi {
        std::mem::swap(&mut chi_lo, &mut chi_hi);
        std::mem::swap(&mut f_lo, &mut f_hi);
    }

    // Start from mid or chi0 clamped to bracket
    let mut chi = chi0.clamp(chi_lo, chi_hi);
    let mut f_chi = f(chi, r0, vr0, alpha, mu, dt);

    // --- Iteration ---
    for _iter in 0..max_iter {
        let err = f_chi.abs();
        if err <= tol {
            return Ok(chi);
        }

        let fp = df_dchi(chi, r0, vr0, alpha, mu);

        // Default: bisection step
        let mut chi_new = 0.5 * (chi_lo + chi_hi);

        // Try Halley/Newton step if derivative is usable
        if fp.is_finite() && fp.abs() > 1e-16 {
            let fpp = d2f_dchi2(chi, r0, vr0, alpha, mu);
            let mut step = f_chi / fp; // Newton

            // Halley correction when fpp is sane
            if fpp.is_finite() {
                let denom = 1.0 - 0.5 * step * (fpp / fp);
                if denom.is_finite() && denom.abs() > 1e-12 {
                    step /= denom;
                }
            }

            // Trust region: prevent ridiculous jumps
            let max_step = 0.5 * (chi_hi - chi_lo);
            if step.is_finite() {
                step = step.clamp(-max_step, max_step);
                let candidate = chi - step;

                // Accept candidate only if it stays in bracket
                if candidate > chi_lo && candidate < chi_hi {
                    chi_new = candidate;
                }
            }
        }

        let f_new = f(chi_new, r0, vr0, alpha, mu, dt);

        // Update bracket
        if f_lo.signum() == f_new.signum() {
            chi_lo = chi_new;
            f_lo = f_new;
        } else {
            chi_hi = chi_new;
            f_hi = f_new;
        }

        chi = chi_new;
        f_chi = f_new;
    }

    Err("Universal Kepler solver did not converge within max_iter".into())
}



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
// pub fn solve_for_universal_anomaly(r0: f64, vr0: f64, alpha: f64, mu: f64, dt: f64, tol: f64, max_iter: usize) -> Result<f64, Box<dyn std::error::Error>> {
//     let mut chi = mu.sqrt() * alpha.abs() * dt;
//     let mut iter = 0;
//     let mut error = f(chi, r0, vr0, alpha, mu, dt).abs();
//     let n = 3.0; // Polynomial degree

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
        
//         let G = df_dchi(chi, r0, vr0, alpha, mu) / error; 
//         let H = G.powi(2) - d2f_dchi2(chi, r0, vr0, alpha, mu) / error;

//         let denom_plus = G + ((n-1.0) * (n*H - G.powi(2))).sqrt();
//         let denom_minus = G - ((n-1.0) * (n*H - G.powi(2))).sqrt();
//         let denom = if denom_plus.abs() > denom_minus.abs() { denom_plus } else { denom_minus };

//         let a = n / denom;
//         chi = chi - a;
//         error = f(chi, r0, vr0, alpha, mu, dt).abs();
//         iter += 1;
//     }
   
//     Ok(chi)
// }

