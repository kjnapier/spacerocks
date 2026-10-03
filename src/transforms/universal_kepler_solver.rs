use crate::transforms::stumpff::{stumpff_c, stumpff_s};

use nalgebra::Vector3;

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
fn ds_dz(z: f64) -> f64 {
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

fn dc_dz(z: f64) -> f64 {
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
    let ds = ds_dz(z);
    let dc = dc_dz(z);

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
        }

        chi = chi_new;
        f_chi = f_new;
    }

    Err("Universal Kepler solver did not converge within max_iter".into())
}


/// Stumpff functions c0..c3 at `z`, by series at a reduced argument and the
/// quadruple-argument identities (no trigonometric calls).
#[inline]
fn stumpff_c0123(z: f64) -> [f64; 4] {
    let mut z = z;
    let mut n = 0;
    while z.abs() > 0.1 {
        z *= 0.25;
        n += 1;
    }
    // c2 = sum (-z)^j / (2j+2)!, c3 = sum (-z)^j / (2j+3)!, to z^6 (error < 1e-21 at |z| = 0.1)
    let c2 = 1.0 / 2.0 - z * (1.0 / 24.0 - z * (1.0 / 720.0 - z * (1.0 / 40320.0 - z * (1.0 / 3628800.0 - z * (1.0 / 479001600.0 - z / 87178291200.0)))));
    let c3 = 1.0 / 6.0 - z * (1.0 / 120.0 - z * (1.0 / 5040.0 - z * (1.0 / 362880.0 - z * (1.0 / 39916800.0 - z * (1.0 / 6227020800.0 - z / 1307674368000.0)))));
    let (mut c0, mut c1, mut c2, mut c3) = (1.0 - z * c2, 1.0 - z * c3, c2, c3);
    for _ in 0..n {
        c3 = (c2 + c0 * c3) * 0.25;
        c2 = 0.5 * c1 * c1;
        c1 *= c0;
        c0 = 2.0 * c0 * c0 - 1.0;
    }
    [c0, c1, c2, c3]
}

/// Advance a two-body orbit by `dt` from position `r0` and velocity `v0` about a body of
/// gravitational parameter `mu`, with Gauss's f and g functions in universal variables
/// (Stumpff-function form, as in Wisdom & Hernandez 2015). Works for any conic and any `dt`.
///
/// Kepler's equation is solved with Halley's method from a starting guess that fits the step:
/// the third-order Taylor series of the universal anomaly for short steps, Danby's starter in
/// the eccentric anomaly for long steps on bound orbits. Halley's method converges cubically,
/// so once its error estimate for the last correction is below round-off the solver stops,
/// carrying the Stumpff functions to the final anomaly by their Taylor series instead of
/// evaluating them again. If the iteration does not converge it falls back to the bracketing
/// [`solve_for_universal_anomaly`].
///
/// # Returns
/// * The position and velocity after `dt`.
pub fn universal_kepler_step(r0: &Vector3<f64>, v0: &Vector3<f64>, mu: f64, dt: f64) -> Result<(Vector3<f64>, Vector3<f64>), Box<dyn std::error::Error>> {
    let r0n = r0.norm();
    let eta0 = r0.dot(v0);
    let beta = 2.0 * mu / r0n - v0.norm_squared();
    let zeta0 = mu - beta * r0n;

    // Whole periods of a bound orbit change nothing. (The period is 2 pi mu / beta^1.5.)
    let mut dt_red = dt;
    let tau_mu = std::f64::consts::TAU * mu;
    if beta > 0.0 && dt * dt * beta * beta * beta > tau_mu * tau_mu {
        let period = tau_mu / (beta * beta.sqrt());
        dt_red -= period * (dt / period).trunc();
    }

    // Initial guess. When the body moves a small fraction of its distance, the Taylor series of
    // x(t) (dx/dt = 1/r) to third order. Otherwise, on a bound orbit, Danby's starter for
    // Kepler's equation in the eccentric anomaly (x = dE / sqrt(beta)), and on an unbound one
    // the series to second order.
    let (t, ir) = (dt_red, 1.0 / r0n);
    let mut x = if t * t * v0.norm_squared() < 0.09 * r0n * r0n {
        t * ir * (1.0 - t * ir * ir * (0.5 * eta0 - t * ir * (0.5 * eta0 * eta0 * ir - zeta0 / 6.0)))
    } else if beta > 0.0 {
        let sb = beta.sqrt();
        // e cos E and e sin E at the start, and the mean anomaly at the end.
        let (ecos, esin) = (1.0 - r0n * beta / mu, eta0 * sb / mu);
        let e0 = esin.atan2(ecos);
        let n_dt = beta * sb / mu * dt_red;
        let m = e0 - esin + n_dt;
        let m = m - std::f64::consts::TAU * (m / std::f64::consts::TAU).round();
        let e1 = m + 0.85 * ecos.hypot(esin) * m.sin().signum();
        // E changes by n dt plus less than 2e, so pick the branch of E1 - E0 closest to n dt.
        let mut de = e1 - e0;
        de += std::f64::consts::TAU * ((n_dt - de) / std::f64::consts::TAU).round();
        de / sb
    } else {
        t * ir - 0.5 * eta0 * t * t * ir * ir * ir
    };

    let mut converged = false;
    let mut last_dx = f64::INFINITY;
    let (mut g1, mut g2, mut g3) = (0.0, 0.0, 0.0);
    for _ in 0..30 {
        let [c0, c1, c2, c3] = stumpff_c0123(beta * x * x);
        let g0 = c0;
        (g1, g2, g3) = (x * c1, x * x * c2, x * x * x * c3);
        let f = r0n * x + eta0 * g2 + zeta0 * g3 - dt_red;
        let fp = r0n + eta0 * g1 + zeta0 * g2;
        let fpp = eta0 * g0 + zeta0 * g1;
        let ifp = 1.0 / fp;
        let dx = f / (fp - 0.5 * f * fpp * ifp);
        x -= dx;
        // Halley's error after this correction is about K dx^3. Once that is far below
        // round-off, carry the g-functions to the new x (g_k' = g_(k-1), g_0' = -beta g_1).
        if dx.abs() <= 1e-4 * x.abs() {
            let fppp = zeta0 * g0 - beta * eta0 * g1;
            let k = (ifp * (0.25 * fpp * fpp * ifp - fppp / 6.0)).abs();
            if 10.0 * k * dx.abs().powi(3) <= 1e-17 * x.abs() {
                let d = -dx;
                let (d2, d3) = (0.5 * d * d, d * d * d / 6.0);
                (g1, g2, g3) = (
                    g1 + d * g0 - d2 * beta * g1 - d3 * beta * g0,
                    g2 + d * g1 + d2 * g0 - d3 * beta * g1,
                    g3 + d * g2 + d2 * g1 + d3 * g0,
                );
                converged = true;
                break;
            }
        }
        if dx.abs() <= 2e-16 * x.abs() || (dx.abs() >= last_dx && dx.abs() <= 1e-12 * x.abs()) {
            let [_, c1, c2, c3] = stumpff_c0123(beta * x * x);
            (g1, g2, g3) = (x * c1, x * x * c2, x * x * x * c3);
            converged = true;
            break;
        }
        last_dx = dx.abs();
    }

    if converged && x.is_finite() {
        let irn = 1.0 / (r0n + eta0 * g1 + zeta0 * g2);
        let f = 1.0 - mu * ir * g2;
        let g = dt_red - mu * g3;
        let fdot = -mu * ir * irn * g1;
        let gdot = 1.0 - mu * irn * g2;
        let (r, v) = (f * r0 + g * v0, fdot * r0 + gdot * v0);
        if r.iter().chain(v.iter()).all(|c| c.is_finite()) {
            return Ok((r, v));
        }
    }
    kepler_step_bracketing(r0, v0, mu, dt)
}

/// [`universal_kepler_step`] with the robust bracketing universal-anomaly solver.
fn kepler_step_bracketing(r0: &Vector3<f64>, v0: &Vector3<f64>, mu: f64, dt: f64) -> Result<(Vector3<f64>, Vector3<f64>), Box<dyn std::error::Error>> {
    let r0n = r0.norm();
    let vr0 = r0.dot(v0) / r0n;
    let alpha = 2.0 / r0n - v0.norm_squared() / mu;
    let sqrt_mu = mu.sqrt();
    let tol = 1e-15 * sqrt_mu * dt.abs();
    let chi = solve_for_universal_anomaly(r0n, vr0, alpha, mu, dt, tol, 100)
        .or_else(|_| solve_for_universal_anomaly(r0n, vr0, alpha, mu, dt, 1e3 * tol, 200))?;

    let z = alpha * chi * chi;
    let (c, s) = (stumpff_c(z), stumpff_s(z));
    let f = 1.0 - chi * chi / r0n * c;
    let g = dt - chi * chi * chi / sqrt_mu * s;
    let r = f * r0 + g * v0;
    let rn = r.norm();
    let fdot = sqrt_mu / (rn * r0n) * chi * (z * s - 1.0);
    let gdot = 1.0 - chi * chi / rn * c;
    Ok((r, fdot * r0 + gdot * v0))
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


#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::GRAVITATIONAL_CONSTANT;

    /// States on random orbits of several families, with a step for each: (a range, e range, dt).
    fn orbits(n_per_family: usize) -> Vec<(Vector3<f64>, Vector3<f64>, f64)> {
        // A small LCG, so the test needs no RNG.
        let mut seed: u64 = 1;
        let mut rnd = || {
            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            (seed >> 11) as f64 / (1u64 << 53) as f64
        };
        let mu = GRAVITATIONAL_CONSTANT;
        let families = [
            (38.0, 50.0, 0.0, 0.2, 120.0),   // KBOs
            (5.0, 30.0, 0.0, 0.06, 120.0),   // giants
            (0.8, 3.0, 0.0, 0.6, 5.0),       // NEOs and main belt
            (2.0, 50.0, 0.9, 0.999, 20.0),   // comets
            (-5.0, -1.0, 1.1, 3.0, 50.0),    // hyperbolic
            (1.0, 3.0, 0.0, 0.3, 3000.0),    // many periods
        ];
        let mut out = vec![];
        for (a0, a1, e0, e1, dt) in families {
            for _ in 0..n_per_family {
                let a: f64 = a0 + (a1 - a0) * rnd();
                let e: f64 = e0 + (e1 - e0) * rnd();
                let p = a * (1.0 - e * e);
                let nu: f64 = if e < 1.0 { std::f64::consts::TAU * rnd() } else { 0.95 * (2.0 * rnd() - 1.0) * (-1.0 / e).acos() };
                let r = p / (1.0 + e * nu.cos());
                let k = (mu / p).sqrt();
                let sign = if rnd() < 0.5 { -1.0 } else { 1.0 };
                out.push((Vector3::new(r * nu.cos(), r * nu.sin(), 0.0), Vector3::new(-k * nu.sin(), k * (e + nu.cos()), 0.0), sign * dt * (0.5 + rnd())));
            }
        }
        out
    }

    #[test]
    fn kepler_step_matches_bracketing_solver() {
        let mu = GRAVITATIONAL_CONSTANT;
        let cases = [
            (Vector3::new(1.0, 0.1, 0.0), Vector3::new(0.001, 0.017, 0.002), 30.0),   // near circular
            (Vector3::new(0.05, 0.0, 0.0), Vector3::new(0.0, 0.1, 0.0), 5.0),         // e ~ 0.98 at pericenter
            (Vector3::new(-4.9, 0.3, 0.0), Vector3::new(0.0, -0.0006, 0.0), 2000.0),  // near aphelion, long step
            (Vector3::new(1.0, 0.0, 0.0), Vector3::new(0.0, 0.03, 0.0), 400.0),       // hyperbolic
            (Vector3::new(1.0, 0.0, 0.0), Vector3::new(0.0, 0.0172, 0.0), 3650.25),   // many periods
            (Vector3::new(1.0, 0.0, 0.0), Vector3::new(0.0, 0.0172, 0.0), -123.4),    // backwards
            // hyperbolic near pericenter with a long step
            (Vector3::new(0.12550094620688393, 0.027550727404811702, 0.0), Vector3::new(-0.0071336265344599335, 0.06930541831091576, 0.0), 54.536383668498914),
        ];
        let orbits = orbits(500);
        for (r0, v0, dt) in cases.into_iter().chain(orbits) {
            let (r, v) = universal_kepler_step(&r0, &v0, mu, dt).unwrap();
            let (rb, vb) = kepler_step_bracketing(&r0, &v0, mu, dt).unwrap();
            assert!((r - rb).norm() < 1e-10 * rb.norm(), "r {r:?} vs {rb:?} ({r0:?}, {v0:?}, {dt})");
            assert!((v - vb).norm() < 1e-10 * vb.norm(), "v {v:?} vs {vb:?} ({r0:?}, {v0:?}, {dt})");
        }
    }

    #[test]
    fn kepler_step_is_reversible() {
        let mu = GRAVITATIONAL_CONSTANT;
        for (r0, v0, dt) in orbits(200) {
            let (r1, v1) = universal_kepler_step(&r0, &v0, mu, dt).unwrap();
            let (r2, v2) = universal_kepler_step(&r1, &v1, mu, -dt).unwrap();
            assert!((r2 - r0).norm() < 1e-11 * r0.norm(), "{r0:?} {v0:?} {dt}: back at {r2:?}");
            assert!((v2 - v0).norm() < 1e-11 * v0.norm(), "{r0:?} {v0:?} {dt}: back at {v2:?}");
        }
    }
}
