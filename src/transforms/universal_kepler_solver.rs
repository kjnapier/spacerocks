use crate::transforms::stumpff::{stumpff_c, stumpff_s};

fn f(chi: f64, r0: f64, vr0: f64, alpha: f64, mu: f64, dt: f64) -> f64 {
    let z = alpha * chi.powi(2);
    let first_term = r0 * vr0 / mu.sqrt() * chi.powi(2) * stumpff_c(z);
    let second_term = (1.0 - alpha * r0) * chi.powi(3) * stumpff_s(z);
    let third_term = r0 * chi;
    let fourth_term = dt * mu.sqrt();
    first_term + second_term + third_term - fourth_term
}

fn df_dchi(chi: f64, r0: f64, vr0: f64, alpha: f64, mu: f64) -> f64 {
    let z = alpha * chi.powi(2);
    let first_term = r0 * vr0 / mu.sqrt() * chi * (1.0 - z * stumpff_s(z));
    let second_term = (1.0 - alpha * r0) * chi.powi(2) * stumpff_c(z);
    let third_term = r0;
    first_term + second_term + third_term
}

pub fn solve_for_universal_anomaly(r0: f64, vr0: f64, alpha: f64, mu: f64, dt: f64, tol: f64, max_iter: usize) -> Result<f64, Box<dyn std::error::Error>> {
    let mut chi = mu.sqrt() * alpha.abs() * dt;
    let mut iter = 0;
    let mut error = f(chi, r0, vr0, alpha, mu, dt).abs();

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
        let f_val = f(chi, r0, vr0, alpha, mu, dt);
        let df_val = df_dchi(chi, r0, vr0, alpha, mu);
        let delta_chi = f_val / df_val;
        chi -= delta_chi;
        error = f(chi, r0, vr0, alpha, mu, dt).abs();
        iter += 1;
    }
    Ok(chi)
}