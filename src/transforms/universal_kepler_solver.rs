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


pub fn newton_raphson(chi0: f64, r0: f64, vr0: f64, alpha: f64, mu: f64, dt: f64, tol: f64, max_iter: usize) -> f64 {
    let mut chi = chi0;
    let mut iter = 0;
    let mut error = f(chi, r0, vr0, alpha, mu, dt).abs();
    while error > tol && iter < max_iter {
        let delta_chi = f(chi, r0, vr0, alpha, mu, dt) / df_dchi(chi, r0, vr0, alpha, mu);
        chi -= delta_chi;
        error = f(chi, r0, vr0, alpha, mu, dt).abs();
        iter += 1;
    }
    chi
}

fn lagrange_f_and_g(chi: f64, r0: f64, vr0: f64, alpha: f64, mu: f64, dt: f64) -> (f64, f64) {
    let z = alpha * chi.powi(2);
    let f = 1.0 - chi.powi(2) / r0 * stumpff_c(z);
    let g = dt - chi.powi(3) / mu.sqrt() * stumpff_s(z);
    (f, g)
}