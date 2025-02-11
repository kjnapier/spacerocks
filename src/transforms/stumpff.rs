//! Stumpff functions for the universal variable formulation.

/// Calculate the Stumpff S function.
///
/// # Arguments
/// * `z` - A real number.
///
/// # Returns
/// * The Stumpff S function evaluated at `z`.
pub fn stumpff_s(z: f64) -> f64 {
    if z == 0.0 {
        1.0 / 6.0
    } else if z > 0.0 {
        let rootz = z.sqrt();
        (rootz - rootz.sin()) / rootz.powi(3)
    } else {
        let rootz = (-z).sqrt();
        (rootz.sinh() - rootz) / rootz.powi(3)
    }
}

/// Calculate the Stumpff C function.
///
/// # Arguments
/// * `z` - A real number.
///
/// # Returns
/// * The Stumpff C function evaluated at `z`.
pub fn stumpff_c(z: f64) -> f64 {
    if z == 0.0 {
        1.0 / 2.0
    } else if z > 0.0 {
        let rootz = z.sqrt();
        (1.0 - rootz.cos()) / z
    } else {
        let rootz = (-z).sqrt();
        (rootz.cosh() - 1.0) / z
    }
}