#[allow(non_snake_case)]
pub fn calc_M_from_E(e: f64, E: f64) -> f64 {
    let mut M;

    if e < 1.0 {
        M = E - e * E.sin();
    }

    else {
        M = e * E.sinh() - E;
    }

    return M;
}