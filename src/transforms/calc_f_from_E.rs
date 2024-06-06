#[allow(non_snake_case)]
pub fn calc_f_from_E(e: f64, E: f64) -> f64 {
    let mut f;

    if e < 1.0 {
        f = 2.0 * ((1.0 + e).sqrt() * (E / 2.0).sin()).atan2((1.0 - e).sqrt() * (E / 2.0).cos())
    }

    else {
        //f = 2.0 * atan2(sqrt((e + 1)/(e - 1)) * tanh(E / 2), 1);
        f = 2.0 * (((e + 1.0) / (e - 1.0)).sqrt() * (E / 2.0).tanh()).atan2(1.0);
        // f = 2.0 * (E.cosh().acos() - e).atan2((1.0 + e) * (E.cosh()).sqrt());
    }

    return f;
}