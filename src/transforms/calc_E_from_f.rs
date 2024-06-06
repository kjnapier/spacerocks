
#![allow(non_snake_case)]
pub fn calc_E_from_f(e: f64, f: f64) -> f64 {
    let mut E;

    if e < 1.0 {
        E = 2.0 * ((1.0 - e).sqrt() * (f / 2.0).sin()).atan2((1.0 + e).sqrt() * (f / 2.0).cos())
    }

    else {
        // let cta = f.cos();
        // E = ((cta + e) / (1.0 + e * cta)).acosh();
        // if f < 0.0 {
        //     E *= -1.0;
        // }

        E = 2.0 * (((e - 1.0) / (e + 1.0)).sqrt() * (f / 2.0).tan()).atanh();
        


    }
    return E;
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_calc_E_from_f() {
        let e = 0.0;
        let f = 0.0;
        let result = calc_E_from_f(e, f);
        assert_eq!(result, 0.0);
    }
}
