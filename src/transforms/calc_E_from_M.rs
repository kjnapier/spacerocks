#![allow(non_snake_case)]
pub fn calc_E_from_M(e: f64, M: f64) -> f64 {

    let mut E;

    // Handle trivial cases
    if e == 0.0 {
        E = M;
        return E;
    }

    if M == 0.0 {
        E = 0.0;
        return E;
    }
    
    if e < 1.0 {

        // Define initial estimate
        let sinM = M.sin();
        if sinM < 0.0 {
            E = M - e * sinM;
        }
        else {
            E = M + e * sinM;
        }

        // Perform Newton-Raphson estimate
        for _ in 0..10 {

            // Compute f(E), f'(E), f''(E) and f'''(E), avoiding recomputation of sine and cosine.
            let esinE = e * E.sin();
            let ecosE = e * E.cos();
            
            let f_E = E - esinE - M;

            if f_E.abs() < 1.0e-16 {
                return E;
            }

            let fP_E = 1.0 - ecosE;
            let fPP_E = esinE;
            let fPPP_E = ecosE;

            let delta_i1 = -f_E / fP_E;
            let delta_i2 = -f_E / (fP_E + 0.5 * delta_i1 * fPP_E);
            let delta_i3 = -f_E / (fP_E + 0.5 * delta_i2 * fPP_E + 1.0/6.0 * fPPP_E * delta_i2 * delta_i2);
            
            // Update E
            E += delta_i3;
        }
    }

    else {

        E = M / M.abs() * (2.0 * M.abs() / e + 1.8).ln();
  		let mut F = E - e * E.sinh() + M;
  		for _ in 1..100 {

  			E = E - F / (1.0 - e * E.cosh());
  			F = E - e * E.sinh() + M;
  			if F.abs() < 1.0e-16 {
  				break;
        }
      }
    }
    return E;
}



// pub fn calc_E_from_M(e: f64, M: f64) -> f64 {

//     let mut E;

//     // Handle trivial cases
//     if e == 0.0 {
//         E = M;
//         return E;
//     }

//     if M == 0.0 {
//         E = 0.0;
//         return E;
//     }
    
//     if e < 1.0 {

//         let mut flag = false;
//         let mut M = M;
//         if M > std::f64::consts::PI {
//             M = 2.0 * std::f64::consts::PI - M;
//             flag = true;
//         }

//         // Define initial estimate
//         let sinM = M.sin();
//         E = e * sinM + f64::max(M, e * (sinM + 0.591));
        

//         // Perform Newton-Raphson estimate
//         for _ in 0..10 {

//             // Compute f(E), f'(E), f''(E) and f'''(E), avoiding recomputation of sine and cosine.
//             let esinE = e * E.sin();
//             let ecosE = e * E.cos();
            
//             let f_E = E - esinE - M;

//             if f_E.abs() < 1.0e-15 {
//                 if flag {
//                     E = 2.0 * std::f64::consts::PI - E;
//                 }
//                 return E;
//             }

//             let fP_E = 1.0 - ecosE;
//             let fPP_E = esinE;
//             let fPPP_E = ecosE;

//             let delta_i1 = -f_E / fP_E;
//             let delta_i2 = -f_E / (fP_E + 0.5 * delta_i1 * fPP_E);
//             let delta_i3 = -f_E / (fP_E + 0.5 * delta_i2 * fPP_E + 1.0/6.0 * fPPP_E * delta_i2 * delta_i2);
            
//             // Update E
//             E += delta_i3;
//         }
//     }

//     else {

//         E = M / M.abs() * (2.0 * M.abs() / e + 1.8).ln();
//   		let mut F = E - e * E.sinh() + M;
//   		for _ in 1..100 {

//   			E = E - F / (1.0 - e * E.cosh());
//   			F = E - e * E.sinh() + M;
//   			if F.abs() < 1.0e-16 {
//   				break;
//         }
//       }
//     }
//     return E;
// }


// imlement test for calc_E_from_M
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_calc_E_from_M() {
        let e = 0.0;
        let M = 0.1;
        let E = calc_E_from_M(e, M);
        assert_eq!(E, 0.1);
    }

   
}