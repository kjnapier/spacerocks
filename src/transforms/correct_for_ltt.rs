use crate::constants::*;
use crate::StateVector;
use crate::SpaceRock;
use crate::Observer;

/// Calculates the observer-centric state vector of a rock, accounting for light-time travel.
///
/// # Arguments
/// * `rock` - A SpaceRock object representing the rock.
/// * `observer` - An Observer object representing the observer.
///
/// # Returns
/// * A StateVector object representing the observer-centric state vector of the rock.
// pub fn correct_for_ltt(rock: &SpaceRock, observer: &Observer) -> StateVector {
//     // calculates the observer-centric state vector of a rock, accounting for light-time travel

//     let mut temp = StateVector::new(rock.position, rock.velocity);

//     let r = rock.position.norm();
//     let xi = MU_BARY / (r * r * r);    
//     let mut ltt0: f64;

//     let mut d_pos = temp.position - observer.position;
//     let mut delta = d_pos.norm();
//     let mut ltt = delta / SPEED_OF_LIGHT;
//     let mut acc = xi * ltt;

//     for _ in 0..10 {

//         ltt0 = ltt;
//         acc = xi * ltt;
//         temp.position = rock.position - (0.5 * acc * rock.position + rock.velocity) * ltt;
//         d_pos = temp.position - observer.position;
//         delta = d_pos.norm();
//         ltt = delta / SPEED_OF_LIGHT;
//         let dltt = (ltt - ltt0).abs();
        
//         // if dltt < 1.0e-6 {
//         //     break;
//         // }

//         if dltt < 1.0e-6 {
//             break;
//         }

//         // acc = xi * ltt;
//     }

//     temp.velocity = rock.velocity + acc * rock.position;
//     let d_vel = temp.velocity - observer.velocity.unwrap();


//     return StateVector::new(d_pos, d_vel);

// }

pub fn correct_for_ltt(rock: &SpaceRock, observer: &Observer) -> StateVector {
    let r0 = rock.position;
    let v0 = rock.velocity;
    let obs_pos = observer.position;
    let obs_vel = observer.velocity.expect("observer velocity required");

    let xi = MU_BARY / r0.norm().powi(3);
    let inv_c = 1.0 / SPEED_OF_LIGHT;

    let mut d_pos = r0 - obs_pos;
    let mut ltt = d_pos.norm() * inv_c;

    // iteration 1
    let mut acc = xi * ltt;
    let pos1 = r0 - v0 * ltt - r0 * (0.5 * acc * ltt);
    d_pos = pos1 - obs_pos;
    ltt = d_pos.norm() * inv_c;

    // iteration 2
    acc = xi * ltt;
    let pos2 = r0 - v0 * ltt - r0 * (0.5 * acc * ltt);
    d_pos = pos2 - obs_pos;
    ltt = d_pos.norm() * inv_c;

    acc = xi * ltt;
    let d_vel = (v0 - r0 * acc) - obs_vel;

    StateVector::new(d_pos, d_vel)
}

// pub fn correct_for_ltt(rock: &SpaceRock, observer: &Observer) -> StateVector {
//     // calculates the observer-centric state vector of a rock, accounting for light-time travel

//     let mut d_pos = rock.position - observer.position;
//     let mut delta = d_pos.norm();
//     let mut ltt = delta / SPEED_OF_LIGHT;
//     let mut ltt0 = ltt.clone();
//     let mut dltt = 1000.0;

//     let t0 = rock.epoch.clone();

//     let mut temp = rock.analytic_at(&(t0.clone() - ltt)).unwrap();

//     for _ in 0..3 {

//         ltt = delta / SPEED_OF_LIGHT;
//         dltt = (ltt - ltt0).abs();
        
//         if dltt < 1.0e-6 {
//             break;
//         }

//         ltt0 = ltt;
//         temp = rock.analytic_at(&(t0.clone() - ltt)).unwrap();
//         d_pos = temp.position - observer.position;
//         delta = d_pos.norm();

//         // acc = xi * ltt;
//     }

//     let d_vel = temp.velocity - observer.velocity.unwrap();
//     return StateVector::new(d_pos, d_vel);

// }