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
pub fn correct_for_ltt(rock: &SpaceRock, observer: &Observer) -> StateVector {
    // calculates the observer-centric state vector of a rock, accounting for light-time travel

    let mut temp = StateVector::new(rock.position, rock.velocity);

    let r = rock.position.norm();
    let xi = MU_BARY / (r * r * r);    
    let mut ltt0: f64;

    let mut d_pos = temp.position - observer.position();
    let mut delta = d_pos.norm();
    let mut ltt = delta / SPEED_OF_LIGHT;
    let mut acc = xi * ltt;

    for _ in 0..10 {

        ltt0 = ltt;
        acc = xi * ltt;
        temp.position = rock.position - (0.5 * acc * rock.position + rock.velocity) * ltt;
        d_pos = temp.position - observer.position();
        delta = d_pos.norm();
        ltt = delta / SPEED_OF_LIGHT;
        let dltt = (ltt - ltt0).abs();
        
        // if dltt < 1.0e-6 {
        //     break;
        // }

        if dltt < 1.0e-10 {
            break;
        }

        // acc = xi * ltt;
    }

    temp.velocity = rock.velocity + acc * rock.position;
    let d_vel = temp.velocity - observer.velocity();


    return StateVector::new(d_pos, d_vel);

}