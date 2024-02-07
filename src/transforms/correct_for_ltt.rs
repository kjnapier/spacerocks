use crate::constants::*;
use crate::statevector::StateVector;
use crate::spacerock::SpaceRock;

pub fn correct_for_ltt(rock: &SpaceRock, observer: &SpaceRock) -> StateVector {

    let mut temp = StateVector::new(rock.position.x, rock.position.y, rock.position.z, 
                                    rock.velocity.x, rock.velocity.y, rock.velocity.z);

    let r = rock.position.norm();
    let xi = MU_BARY / (r * r * r);    
    let mut ltt0 = 0.0;

    let mut d_pos = temp.position - observer.position;
    let mut delta = d_pos.norm();
    let mut ltt = delta / SPEED_OF_LIGHT;
    let mut acc = xi * ltt;

    for _ in 0..3 {

        ltt0 = ltt;
        temp.position = rock.position - (0.5 * acc * rock.position + rock.velocity) * ltt;
        d_pos = temp.position - observer.position;
        delta = d_pos.norm();
        ltt = delta / SPEED_OF_LIGHT;
        let dltt = (ltt - ltt0).abs();
        
        if dltt < 1.0e-6 {
            break;
        }

        acc = xi * ltt;
    }

    temp.velocity = rock.velocity + acc * rock.position;
    let d_vel = temp.velocity - observer.velocity;

    return StateVector::new(d_pos.x, d_pos.y, d_pos.z, d_vel.x, d_vel.y, d_vel.z);

}