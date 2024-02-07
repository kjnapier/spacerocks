
#[derive(Debug, Clone, PartialEq)]
pub struct SphericalState {
    pub phi: f64,
    pub theta: f64,
    pub r: f64,
    pub vr: f64,
    pub vo: f64,
    pub inc: f64
}

impl SphericalState {
    pub fn new(phi: f64, theta: f64, r: f64, vr: f64, vo: f64, inc: f64) -> Self {
        SphericalState {
            phi: phi,
            theta: theta,
            r: r,
            vr: vr,
            vo: vo,
            inc: inc
        }
    }
}