//! Physics on plain state vectors.
//!
//! The functions here take a Cartesian state `[x, y, z, vx, vy, vz]` (AU, AU/day) and plain
//! numbers (μ in AU³/day², times in days) and know nothing about names, epochs, origins or
//! reference planes. Keeping that metadata, and checking it is consistent, is the job of the
//! types built on top ([`SpaceRock`](crate::SpaceRock), [`Population`](crate::Population)),
//! whose methods are thin wrappers around these functions. A collection applies the same
//! function over its state array, so single-rock and batch code share one implementation.

use nalgebra::{Matrix3, Vector3};

use crate::transforms::{calc_conic_anomaly_from_true_anomaly, calc_mean_anomaly_from_conic_anomaly, universal_kepler_step};
use crate::ReferencePlane;

/// Cartesian state `[x, y, z, vx, vy, vz]` in AU and AU/day.
pub type State = [f64; 6];

/// Build a state from a position and a velocity.
#[inline]
pub fn from_pv(position: &Vector3<f64>, velocity: &Vector3<f64>) -> State {
    [position.x, position.y, position.z, velocity.x, velocity.y, velocity.z]
}

#[inline]
pub fn position(s: &State) -> Vector3<f64> {
    Vector3::new(s[0], s[1], s[2])
}

#[inline]
pub fn velocity(s: &State) -> Vector3<f64> {
    Vector3::new(s[3], s[4], s[5])
}

/// Split a state into position and velocity.
#[inline]
pub fn pv(s: &State) -> (Vector3<f64>, Vector3<f64>) {
    (position(s), velocity(s))
}

/// Move a state along its two-body orbit about a body with gravitational parameter `mu` by
/// `dt` days (universal-variable Kepler solver).
pub fn kepler_step(s: &State, mu: f64, dt: f64) -> Result<State, Box<dyn std::error::Error>> {
    let (r, v) = pv(s);
    let (r1, v1) = universal_kepler_step(&r, &v, mu, dt)?;
    Ok(from_pv(&r1, &v1))
}

/// Matrix rotating vectors from reference plane `from` to `to`.
pub fn rotation(from: &ReferencePlane, to: &ReferencePlane) -> Result<Matrix3<f64>, Box<dyn std::error::Error>> {
    // `from`'s inverse takes us back to J2000. Not every plane matrix is orthonormal to machine
    // precision, so use the inverse rather than the transpose.
    let inv = from.get_rotation_matrix().try_inverse().ok_or("Could not invert rotation matrix")?;
    Ok(to.get_rotation_matrix() * inv)
}

/// Apply a rotation to a state's position and velocity.
#[inline]
pub fn rotate(s: &State, rot: &Matrix3<f64>) -> State {
    let (r, v) = pv(s);
    from_pv(&(rot * r), &(rot * v))
}

/// State from Keplerian elements: perihelion distance `q` (AU), eccentricity `e`, and the
/// angles `inc`, `arg`, `node` and true anomaly `f` (radians), about a body with
/// gravitational parameter `mu`.
pub fn from_kepler(q: f64, e: f64, inc: f64, arg: f64, node: f64, f: f64, mu: f64) -> Result<State, Box<dyn std::error::Error>> {
    if e >= 1.0 {
        let max_true_anomaly = (-1.0 / e).acos();
        if f.abs() > max_true_anomaly {
            return Err("True anomaly is not commensurate with eccentricity".into());
        }
    }

    let p = q * (1.0 + e);
    let h = (p * mu).sqrt();
    let r = p / (1.0 + e * f.cos());
    let vr = mu * f.sin() * e / h;

    let rot_x = node.cos() * (arg + f).cos() - node.sin() * (arg + f).sin() * inc.cos();
    let rot_y = node.sin() * (arg + f).cos() + node.cos() * (arg + f).sin() * inc.cos();
    let rot_z = (arg + f).sin() * inc.sin();

    let rot_x2 = node.cos() * (arg + f).sin() + node.sin() * (arg + f).cos() * inc.cos();
    let rot_y2 = node.sin() * (arg + f).sin() - node.cos() * (arg + f).cos() * inc.cos();
    let rot_z2 = (arg + f).cos() * inc.sin();

    let nudot = h / r.powi(2);
    Ok([
        r * rot_x,
        r * rot_y,
        r * rot_z,
        vr * rot_x - r * nudot * rot_x2,
        vr * rot_y - r * nudot * rot_y2,
        vr * rot_z + r * nudot * rot_z2,
    ])
}

/// Osculating orbital elements of one state. Angles are in radians, distances in AU.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Elements {
    /// Semi-major axis (negative for unbound orbits).
    pub a: f64,
    pub e: f64,
    /// Perihelion distance.
    pub q: f64,
    pub inc: f64,
    /// Longitude of the ascending node (0 for equatorial orbits).
    pub node: f64,
    /// Argument of perihelion (0 for circular orbits).
    pub arg: f64,
    pub true_anomaly: f64,
    /// Eccentric, parabolic or hyperbolic anomaly (NaN if it cannot be computed).
    pub conic_anomaly: f64,
    /// NaN if it cannot be computed.
    pub mean_anomaly: f64,
}

/// The quantities every element is derived from, computed once.
struct Orbit {
    r: Vector3<f64>,
    hvec: Vector3<f64>,
    evec: Vector3<f64>,
    nvec: Vector3<f64>,
}

impl Orbit {
    #[inline]
    fn new(s: &State, mu: f64) -> Orbit {
        let (r, v) = pv(s);
        let hvec = r.cross(&v);
        let evec = v.cross(&hvec) / mu - r / r.norm();
        let nvec = Vector3::new(-hvec.y, hvec.x, 0.0);
        Orbit { r, hvec, evec, nvec }
    }

    /// The ascending node, or the x-axis for (near-)equatorial orbits, where the node is
    /// undefined and taken to be 0.
    #[inline]
    fn node_direction(&self) -> Vector3<f64> {
        let h = self.hvec.norm();
        if self.nvec.norm() <= 1e-11 * h {
            Vector3::new(1.0, 0.0, 0.0)
        } else {
            self.nvec / self.nvec.norm()
        }
    }

    /// Signed angle from `a` to `b` about the orbit normal, in [0, 2π).
    #[inline]
    fn angle_in_plane(&self, a: &Vector3<f64>, b: &Vector3<f64>) -> f64 {
        let hhat = self.hvec.normalize();
        let ang = a.cross(b).dot(&hhat).atan2(a.dot(b));
        ang.rem_euclid(2.0 * std::f64::consts::PI)
    }

    #[inline]
    fn e(&self) -> f64 {
        self.evec.norm()
    }

    #[inline]
    fn q(&self, mu: f64) -> f64 {
        self.hvec.norm().powi(2) / (mu * (1.0 + self.e()))
    }

    #[inline]
    fn inc(&self) -> f64 {
        (self.hvec.x.hypot(self.hvec.y)).atan2(self.hvec.z)
    }

    #[inline]
    fn node(&self) -> f64 {
        if self.nvec.norm() <= 1e-11 * self.hvec.norm() {
            return 0.0;
        }
        self.nvec.y.atan2(self.nvec.x).rem_euclid(2.0 * std::f64::consts::PI)
    }

    #[inline]
    fn arg(&self) -> f64 {
        if self.evec.norm() < 1e-10 {
            return 0.0;
        }
        self.angle_in_plane(&self.node_direction(), &self.evec)
    }

    /// For circular orbits (where perihelion is undefined) this is the argument of latitude,
    /// consistent with `arg == 0`.
    #[inline]
    fn true_anomaly(&self) -> f64 {
        let reference = if self.evec.norm() < 1e-10 { self.node_direction() } else { self.evec };
        self.angle_in_plane(&reference, &self.r)
    }
}

/// Semi-major axis (AU) from the specific orbital energy; negative for unbound orbits.
#[inline]
pub fn semi_major_axis(s: &State, mu: f64) -> f64 {
    -mu / (2.0 * specific_energy(s, mu))
}

/// Specific orbital energy (AU²/day²).
#[inline]
pub fn specific_energy(s: &State, mu: f64) -> f64 {
    let (r, v) = pv(s);
    v.dot(&v) / 2.0 - mu / r.norm()
}

/// Specific angular momentum vector r × v.
#[inline]
pub fn angular_momentum(s: &State) -> Vector3<f64> {
    position(s).cross(&velocity(s))
}

/// Eccentricity vector (points at perihelion, length e).
#[inline]
pub fn eccentricity_vector(s: &State, mu: f64) -> Vector3<f64> {
    Orbit::new(s, mu).evec
}

#[inline]
pub fn eccentricity(s: &State, mu: f64) -> f64 {
    Orbit::new(s, mu).e()
}

/// Perihelion distance (AU).
#[inline]
pub fn perihelion(s: &State, mu: f64) -> f64 {
    Orbit::new(s, mu).q(mu)
}

/// Inclination (radians). Does not depend on μ.
#[inline]
pub fn inclination(s: &State) -> f64 {
    Orbit::new(s, 1.0).inc()
}

/// Longitude of the ascending node (radians, 0 for equatorial orbits). Does not depend on μ.
#[inline]
pub fn node(s: &State) -> f64 {
    Orbit::new(s, 1.0).node()
}

/// Argument of perihelion (radians, 0 for circular orbits; for equatorial orbits this is the
/// longitude of perihelion).
#[inline]
pub fn argument_of_perihelion(s: &State, mu: f64) -> f64 {
    Orbit::new(s, mu).arg()
}

/// True anomaly (radians). For circular orbits this is the argument of latitude.
#[inline]
pub fn true_anomaly(s: &State, mu: f64) -> f64 {
    Orbit::new(s, mu).true_anomaly()
}

/// Conic (eccentric, parabolic or hyperbolic) anomaly (radians).
pub fn conic_anomaly(s: &State, mu: f64) -> Result<f64, Box<dyn std::error::Error>> {
    let o = Orbit::new(s, mu);
    Ok(calc_conic_anomaly_from_true_anomaly(o.e(), o.true_anomaly())?)
}

/// Mean anomaly (radians).
pub fn mean_anomaly(s: &State, mu: f64) -> Result<f64, Box<dyn std::error::Error>> {
    let o = Orbit::new(s, mu);
    let e = o.e();
    let conic = calc_conic_anomaly_from_true_anomaly(e, o.true_anomaly())?;
    Ok(calc_mean_anomaly_from_conic_anomaly(e, conic)?)
}

/// All osculating elements of a state in one pass. Gives the same values as the individual
/// functions in this module, but computes the angular momentum and eccentricity vectors once.
pub fn elements(s: &State, mu: f64) -> Elements {
    let o = Orbit::new(s, mu);
    let e = o.e();
    let true_anomaly = o.true_anomaly();
    let conic_anomaly = calc_conic_anomaly_from_true_anomaly(e, true_anomaly).unwrap_or(f64::NAN);
    let mean_anomaly = if conic_anomaly.is_nan() {
        f64::NAN
    } else {
        calc_mean_anomaly_from_conic_anomaly(e, conic_anomaly).unwrap_or(f64::NAN)
    };
    Elements {
        a: semi_major_axis(s, mu),
        e,
        q: o.q(mu),
        inc: o.inc(),
        node: o.node(),
        arg: o.arg(),
        true_anomaly,
        conic_anomaly,
        mean_anomaly,
    }
}
