//! Origin of chosen reference frame.

use std::collections::HashMap;
use std::sync::{OnceLock, RwLock};

use crate::errors::OriginError;
use crate::spice::bodies::{builtin_body_id, builtin_body_name};

/// The center of a coordinate system: a body id and its gravitational parameter μ (AU³/day²).
///
/// The id is the body's NAIF code when it has one (10 for the Sun, 0 for the Solar System
/// Barycenter, 399 for the Earth, ...). Origins that are not SPICE bodies, such as a particle of
/// a [`crate::nbody::Simulation`] or its barycenter, get an id of their own from
/// [`Origin::new_custom`], which remembers the name. `Origin` is `Copy`, so populations and
/// rocks carry it by value.
///
/// μ is used for the two-body (Keplerian) calculations about the origin.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Origin {
    /// NAIF id of the body, or an id from [`Origin::new_custom`] for other origins.
    pub id: i32,
    /// Gravitational parameter (AU³/day²).
    pub mu: f64,
}

impl Default for Origin {
    fn default() -> Self {
        Origin::SSB
    }
}

/// Names of origins that are not in the NAIF table: id `CUSTOM_BASE + k` is `names[k]`.
struct Registry {
    names: Vec<&'static str>,
    ids: HashMap<&'static str, i32>,
    /// Names made up for unnamed NAIF ids.
    numeric: HashMap<i32, &'static str>,
}

const CUSTOM_BASE: i32 = i32::MIN;

fn registry() -> &'static RwLock<Registry> {
    static REG: OnceLock<RwLock<Registry>> = OnceLock::new();
    REG.get_or_init(|| RwLock::new(Registry { names: Vec::new(), ids: HashMap::new(), numeric: HashMap::new() }))
}

impl Origin {
    /// The Sun.
    pub const SUN: Origin = Origin { id: 10, mu: 0.000_295_912_208_284_119_5 };
    /// The Solar System Barycenter.
    pub const SSB: Origin = Origin { id: 0, mu: 2.9630927493968080e-04 };

    /// An origin with NAIF id `id` and gravitational parameter `mu`.
    pub fn new(id: i32, mu: f64) -> Origin {
        Origin { id, mu }
    }

    /// An origin named `name` with gravitational parameter `mu`. A name in the NAIF table (or
    /// an integer id) gets that body's id; any other name gets an id of its own, the same one
    /// every time it is used.
    ///
    /// # Example
    /// ```
    /// use spacerocks::coordinates::Origin;
    /// let earth = Origin::new_custom(0.000_000_000_889_954, "EARTH");
    /// assert_eq!(earth.id, 399);
    /// let other = Origin::new_custom(1e-10, "my barycenter");
    /// assert_eq!(other.name(), "my barycenter");
    /// ```
    pub fn new_custom(mu: f64, name: &str) -> Origin {
        if let Some(id) = builtin_body_id(name) {
            return Origin { id, mu };
        }
        if let Some(&id) = registry().read().unwrap().ids.get(name) {
            return Origin { id, mu };
        }
        let mut reg = registry().write().unwrap();
        if let Some(&id) = reg.ids.get(name) {
            return Origin { id, mu };
        }
        let id = CUSTOM_BASE + reg.names.len() as i32;
        let name: &'static str = Box::leak(name.to_string().into_boxed_str());
        reg.names.push(name);
        reg.ids.insert(name, id);
        Origin { id, mu }
    }

    /// SUN or SSB, by name.
    ///
    /// # Example
    /// ```
    /// # use spacerocks::coordinates::Origin;
    /// let origin = Origin::from_str("SUN").unwrap();
    /// ```
    pub fn from_str(s: &str) -> Result<Origin, OriginError> {
        match s.to_uppercase().as_str() {
            "SUN" => Ok(Origin::SUN),
            "SSB" => Ok(Origin::SSB),
            _ => Err(OriginError::InvalidOrigin(s.to_string())),
        }
    }

    /// Whether this is neither the Sun nor the Solar System Barycenter.
    pub fn is_custom(&self) -> bool {
        *self != Origin::SUN && *self != Origin::SSB
    }

    /// The Solar System Barycenter.
    pub fn ssb() -> Origin {
        Origin::SSB
    }

    /// The Sun.
    pub fn sun() -> Origin {
        Origin::SUN
    }

    /// The gravitational parameter μ (AU³/day²).
    ///
    /// # Example
    /// ```
    /// # use spacerocks::coordinates::Origin;
    /// let origin = Origin::from_str("SUN").unwrap();
    /// assert_eq!(origin.mu(), 0.000_295_912_208_284_119_5);
    /// ```
    pub fn mu(&self) -> f64 {
        self.mu
    }

    /// The origin's name: SUN and SSB, the NAIF name of other SPICE bodies, or the name given
    /// to [`Origin::new_custom`].
    ///
    /// # Example
    /// ```
    /// # use spacerocks::coordinates::Origin;
    /// assert_eq!(Origin::ssb().name(), "SSB");
    /// assert_eq!(Origin::sun().name(), "SUN");
    /// ```
    pub fn name(&self) -> &'static str {
        match self.id {
            0 => return "SSB",
            10 => return "SUN",
            _ => {}
        }
        if self.id >= CUSTOM_BASE && self.id < CUSTOM_BASE + (1 << 30) {
            if let Some(&n) = registry().read().unwrap().names.get((self.id - CUSTOM_BASE) as usize) {
                return n;
            }
        }
        if let Some(n) = builtin_body_name(self.id) {
            return n;
        }
        if let Some(&n) = registry().read().unwrap().numeric.get(&self.id) {
            return n;
        }
        let mut reg = registry().write().unwrap();
        let n: &'static str = Box::leak(self.id.to_string().into_boxed_str());
        *reg.numeric.entry(self.id).or_insert(n)
    }

    /// Same as [`Origin::name`].
    pub fn as_str(&self) -> &'static str {
        self.name()
    }
}

impl std::fmt::Display for Origin {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", self.name())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ids_and_names() {
        assert_eq!(Origin::SUN.name(), "SUN");
        assert_eq!(Origin::SSB.to_string(), "SSB");
        let e = Origin::new_custom(1e-9, "earth");
        assert_eq!((e.id, e.name()), (399, "EARTH"));
        let a = Origin::new_custom(1e-9, "simulation_barycenter");
        let b = Origin::new_custom(2e-9, "simulation_barycenter");
        assert_eq!(a.id, b.id);
        assert_ne!(a, b);
        assert_eq!(b.name(), "simulation_barycenter");
        assert!(a.is_custom() && !Origin::SUN.is_custom());
        assert_eq!(Origin::new(2000001, 0.0).name(), "CERES");
        assert_eq!(Origin::new(123456789, 0.0).name(), "123456789");
    }
}
