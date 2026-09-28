//! Body name <-> NAIF ID mapping.

use std::collections::HashMap;
use std::sync::OnceLock;

use super::builtin_bodies::BUILTIN_BODIES;
use super::error::{Result, SpiceError};
use crate::constants::MASSES;

/// Normalize a body name the way SPICE does: uppercase, trimmed, internal whitespace compressed.
pub fn normalize_name(name: &str) -> String {
    name.split_whitespace().collect::<Vec<_>>().join(" ").to_uppercase()
}

struct BuiltinMaps {
    by_name: HashMap<String, i32>,
    by_code: HashMap<i32, &'static str>,
}

fn builtin() -> &'static BuiltinMaps {
    static MAPS: OnceLock<BuiltinMaps> = OnceLock::new();
    MAPS.get_or_init(|| {
        let mut by_name = HashMap::with_capacity(BUILTIN_BODIES.len());
        let mut by_code = HashMap::new();
        for &(name, code) in BUILTIN_BODIES {
            by_name.insert(normalize_name(name), code);
            by_code.insert(code, name); // last one wins = preferred name
        }
        BuiltinMaps { by_name, by_code }
    })
}

/// Look up a body ID from the built-in table (or parse an integer ID).
pub fn builtin_body_id(name: &str) -> Option<i32> {
    let n = normalize_name(name);
    if let Some(&c) = builtin().by_name.get(&n) {
        return Some(c);
    }
    n.parse::<i32>().ok()
}

/// Preferred built-in name for a body ID.
pub fn builtin_body_name(code: i32) -> Option<&'static str> {
    builtin().by_code.get(&code).copied()
}

/// A body known to SPICE, with the mass information the n-body code needs.
#[derive(Debug, Clone, PartialEq)]
pub struct SpiceBody {
    /// NAIF integer ID.
    pub code: i32,
    /// Name as supplied (uppercased).
    pub name: String,
    /// GM in km^3/s^2 if known (0 otherwise).
    pub gm: f64,
    /// Mass in solar masses from `constants::MASSES` (0 if unknown).
    pub mass: f64,
}

impl SpiceBody {
    pub fn new(code: i32, name: &str, gm: f64) -> Self {
        let mass = MASSES.get(&name.to_lowercase()).copied().unwrap_or(0.0);
        SpiceBody {
            code,
            name: name.to_uppercase(),
            gm,
            mass,
        }
    }

    /// Resolve a body from the built-in NAIF name table (or an integer ID string).
    /// Use [`SpiceKernel::body`](super::SpiceKernel::body) to also honour names defined in
    /// loaded text kernels.
    pub fn from_name(name: &str) -> Result<Self> {
        let code = builtin_body_id(name).ok_or_else(|| SpiceError::UnknownBody(name.to_string()))?;
        Ok(SpiceBody::new(code, &normalize_name(name), 0.0))
    }

    pub fn from_code(code: i32) -> Self {
        let name = builtin_body_name(code).map(|s| s.to_string()).unwrap_or_else(|| code.to_string());
        SpiceBody::new(code, &name, 0.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn builtin_names() {
        assert_eq!(builtin_body_id("earth"), Some(399));
        assert_eq!(builtin_body_id("  Earth   Barycenter "), Some(3));
        assert_eq!(builtin_body_id("jwst"), Some(-170));
        assert_eq!(builtin_body_id("2000001"), Some(2000001));
        assert_eq!(builtin_body_id("ssb"), Some(0));
        assert_eq!(builtin_body_name(399), Some("EARTH"));
        assert_eq!(builtin_body_id("not a body"), None);
    }
}
