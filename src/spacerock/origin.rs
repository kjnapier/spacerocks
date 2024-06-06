use std::borrow::Cow;

use serde::{Serialize, Deserialize};

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum Origin {
    SUN,
    SSB,
    Custom {name: String, mu: f64},
}

impl Origin {

    pub fn new_custom(mu: f64, name: &str) -> Origin {
        Origin::Custom { mu: mu, name: name.to_string() }
    }

    pub fn from_string(s: &str) -> Origin {
        // s.to_upper
        s.to_uppercase();
        match s {
            "SUN" => Origin::SUN,
            "SSB" => Origin::SSB,
            _ => panic!("Unknown origin: {}", s),
        }
    }

    pub fn ssb() -> Origin {
        Origin::SSB
    }

    pub fn sun() -> Origin {
        Origin::SUN
    }

    pub fn mu(&self) -> f64 {
        match self {
            Origin::SUN => 0.00029591220828411951,
            Origin::SSB => 0.00029630927493457475,
            Origin::Custom { mu, .. } => *mu,
        }
    }

    pub fn name(&self) -> &str {
        match self {
            Origin::SUN => "SUN",
            Origin::SSB => "SSB",
            Origin::Custom { name, .. } => name,
        }
    }

    pub fn to_string(&self) -> String {
        self.name().to_string()
    }
}


impl Default for Origin {
    fn default() -> Origin {
        Origin::SSB
    }
}

impl std::fmt::Display for Origin {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", self.name())
    }
}