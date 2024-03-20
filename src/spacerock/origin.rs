use std::borrow::Cow;

#[derive(Clone)]
pub enum Origin {
    Sun,
    Barycenter,
    Custom {name: Cow<'static, str>, mu: f64},
}

impl Origin {

    pub fn new_custom(mu: f64, name: &'static str) -> Origin {
        Origin::Custom { mu: mu, name: Cow::Borrowed(name) }
    }

    pub fn from_string(s: &str) -> Origin {
        match s {
            "Sun" => Origin::Sun,
            "Barycenter" => Origin::Barycenter,
            _ => panic!("Unknown origin: {}", s),
        }
    }

    pub fn mu(&self) -> f64 {
        match self {
            Origin::Sun => 0.00029591220828411951,
            Origin::Barycenter => 0.00029630927493457475,
            Origin::Custom { mu, .. } => *mu,
        }
    }

    pub fn name(&self) -> &str {
        match self {
            Origin::Sun => "Sun",
            Origin::Barycenter => "Barycenter",
            Origin::Custom { name, .. } => name,
        }
    }
}


impl Default for Origin {
    fn default() -> Origin {
        Origin::Barycenter
    }
}