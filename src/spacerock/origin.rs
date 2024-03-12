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

    pub fn mu(&self) -> f64 {
        match self {
            Origin::Sun => 0.999,
            Origin::Barycenter => 1.000,
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