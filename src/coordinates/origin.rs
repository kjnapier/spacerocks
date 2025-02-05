use serde::{Serialize, Deserialize};
use crate::errors::OriginError;

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[derive(Default)]
pub enum Origin {
    SUN,
    #[default]
    /// Solar System Barycenter
    SSB,
    /// Define a custom origin with a specified gravitational parameter mu.
    Custom {name: String, mu: f64},
}


/// The coordinates of a SpaceRock must be specified with respect to an origin.
/// The origin can be the Sun, the Solar System Barycenter (SSB) or a custom origin.
/// The custom origin is specified by the gravitational parameter mu, which is the product of 
/// the gravitational constant G and the mass of the origin.
impl Origin {

    /// Create a new custom origin with the specified gravitational parameter mu.
    ///
    /// # Arguments
    /// * `mu` - The gravitational parameter of the custom origin.
    /// * `name` - The name of the custom origin.
    ///
    /// # Example
    /// ```
    /// use spacerocks::coordinates::Origin;
    /// let earth_origin = Origin::new_custom(0.000_000_000_889_954, "EARTH");
    /// ```
    pub fn new_custom(mu: f64, name: &str) -> Origin {
        Origin::Custom { mu, name: name.to_string() }
    }

    
    /// Create an Origin from a string.
    ///
    /// # Arguments
    /// * `s` - The string representation of the Origin (SUN, SSB, or a custom origin).
    ///
    /// # Example
    /// ```
    /// let origin = Origin::from_str("SUN").unwrap();
    /// ```
    pub fn from_str(s: &str) -> Result<Origin, OriginError> {
        match s.to_uppercase().as_str() {
            "SUN" => Ok(Origin::SUN),
            "SSB" => Ok(Origin::SSB),
            _ => Err(OriginError::InvalidOrigin(s.to_string())),
        }
    }

    /// Return the string representation of the Origin.
    ///
    /// # Example
    /// ```
    /// let origin = Origin::from_str("SUN").unwrap();
    /// assert_eq!(origin.as_str(), "SUN");
    /// ```
    pub fn as_str(&self) -> &str {
        match self {
            Origin::SUN => "SUN",
            Origin::SSB => "SSB",
            Origin::Custom { name, .. } => name,
        }
    }

    /// Generate the Origin of the Solar System Barycenter (SSB).
    ///
    /// # Example
    /// ```
    /// let origin = Origin::ssb();
    /// ```
    pub fn ssb() -> Origin {
        Origin::SSB
    }

    /// Generate the Origin of the Sun.
    ///
    /// # Example
    /// ```
    /// let origin = Origin::sun();
    /// ```
    pub fn sun() -> Origin {
        Origin::SUN
    }

    /// Return the gravitational parameter mu of the Origin.
    ///
    /// # Example
    /// ```
    /// let origin = Origin::from_str("SUN").unwrap();
    /// assert_eq!(origin.mu(), 0.000_295_912_208_284_119_5);
    /// ```
    pub fn mu(&self) -> f64 {
        match self {
            Origin::SUN => 0.000_295_912_208_284_119_5,
            Origin::SSB => 2.9630927493968080e-04, //0.00029630927493457475, ,
            Origin::Custom { mu, .. } => *mu,
        }
    }

    /// Return the name of the Origin.
    ///
    /// # Example
    /// ```
    /// let origin = Origin::SSB();
    /// assert_eq!(origin.name(), "SSB");
    /// ```
    pub fn name(&self) -> &str {
        match self {
            Origin::SUN => "SUN",
            Origin::SSB => "SSB",
            Origin::Custom { name, .. } => name,
        }
    }

    /// Return the string representation of the Origin.
    ///
    /// # Example
    /// ```
    /// let origin = Origin::SUN();
    /// assert_eq!(origin.to_string(), "SUN");
    /// ```
    pub fn to_string(&self) -> String {
        self.name().to_string()
    }
}



impl std::fmt::Display for Origin {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", self.name())
    }
}