#[derive(Debug, Clone, PartialEq)]
pub enum TimeFormat {
    JD,
}

impl TimeFormat {
    pub fn from_str(s: &str) -> Option<Self> {
        match s {
            "JD" => Some(TimeFormat::JD),
            _ => None,
        }
    }

    pub fn to_str(&self) -> &str {
        match self {
            TimeFormat::JD => "JD",
        }
    }
}

impl std::fmt::Display for TimeFormat {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            TimeFormat::JD => write!(f, "JD"),
        }
    }
}

impl Default for TimeFormat {
    fn default() -> Self {
        TimeFormat::JD
    }
}