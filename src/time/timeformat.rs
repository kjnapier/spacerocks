#[derive(Debug, Clone, PartialEq)]
pub enum TimeFormat {
    JD,
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