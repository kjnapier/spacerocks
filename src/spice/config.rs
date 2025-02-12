use std::path::PathBuf;
use serde::{Deserialize, Serialize};
use crate::errors::KernelError;
use std::fs;
use toml;

/// Specification for a SPICE kernel file.
#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct KernelSpec {
    pub name: String,
    pub kernel_type: String,
}

/// Configuration for SPICE kernel management and downloading.
/// 
/// Holds settings for kernel paths, auto-download behavior,
/// and default kernel specifications.
#[derive(Debug, Serialize, Deserialize)]
pub struct Config {
    /// List of default kernels to load on initialization
    #[serde(default)]
    pub default_kernels: Vec<KernelSpec>,
    /// List of directories to search for kernels
    #[serde(default)]
    pub kernel_paths: Vec<PathBuf>,
    /// Whether to automatically download missing kernels
    #[serde(default = "default_download_setting")]
    pub auto_download: bool,
    /// Directory where downloaded kernels are stored
    #[serde(default = "default_download_dir")]
    pub download_dir: PathBuf,
}

impl Config {
    /// Creates a default configuration with specified download behavior.
    /// 
    /// # Arguments
    /// * `download` - Whether to enable automatic downloading of missing kernels
    /// 
    /// # Returns
    /// A new Config instance with default settings and kernel specifications
    pub fn default_with_download(download: bool) -> Self {
        let default_path = dirs::home_dir()
            .unwrap_or_else(|| PathBuf::from("."))
            .join(".spacerocks")
            .join("spice");

        Config {
            default_kernels: vec![
                KernelSpec {
                    name: "latest_leapseconds.tls".to_string(),
                    kernel_type: "lsk".to_string(),
                },
                KernelSpec {
                    name: "de440s.bsp".to_string(),
                    kernel_type: "spk/planets".to_string(),
                },
                KernelSpec {
                    name: "earth_1962_240827_2124_combined.bpc".to_string(),
                    kernel_type: "pck".to_string(),
                },
                KernelSpec {
                    name: "codes_300ast_20100725.bsp".to_string(),
                    kernel_type: "spk/asteroids".to_string(),
                },
                KernelSpec {
                    name: "codes_300ast_20100725.tf".to_string(),
                    kernel_type: "spk/asteroids".to_string(),
                },
            ],
            kernel_paths: vec![default_path.clone()],
            auto_download: download,
            download_dir: default_path,
        }
    }

    /// Loads configuration from a TOML file.
    /// 
    /// # Arguments
    /// * `path` - Path to the configuration file
    /// 
    /// # Returns
    /// Result containing Config on success or KernelError on failure
    /// 
    /// # Errors
    /// Returns KernelError if:
    /// - File cannot be read
    /// - TOML parsing fails

    pub fn from_file(path: &str) -> Result<Self, KernelError> {
        let content = fs::read_to_string(path)
            .map_err(|e| KernelError::IoError(e.to_string()))?;
            
        toml::from_str(&content)
            .map_err(|e| KernelError::ConfigError(e.to_string()))
    }
}

/// Default setting for auto-download behavior
fn default_download_setting() -> bool {
    true
}

/// Default directory for storing downloaded kernels
fn default_download_dir() -> PathBuf {
    dirs::home_dir()
        .unwrap_or_else(|| PathBuf::from("."))
        .join(".spacerocks")
        .join("spice")
}