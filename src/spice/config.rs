use std::path::PathBuf;
use serde::{Deserialize, Serialize};
use crate::errors::KernelError;
use std::fs;
use toml;

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct KernelSpec {
    pub name: String,
    pub kernel_type: String,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct Config {
    #[serde(default)]
    pub default_kernels: Vec<KernelSpec>,
    #[serde(default)]
    pub kernel_paths: Vec<PathBuf>,
    #[serde(default = "default_download_setting")]
    pub auto_download: bool,
    #[serde(default = "default_download_dir")]
    pub download_dir: PathBuf,
}

impl Config {
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

    pub fn from_file(path: &str) -> Result<Self, KernelError> {
        let content = fs::read_to_string(path)
            .map_err(|e| KernelError::IoError(e.to_string()))?;
            
        toml::from_str(&content)
            .map_err(|e| KernelError::ConfigError(e.to_string()))
    }
}

fn default_download_setting() -> bool {
    true
}

fn default_download_dir() -> PathBuf {
    dirs::home_dir()
        .unwrap_or_else(|| PathBuf::from("."))
        .join(".spacerocks")
        .join("spice")
}