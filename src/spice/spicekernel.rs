use std::path::PathBuf;
use std::fs::{self, File};
use std::io::Write;
use std::collections::HashMap;
use std::time::SystemTime;

use crate::errors::KernelError;
use crate::constants::SPICE_URL;
use super::config::{Config, KernelSpec};

/// Metadata for a loaded SPICE kernel
#[derive(Debug)]
struct KernelMetadata {
    path: PathBuf,
    kernel_type: String,
    load_time: SystemTime,
}

/// Manages SPICE kernel loading, downloading, and configuration.
/// 
/// Handles:
/// - Loading and unloading SPICE kernel files
/// - Downloading missing kernels
/// - Tracking loaded kernel state
pub struct SpiceKernel {
    /// Currently loaded kernel file paths
    loaded_files: Vec<String>,
    /// Optional configuration for kernel management
    config: Option<Config>,
}

impl SpiceKernel {
    /// Creates a new empty SpiceKernel instance.
    pub fn new() -> Self {
        SpiceKernel {
            loaded_files: vec![],
            config: None,
        }
    }

    /// Creates a SpiceKernel with default configuration.
    /// 
    /// # Arguments
    /// * `force_download` - Optional flag to force download of kernels
    /// 
    /// # Returns
    /// * [`SpiceKernel`]
    ///
    /// # Errors
    /// Returns [`KernelError`] if:
    /// - Configuration file cannot be read
    /// - Kernel loading fails
    /// - Downloads fail when required

    pub fn defaults(force_download: Option<bool>) -> Result<Self, KernelError> {
        let config = Config::default_with_download(true);
        let mut kernel = SpiceKernel {
            loaded_files: vec![],
            config: Some(config),
        };
        
        println!("\nUsing default configuration:");
        if let Some(cfg) = &kernel.config {
            println!("  Kernel paths: {:?}", cfg.kernel_paths);
            println!("  Download directory: {:?}", cfg.download_dir);
            println!("  Auto-download: {}", cfg.auto_download);
        }
        kernel.load_kernels(force_download)?;
        Ok(kernel)
    }

    /// Creates a SpiceKernel from a configuration file.
    /// 
    /// # Arguments
    /// * `path` - Path to configuration file
    /// * `force_download` - Optional flag to force download of kernels
    /// 
    /// # Returns
    /// Result containing SpiceKernel on success or KernelError on failure
    /// 
    /// # Errors
    /// Returns KernelError if:
    /// - Configuration file cannot be read
    /// - Kernel loading fails
    /// - Downloads fail when required
    pub fn from_config(path: &str, force_download: Option<bool>) -> Result<Self, KernelError> {
        println!("Loading configuration from {}", path);
        
        let config = Config::from_file(path)?;
        let mut kernel = SpiceKernel {
            loaded_files: vec![],
            config: Some(config),
        };
        
        println!("\nConfiguration loaded:");
        if let Some(cfg) = &kernel.config {
            println!("  Kernel paths: {:?}", cfg.kernel_paths);
            println!("  Download directory: {:?}", cfg.download_dir);
            println!("  Auto-download: {}", cfg.auto_download);
        }
        
        kernel.load_kernels(force_download)?;
        Ok(kernel)
    }

    /// Processes a single kernel specification, attempting to load or download it.
    /// 
    /// # Arguments
    /// * `kernel_spec` - Specification of the kernel to process
    /// * `force_download` - Whether to force download even if kernel exists
    /// 
    /// # Returns
    /// Result indicating success or failure
    fn process_kernel(&mut self, kernel_spec: &KernelSpec, force_download: Option<bool>) -> Result<(), KernelError> {
        // If we have kernel paths, check them first
        if let Some(config) = &self.config {
            println!("\nProcessing kernel: {}", kernel_spec.name);

             // Always download earth orientation file 
             let force_download = force_download.unwrap_or(false);
             let is_earth_file = kernel_spec.name.starts_with("earth_") && 
                               kernel_spec.name.ends_with("_combined.bpc");
             
             if force_download || is_earth_file {
                 println!("➜ Downloading kernel...");
                 fs::create_dir_all(&config.download_dir)
                     .map_err(|e| KernelError::IoError(e.to_string()))?;
                     
                 let path = self.download_kernel(&kernel_spec.kernel_type, &kernel_spec.name)?;
                 return self.load(path.to_str().unwrap());
             }

            // Check each path for existence
            for path in &config.kernel_paths {
                let kernel_path = path.join(&kernel_spec.name);
                if kernel_path.exists() {
                    println!("✓ Found existing kernel at: {}", kernel_path.display());
                    return self.load(kernel_path.to_str().unwrap());
                }
            }
            
            // Not found in paths - try downloading
            println!("➜ Downloading kernel...");
            fs::create_dir_all(&config.download_dir)
                .map_err(|e| KernelError::IoError(e.to_string()))?;
                
            let path = self.download_kernel(&kernel_spec.kernel_type, &kernel_spec.name)?;
            return self.load(path.to_str().unwrap());
        }
        
        Err(KernelError::InvalidConfig("No configuration provided".to_string()))
    }
    
    /// Load all kernels specified in the configuration.
    /// 
    /// # Arguments
    /// * `force_download` - Whether to force download of all kernels
    /// 
    /// # Returns
    /// Result indicating success or failure
    fn load_kernels(&mut self, force_download: Option<bool>) -> Result<(), KernelError> {
        let kernels = self.config.as_ref()
            .map(|c| c.default_kernels.clone())
            .unwrap_or_default();
            
        for kernel in kernels {
            self.process_kernel(&kernel, force_download)?;
        }
        Ok(())
    }

    /// Downloads a kernel file from the SPICE server.
    /// 
    /// # Arguments
    /// * `kernel_type` - Type of kernel (e.g., "spk/planets")
    /// * `filename` - Name of kernel file to download
    /// 
    /// # Returns
    /// Result containing PathBuf of downloaded file or KernelError
    fn download_kernel(&self, kernel_type: &str, filename: &str) -> Result<PathBuf, KernelError> {
        let config = self.config.as_ref()
            .ok_or_else(|| KernelError::InvalidConfig("No configuration provided".to_string()))?;
            
        let url = format!("{}/{}/{}", SPICE_URL, kernel_type, filename);
        let path = config.download_dir.join(filename);
        
        println!("    Downloading from {}", url);
        println!("    Saving to {}", path.display());
        
        let response = reqwest::blocking::get(&url)
            .map_err(|e| KernelError::DownloadError(e.to_string()))?;
            
        if !response.status().is_success() {
            return Err(KernelError::DownloadError(
                format!("Download failed with status: {}", response.status())
            ));
        }
        let content = response.bytes()
            .map_err(|e| KernelError::DownloadError(e.to_string()))?;
            
        File::create(&path)
            .map_err(|e| KernelError::IoError(e.to_string()))?
            .write_all(&content)
            .map_err(|e| KernelError::IoError(e.to_string()))?;
        
        Ok(path)
    }

    /// Loads a SPICE kernel file.
    /// 
    /// # Arguments
    /// * `path` - Path to kernel file
    /// 
    /// # Returns
    /// Result containing unit on success or KernelError on failure
    pub fn load(&mut self, path: &str) -> Result<(), KernelError> {
        if self.loaded_files.contains(&path.to_string()) {
            println!("Kernel already loaded: {}", path);
            return Ok(());
        }
        
        println!("Loading kernel: {}", path);
        spice::furnsh(path);
        self.loaded_files.push(path.to_string());
        Ok(())
    }

    /// Unloads all currently loaded kernels.
    pub fn unload(&mut self) {
        println!("Unloading all kernels");
        spice::kclear();
        self.loaded_files.clear();
    }
    
    /// Returns slice of currently loaded kernel file paths.
    pub fn loaded_kernels(&self) -> &[String] {
        &self.loaded_files
    }
}