use std::path::PathBuf;
use std::fs::{self, File};
use std::io::Write;
use std::collections::HashMap;
use std::time::SystemTime;

use crate::errors::KernelError;
use crate::constants::SPICE_URL;
use super::config::{Config, KernelSpec};

#[derive(Debug)]
struct KernelMetadata {
    path: PathBuf,
    kernel_type: String,
    load_time: SystemTime,
}

pub struct SpiceKernel {
    loaded_files: Vec<String>,
    config: Option<Config>,
}

impl SpiceKernel {
    pub fn new() -> Self {
        SpiceKernel {
            loaded_files: vec![],
            config: None,
        }
    }

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
    
    fn load_kernels(&mut self, force_download: Option<bool>) -> Result<(), KernelError> {
        let kernels = self.config.as_ref()
            .map(|c| c.default_kernels.clone())
            .unwrap_or_default();
            
        for kernel in kernels {
            self.process_kernel(&kernel, force_download)?;
        }
        Ok(())
    }

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

    pub fn unload(&mut self) {
        println!("Unloading all kernels");
        spice::kclear();
        self.loaded_files.clear();
    }
    
    pub fn loaded_kernels(&self) -> &[String] {
        &self.loaded_files
    }
}