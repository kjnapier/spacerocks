use crate::spice::spk::Spk;
use crate::spice::bpc::Bpc;
use crate::spice::error::SpiceError;
use crate::spice::spicebody::SpiceBody;
use std::collections::HashMap;

use std::path::PathBuf;
use std::fs::{self, File};
use std::io::Write;
use std::time::SystemTime;

use crate::errors::KernelError;
use crate::constants::SPICE_URL;
use crate::coordinates::ReferencePlane;
use super::config::{Config, KernelSpec};

use nalgebra::Matrix3;

// use spice;

// ----- SpiceKernel with automatic graph building and caching -----

/// Each connection edge now includes the spk_index from which it was derived.
#[derive(Debug, Clone)]
pub struct ConnectionEdge {
    pub child: i32,
    pub parent: i32,
    pub spk_index: usize,
}

/// The SpiceKernel now holds a collection of SPK files and a
/// pre-built connection cache mapping each body to its chain from itself up to 0.
pub struct SpiceKernel {
    pub spk: Vec<Spk>,
    pub bpc: Option<Bpc>,
    /// Maps a body spiceid (child) to its reference (parent) and the SPK index where that relation came from.
    pub parent_map: HashMap<i32, (i32, usize)>,
    /// Cached connection chains from each body up to 0.
    /// Each chain is stored as a vector of ConnectionEdge, running from the body upward.
    pub connection_cache: HashMap<i32, Vec<ConnectionEdge>>,
    // config: Option<Config>,
}

impl SpiceKernel {
    pub fn new() -> Self {
        SpiceKernel {
            spk: Vec::new(),
            bpc: None,
            parent_map: HashMap::new(),
            connection_cache: HashMap::new(),
            //  config: None,
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
    // pub fn defaults(force_download: Option<bool>) -> Result<Self, KernelError> {
    //     let config = Config::default_with_download(true);
    //     let mut kernel = SpiceKernel::new();
    //     kernel.config = Some(config);
        
    //     println!("\nUsing default configuration:");
    //     if let Some(cfg) = &kernel.config {
    //         println!("  Kernel paths: {:?}", cfg.kernel_paths);
    //         println!("  Download directory: {:?}", cfg.download_dir);
    //         println!("  Auto-download: {}", cfg.auto_download);
    //     }
    //     kernel.load_kernels(force_download)?;
    //     Ok(kernel)
    // }

    /// Processes a single kernel specification, attempting to load or download it.
    /// 
    /// # Arguments
    /// * `kernel_spec` - Specification of the kernel to process
    /// * `force_download` - Whether to force download even if kernel exists
    /// 
    /// # Returns
    /// Result indicating success or failure
    // fn process_kernel(&mut self, kernel_spec: &KernelSpec, force_download: Option<bool>) -> Result<(), Box<dyn std::error::Error>> {
    //     // If we have kernel paths, check them first
    //     if let Some(config) = &self.config {
    //         println!("\nProcessing kernel: {}", kernel_spec.name);

    //          // Always download earth orientation file 
    //          let force_download = force_download.unwrap_or(false);
    //          let is_earth_file = kernel_spec.name.starts_with("earth_") && 
    //                            kernel_spec.name.ends_with("_combined.bpc");
             
    //          if force_download || is_earth_file {
    //              println!("➜ Downloading kernel...");
    //              fs::create_dir_all(&config.download_dir)
    //                  .map_err(|e| KernelError::IoError(e.to_string()))?;
                     
    //              let path = self.download_kernel(&kernel_spec.kernel_type, &kernel_spec.name)?;
    //              return self.load(path.to_str().unwrap());
    //          }

    //         // Check each path for existence
    //         for path in &config.kernel_paths {
    //             let kernel_path = path.join(&kernel_spec.name);
    //             if kernel_path.exists() {
    //                 println!("✓ Found existing kernel at: {}", kernel_path.display());
    //                 return self.load(kernel_path.to_str().unwrap());
    //             }
    //         }
            
    //         // Not found in paths - try downloading
    //         println!("➜ Downloading kernel...");
    //         fs::create_dir_all(&config.download_dir)
    //             .map_err(|e| KernelError::IoError(e.to_string()))?;
                
    //         let path = self.download_kernel(&kernel_spec.kernel_type, &kernel_spec.name)?;
    //         return self.load(path.to_str().unwrap());
    //     }
        
    //     Err(Box::new(KernelError::InvalidConfig("No configuration provided".to_string())))
    // }
    
    /// Load all kernels specified in the configuration.
    /// 
    /// # Arguments
    /// * `force_download` - Whether to force download of all kernels
    /// 
    /// # Returns
    /// Result indicating success or failure
    // fn load_kernels(&mut self, force_download: Option<bool>) -> Result<(), KernelError> {
    //     let kernels = self.config.as_ref()
    //         .map(|c| c.default_kernels.clone())
    //         .unwrap_or_default();
            
    //     for kernel in kernels {
    //         self.process_kernel(&kernel, force_download)?;
    //     }
    //     Ok(())
    // }

    /// Downloads a kernel file from the SPICE server.
    /// 
    /// # Arguments
    /// * `kernel_type` - Type of kernel (e.g., "spk/planets")
    /// * `filename` - Name of kernel file to download
    /// 
    /// # Returns
    /// Result containing PathBuf of downloaded file or KernelError
    // fn download_kernel(&self, kernel_type: &str, filename: &str) -> Result<PathBuf, KernelError> {
    //     let config = self.config.as_ref()
    //         .ok_or_else(|| KernelError::InvalidConfig("No configuration provided".to_string()))?;
            
    //     let url = format!("{}/{}/{}", SPICE_URL, kernel_type, filename);
    //     let path = config.download_dir.join(filename);
        
    //     println!("    Downloading from {}", url);
    //     println!("    Saving to {}", path.display());
        
    //     let response = reqwest::blocking::get(&url)
    //         .map_err(|e| KernelError::DownloadError(e.to_string()))?;
            
    //     if !response.status().is_success() {
    //         return Err(KernelError::DownloadError(
    //             format!("Download failed with status: {}", response.status())
    //         ));
    //     }
    //     let content = response.bytes()
    //         .map_err(|e| KernelError::DownloadError(e.to_string()))?;
            
    //     File::create(&path)
    //         .map_err(|e| KernelError::IoError(e.to_string()))?
    //         .write_all(&content)
    //         .map_err(|e| KernelError::IoError(e.to_string()))?;
        
    //     Ok(path)
    // }

    // pub fn load(&mut self, path: &str) -> Result<(), Box<dyn std::error::Error>> {
    //     let metadata = fs::metadata(path)
    //         .map_err(|e| KernelError::IoError(e.to_string()))?;
    //     if metadata.is_file() {
    //         let ext = path.split('.').last().unwrap_or_default();
    //         match ext {
    //             "bsp" => self.load_spk(path)?,
    //             "spk" => self.load_spk(path)?,
    //             "bpc" => self.load_pck(path)?,
    //             "pck" => self.load_pck(path)?,
    //             _ => return Err(Box::new("Unsupported kernel type".to_string())),
    //         }
    //     } else {
    //         return Err(Box::new("Kernel path is not a file".to_string()));
    //     }
    //     Ok(())
    // }

    pub fn load_spk(&mut self, path: &str) -> Result<(), Box<dyn std::error::Error>> {
        let spk = Spk::open(path)?;
        self.add_spk(spk);
        Ok(())
    }

    pub fn load_bpc(&mut self, path: &str) -> Result<(), Box<dyn std::error::Error>> {
        let bpc = Bpc::open(path)?;
        self.bpc = Some(bpc);
        Ok(())
    }

    /// When adding an SPK file, update the parent_map with each target's (child, parent)
    /// relation (recording the spk index) and add it to the list of spks.
    pub fn add_spk(&mut self, spk: Spk) {
        let spk_index = self.spk.len(); // the new spk's index
        for target in &spk.targets {
            // Insert the connection if not already present.
            self.parent_map.entry(target.code).or_insert((target.cen, spk_index));
        }
        self.spk.push(spk);
        self.precompute_connection_chains();
    }

    /// Precompute the connection chains for all bodies (keys in parent_map)
    /// so that each chain from a body to 0 is cached.
    pub fn precompute_connection_chains(&mut self) {
        for &node in self.parent_map.keys() {
            let chain = self.compute_chain_to_root(node);
            self.connection_cache.insert(node, chain);
        }
        // Ensure that node 0 (the absolute frame) has an empty chain.
        if !self.connection_cache.contains_key(&0) {
            self.connection_cache.insert(0, Vec::new());
        }
    }

    /// Compute the connection chain from a given node up to 0 using parent_map.
    /// Returns a vector of ConnectionEdge (child, parent, spk_index).
    fn compute_chain_to_root(&self, mut node: i32) -> Vec<ConnectionEdge> {
        let mut chain = Vec::new();
        while let Some(&(parent, spk_idx)) = self.parent_map.get(&node) {
            if parent == node { break; }
            chain.push(ConnectionEdge { child: node, parent, spk_index: spk_idx });
            node = parent;
        }
        chain
    }

    /// Build a node path (child → … → root) from a starting node and its cached chain.
    /// The returned vector includes the node itself, then each parent's code.
    #[inline(always)]
    fn node_path_from_chain(node: i32, chain: &[ConnectionEdge]) -> Vec<i32> {
        let mut path = Vec::with_capacity(chain.len() + 1);
        path.push(node);
        for edge in chain {
            path.push(edge.parent);
        }
        path
    }

    /// Find the lowest common ancestor (LCA) of two node paths.
    /// The paths are assumed to run from the node up to the root.
    /// This version iterates backwards without allocating new vectors.
    #[inline(always)]
    fn find_lowest_common_ancestor(path1: &[i32], path2: &[i32]) -> Option<i32> {
        let mut i = path1.len();
        let mut j = path2.len();
        while i > 0 && j > 0 && path1[i - 1] == path2[j - 1] {
            i -= 1;
            j -= 1;
        }
        if i < path1.len() { Some(path1[i]) } else { None }
    }

    /// Compute the connection chain from `start` to `target` as a vector of ConnectionEdge.
    /// The chain is built by combining the cached chain from `start` to 0 and from `target` to 0.
    /// Edges from the target's chain are inverted.
    ///
    /// For example, if the cached chains are:
    ///   start=200001: [(200001, 10, 1), (10, 0, 0)]
    ///   target=10: [(10, 0, 0)]
    /// then connection_chain(200001, 0) returns [(200001, 10, 1), (10, 0, 0)]
    /// and for (301, 10), if
    ///   301: [(301, 3, 1), (3, 0, 0)]
    ///   10: [(10, 0, 0)]
    /// then it returns [(301, 3, 1), (3, 0, 0), (0, 10, 0)].
    pub fn connection_chain(&self, start: i32, target: i32) -> Option<Vec<ConnectionEdge>> {
        let chain_start = self.connection_cache.get(&start)?;
        let binding = Vec::new();
        let chain_target = self.connection_cache.get(&target).unwrap_or(&binding);
        let path_start = Self::node_path_from_chain(start, chain_start);
        let path_target = Self::node_path_from_chain(target, chain_target);
        let lca = Self::find_lowest_common_ancestor(&path_start, &path_target)?;
        
        // Build the upward chain from start to LCA using parent_map (which is fast).
        let mut chain_up = Vec::new();
        let mut node = start;
        while node != lca {
            let &(parent, spk_idx) = self.parent_map.get(&node)?;
            chain_up.push(ConnectionEdge { child: node, parent, spk_index: spk_idx });
            node = parent;
        }
        // Similarly, build the upward chain from target to LCA.
        let mut chain_down = Vec::new();
        let mut node = target;
        while node != lca {
            let &(parent, spk_idx) = self.parent_map.get(&node)?;
            chain_down.push(ConnectionEdge { child: node, parent, spk_index: spk_idx });
            node = parent;
        }
        // Invert the downward chain to run from LCA down to target.
        chain_down.reverse();
        let inverted_chain_down: Vec<ConnectionEdge> = chain_down
            .into_iter()
            .map(|edge| ConnectionEdge {
                child: edge.parent,
                parent: edge.child,
                spk_index: edge.spk_index,
            })
            .collect();
        chain_up.extend(inverted_chain_down);
        Some(chain_up)
    }
}

// --- Now we add a method to compute the overall state via a connection chain ---
impl SpiceKernel {
    /// Compute the relative state (position and velocity) of `from` with respect to `to`
    /// at the given epoch by using the pre-computed connection chain.
    ///
    /// For each edge in the chain, we look up the SPK (via spk_index), use its target_map to
    /// find the target index for the child body, and call that SPK’s state_at() to get the state vector.
    /// We then sum the position and velocity components.
    pub fn compute_state(&self, body: &SpiceBody, origin: &SpiceBody, epoch: f64) -> Result<(f64, f64, f64, f64, f64, f64), SpiceError> {
        let from = origin.code;
        let to = body.code;
        let chain = self.connection_chain(from, to)
            .ok_or_else(|| SpiceError::ParseError(format!("No connection chain found from {} to {}", from, to)))?;

        let mut pos = [0.0; 3];
        let mut vel = [0.0; 3];
        
        for edge in chain {
            let spk = &self.spk[edge.spk_index];
            let target_index = spk.target_map.get(&edge.parent)
                .ok_or_else(|| SpiceError::ParseError(format!("Body {} not found in SPK index {}", edge.child, edge.spk_index)))?;

            let state = spk.state_at(epoch, *target_index)?;
            pos[0] += state.0;
            pos[1] += state.1;
            pos[2] += state.2;
            vel[0] += state.3;
            vel[1] += state.4;
            vel[2] += state.5;
        }
        Ok((pos[0], pos[1], pos[2], vel[0], vel[1], vel[2]))
    }

    pub fn pxform(&self, epoch: f64) -> Result<[[f64; 3]; 3], Box<dyn std::error::Error>> {
        match self.bpc {
            Some(ref bpc) => {
                let m = bpc.rotation_matrix_at_epoch(epoch)?;

                let current_reference_plane = ReferencePlane::from_str("ECLIPJ2000")?;
                let reference_plane = ReferencePlane::from_str("J2000")?;
                let inv = reference_plane.get_rotation_matrix().try_inverse().ok_or("Could not invert rotation matrix")?;
                let rot = current_reference_plane.get_rotation_matrix() * inv;

                // println!("Rotation matrix: {:?}", rot);
                // turn the m matrix into a 3x3 matrix

                let m: Matrix3<f64> = m.into();
                Ok((m * rot).into())
            },
            None => {
                Err(Box::new(SpiceError::ParseError("No BPC loaded".to_string())))
            }
        }
    }
}


/// Metadata for a loaded SPICE kernel
#[derive(Debug)]
struct KernelMetadata {
    path: PathBuf,
    kernel_type: String,
    load_time: SystemTime,
}