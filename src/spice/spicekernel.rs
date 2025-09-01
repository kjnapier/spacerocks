use crate::spice::spk::Spk;
use crate::spice::bpc::Bpc;
use crate::spice::error::SpiceError;
use crate::spice::spicebody::SpiceBody;
use std::collections::{HashMap, HashSet};

use std::path::PathBuf;
use std::fs::{self, File};
use std::io::Write;
use std::time::SystemTime;

use crate::errors::KernelError;
use crate::constants::SPICE_URL;
use crate::coordinates::ReferencePlane;
use super::config::{Config, KernelSpec};

use nalgebra::Matrix3;

use std::sync::{Arc, Mutex};

use rayon::prelude::*;

// ----- SpiceKernel with automatic graph building and caching -----

/// Each connection edge now includes the spk_index from which it was derived.
#[derive(Debug, Clone)]
pub struct ConnectionEdge {
    pub child: i32,
    pub parent: i32,
    pub spk_index: usize,
}

// The SpiceKernel now holds a collection of SPK files and a
// pre-built connection cache mapping each body to its chain from itself up to 0.
// #[derive(Debug)]
// pub struct SpiceKernel {
//     pub spk: Vec<Spk>,
//     pub bpc: Option<Bpc>,
//     /// Maps a body spiceid (child) to its reference (parent) and the SPK index where that relation came from.
//     pub parent_map: HashMap<i32, (i32, usize)>,
//     /// Cached connection chains from each body up to 0.
//     /// Each chain is stored as a vector of ConnectionEdge, running from the body upward.
//     // pub connection_cache: HashMap<i32, Vec<ConnectionEdge>>,
//     pub connection_cache: HashMap<i32, Vec<ConnectionEdge>>,
//     // config: Option<Config>,
// }

#[derive(Debug)]
pub struct SpiceKernel {
    pub spk: Vec<Spk>,
    pub bpc: Option<Bpc>,
    /// Maps a body spiceid (child) to its reference (parent) and the SPK index where that relation came from.
    pub parent_map: HashMap<i32, (i32, usize)>,
    /// Cached connection chains from each body up to 0.
    pub connection_cache: HashMap<i32, Vec<ConnectionEdge>>,
    /// New: Cache for computed connection chains between arbitrary nodes.
    pub connection_chain_cache: Mutex<HashMap<(i32, i32), Vec<ConnectionEdge>>>,
    /// Optionally, cache the full node path for each body to avoid re-allocation.
    pub node_path_cache: Mutex<HashMap<i32, Vec<i32>>>,
}

impl SpiceKernel {
    // pub fn new() -> Self {
    //     SpiceKernel {
    //         spk: Vec::new(),
    //         bpc: None,
    //         parent_map: HashMap::new(),
    //         connection_cache: HashMap::new(),
    //     }
    // }

    pub fn new() -> Self {
        SpiceKernel {
            spk: Vec::new(),
            bpc: None,
            parent_map: HashMap::new(),
            connection_cache: HashMap::new(),
            connection_chain_cache: Mutex::new(HashMap::new()),
            node_path_cache: Mutex::new(HashMap::new()),
        }
    }

    /// Precompute connection chains and optionally node paths.
    pub fn precompute_connection_chains(&mut self) {
        for &node in self.parent_map.keys() {
            let chain = self.compute_chain_to_root(node);
            self.connection_cache.insert(node, chain);
            // Also compute and cache the node path from this node to the root.
            let path = Self::node_path_from_chain(node, self.connection_cache.get(&node).unwrap());
            self.node_path_cache.lock().unwrap().insert(node, path);
        }
        // Ensure that node 0 (the absolute frame) has an empty chain and path.
        if !self.connection_cache.contains_key(&0) {
            self.connection_cache.insert(0, Vec::new());
            self.node_path_cache.lock().unwrap().insert(0, vec![0]);
        }
        // Clear the connection chain cache because any queries will now benefit from the new precomputed data.
        self.connection_chain_cache.lock().unwrap().clear();
    }

    /// Compute the connection chain from a given node up to 0 using parent_map.
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
    /// This version checks the node_path_cache first.
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

    pub fn connection_chain(&self, start: i32, target: i32) -> Option<Vec<ConnectionEdge>> {
        let key = (start, target);
        // Check if we already computed this chain.
        if let Ok(cache) = self.connection_chain_cache.lock() {
            if let Some(cached_chain) = cache.get(&key) {
                return Some(cached_chain.clone());
            }
        }
        
        // Get the node paths from the cache (or compute them if missing).
        let path_start = if let Some(p) = self.node_path_cache.lock().unwrap().get(&start) {
            p.clone()
        } else {
            // Fallback to computing from connection_cache if not pre-cached.
            let chain_start = self.connection_cache.get(&start)?;
            Self::node_path_from_chain(start, chain_start)
        };
        let path_target = if let Some(p) = self.node_path_cache.lock().unwrap().get(&target) {
            p.clone()
        } else {
            // Use a static empty slice instead of a temporary Vec.
            let chain_target = self.connection_cache.get(&target)
                .map(|v| v.as_slice())
                .unwrap_or(&[]);
            Self::node_path_from_chain(target, chain_target)
        };
    
        let lca = Self::find_lowest_common_ancestor(&path_start, &path_target)?;
        
        // Build upward chain from start to LCA.
        let mut chain_up = Vec::new();
        let mut node = start;
        while node != lca {
            let &(parent, spk_idx) = self.parent_map.get(&node)?;
            chain_up.push(ConnectionEdge { child: node, parent, spk_index: spk_idx });
            node = parent;
        }
        // Build upward chain from target to LCA.
        let mut chain_down = Vec::new();
        let mut node = target;
        while node != lca {
            let &(parent, spk_idx) = self.parent_map.get(&node)?;
            chain_down.push(ConnectionEdge { child: node, parent, spk_index: spk_idx });
            node = parent;
        }
        chain_down.reverse();
        // Invert the downward chain to go from LCA down to target.
        let inverted_chain_down: Vec<ConnectionEdge> = chain_down
            .into_iter()
            .map(|edge| ConnectionEdge {
                child: edge.parent,
                parent: edge.child,
                spk_index: edge.spk_index,
            })
            .collect();
        chain_up.extend(inverted_chain_down);
        
        // Cache the computed chain.
        if let Ok(mut cache) = self.connection_chain_cache.lock() {
            cache.insert(key, chain_up.clone());
        }
        Some(chain_up)
    }

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
            // println!("Adding target {} with parent {} to parent_map", target.code, target.cen);
        }
        self.spk.push(spk);
        self.precompute_connection_chains();
    }

    // /// Precompute the connection chains for all bodies (keys in parent_map)
    // /// so that each chain from a body to 0 is cached.
    // pub fn precompute_connection_chains(&mut self) {
    //     for &node in self.parent_map.keys() {
    //         let chain = self.compute_chain_to_root(node);
    //         self.connection_cache.insert(node, chain);
    //     }
    //     // Ensure that node 0 (the absolute frame) has an empty chain.
    //     if !self.connection_cache.contains_key(&0) {
    //         self.connection_cache.insert(0, Vec::new());
    //     }
    // }

    // /// Compute the connection chain from a given node up to 0 using parent_map.
    // /// Returns a vector of ConnectionEdge (child, parent, spk_index).
    // fn compute_chain_to_root(&self, mut node: i32) -> Vec<ConnectionEdge> {
    //     let mut chain = Vec::new();
    //     while let Some(&(parent, spk_idx)) = self.parent_map.get(&node) {
    //         if parent == node { break; }
    //         chain.push(ConnectionEdge { child: node, parent, spk_index: spk_idx });
    //         node = parent;
    //     }
    //     chain
    // }

    // /// Build a node path (child → … → root) from a starting node and its cached chain.
    // /// The returned vector includes the node itself, then each parent's code.
    // #[inline(always)]
    // fn node_path_from_chain(node: i32, chain: &[ConnectionEdge]) -> Vec<i32> {
    //     let mut path = Vec::with_capacity(chain.len() + 1);
    //     path.push(node);
    //     for edge in chain {
    //         path.push(edge.parent);
    //     }
    //     path
    // }

    // /// Find the lowest common ancestor (LCA) of two node paths.
    // /// The paths are assumed to run from the node up to the root.
    // /// This version iterates backwards without allocating new vectors.
    // #[inline(always)]
    // fn find_lowest_common_ancestor(path1: &[i32], path2: &[i32]) -> Option<i32> {
    //     let mut i = path1.len();
    //     let mut j = path2.len();
    //     while i > 0 && j > 0 && path1[i - 1] == path2[j - 1] {
    //         i -= 1;
    //         j -= 1;
    //     }
    //     if i < path1.len() { Some(path1[i]) } else { None }
    // }

    // /// Compute the connection chain from `start` to `target` as a vector of ConnectionEdge.
    // /// The chain is built by combining the cached chain from `start` to 0 and from `target` to 0.
    // /// Edges from the target's chain are inverted.
    // ///
    // /// For example, if the cached chains are:
    // ///   start=200001: [(200001, 10, 1), (10, 0, 0)]
    // ///   target=10: [(10, 0, 0)]
    // /// then connection_chain(200001, 0) returns [(200001, 10, 1), (10, 0, 0)]
    // /// and for (301, 10), if
    // ///   301: [(301, 3, 1), (3, 0, 0)]
    // ///   10: [(10, 0, 0)]
    // /// then it returns [(301, 3, 1), (3, 0, 0), (0, 10, 0)].
    // pub fn connection_chain(&self, start: i32, target: i32) -> Option<Vec<ConnectionEdge>> {

    //     let chain_start = self.connection_cache.get(&start)?;
    //     let binding = Vec::new();
    //     let chain_target = self.connection_cache.get(&target).unwrap_or(&binding);
    //     let path_start = Self::node_path_from_chain(start, chain_start);
    //     let path_target = Self::node_path_from_chain(target, chain_target);
    //     let lca = Self::find_lowest_common_ancestor(&path_start, &path_target)?;
        
    //     // Build the upward chain from start to LCA using parent_map (which is fast).
    //     let mut chain_up = Vec::new();
    //     let mut node = start;
    //     while node != lca {
    //         let &(parent, spk_idx) = self.parent_map.get(&node)?;
    //         chain_up.push(ConnectionEdge { child: node, parent, spk_index: spk_idx });
    //         node = parent;
    //     }
    //     // Similarly, build the upward chain from target to LCA.
    //     let mut chain_down = Vec::new();
    //     let mut node = target;
    //     while node != lca {
    //         let &(parent, spk_idx) = self.parent_map.get(&node)?;
    //         chain_down.push(ConnectionEdge { child: node, parent, spk_index: spk_idx });
    //         node = parent;
    //     }
    //     // Invert the downward chain to run from LCA down to target.
    //     chain_down.reverse();
    //     let inverted_chain_down: Vec<ConnectionEdge> = chain_down
    //         .into_iter()
    //         .map(|edge| ConnectionEdge {
    //             child: edge.parent,
    //             parent: edge.child,
    //             spk_index: edge.spk_index,
    //         })
    //         .collect();
    //     chain_up.extend(inverted_chain_down);
    //     Some(chain_up)
    // }
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
        // let chain = self.connection_chain(from, to)
        //     .ok_or_else(|| SpiceError::ParseError(format!("No connection chain found from {} to {}", from, to)))?;


        // try to get the chain from the cache
        let chain = self.connection_chain(from, to)
            .ok_or_else(|| SpiceError::ParseError(format!("No connection chain found from {} to {}", from, to)))?;


        let mut pos = [0.0; 3];
        let mut vel = [0.0; 3];
        
        for edge in chain {
            let spk = &self.spk[edge.spk_index];
            // println!("Using SPK index {}", edge.spk_index);
            // println!("Using target map {:?}", spk.target_map);
            // println!("Using targets {:?}", spk.targets);
            // println!("Using edge {:?}", edge);
            // println!("parent {:?}", edge.parent);
            // println!("child {:?}", edge.child);

            // if the child is not in the target map, make the parent the child and the child the parent            

            // let target_index = spk.target_map.get(target)
            //     .ok_or_else(|| SpiceError::ParseError(format!("Body {} not found in SPK index {}", edge.child, edge.spk_index)))?;

            // if edge.parent is in the target map:
            if spk.target_map.contains_key(&edge.parent) {
                let target_index = spk.target_map.get(&edge.parent)
                .ok_or_else(|| SpiceError::ParseError(format!("Body {} not found in SPK index {}", edge.child, edge.spk_index)))?;

                let state = spk.state_at(epoch, *target_index)?;
                pos[0] += state.0;
                pos[1] += state.1;
                pos[2] += state.2;
                vel[0] += state.3;
                vel[1] += state.4;
                vel[2] += state.5;
            } else if spk.target_map.contains_key(&edge.child) {
                let target_index = spk.target_map.get(&edge.child)
                .ok_or_else(|| SpiceError::ParseError(format!("Body {} not found in SPK index {}", edge.child, edge.spk_index)))?;

                let state = spk.state_at(epoch, *target_index)?;
                pos[0] -= state.0;
                pos[1] -= state.1;
                pos[2] -= state.2;
                vel[0] -= state.3;
                vel[1] -= state.4;
                vel[2] -= state.5;
            }
            
        }
        Ok((pos[0], pos[1], pos[2], vel[0], vel[1], vel[2]))
    }

    pub fn get_barycentric_states(&self, bodies: &Vec<SpiceBody>, epoch: f64) -> Result<Vec<(f64, f64, f64, f64, f64, f64)>, SpiceError> {
        // get the index and spk index for each body
        let mut indices = Vec::with_capacity(bodies.len());
        let mut spk_indices = Vec::with_capacity(bodies.len());

        for body in bodies {
            let index = body.code;
            let spk_index = self.parent_map.get(&index)
                .ok_or_else(|| SpiceError::ParseError(format!("Body {} not found in parent map", index)))?.1;
            indices.push(index);
            spk_indices.push(spk_index);
        }


        let mut codes = Vec::with_capacity(bodies.len());
        let mut centers = Vec::with_capacity(bodies.len());
        let mut states = Vec::with_capacity(bodies.len());
        for (idx, spk_index) in indices.iter().zip(spk_indices.iter()) {
            let spk = &self.spk[*spk_index];
            let target_map = &spk.target_map;
            let target_index = target_map.get(idx)
                .ok_or_else(|| SpiceError::ParseError(format!("Body {} not found in target map", idx)))?;
            let target = &spk.targets[*target_index];

            let state = spk.state_at(epoch, *target_index)?;
            let state = (state.0, state.1, state.2, state.3, state.4, state.5);
    
            states.push(state);
            codes.push(target.code);
            centers.push(target.cen);
        }
        
        // find all unique non-zero centers
        let mut unique_centers = HashSet::new();
        for center in &centers {
            if *center != 0 {
                unique_centers.insert(*center);
            }
        }

        let mut center_states = HashMap::new();
        let ssb = SpiceBody::new(0, "SSB", 0.0);
        for c in unique_centers {
            let state = self.compute_state(&SpiceBody::new(c, "oops", 0.0), &ssb, epoch)?;
            center_states.insert(c, state);
        }

        for idx in 0..states.len() {
            let center = centers[idx];
            if center == 0 {
                continue;
            }

            let correction = center_states.get(&center)
                .ok_or_else(|| SpiceError::ParseError(format!("Center {} not found in center states", center)))?;


            let state = states[idx];
            states[idx] = (
                state.0 + correction.0,
                state.1 + correction.1,
                state.2 + correction.2,
                state.3 + correction.3,
                state.4 + correction.4,
                state.5 + correction.5,
            );
        }

        Ok((states))
        
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
