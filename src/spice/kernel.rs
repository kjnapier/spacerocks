//! `SpiceKernel`: a self-contained, thread-safe collection of loaded kernels.
//!
//! Unlike the CSPICE toolkit, which keeps a single global kernel pool, every `SpiceKernel`
//! is its own independent context. Kernel files are reference counted, so cloning a
//! `SpiceKernel` is cheap and the clones share the underlying (memory-mapped) data:
//!
//! ```no_run
//! use spacerocks::spice::SpiceKernel;
//! let mut base = SpiceKernel::new();
//! base.load("de440s.bsp").unwrap();
//! let mut with_jwst = base.clone();          // shares de440s.bsp with `base`
//! with_jwst.load("jwst_pred.bsp").unwrap();  // only affects `with_jwst`
//! ```
//!
//! Queries never mutate the kernel, so a `&SpiceKernel` (or `Arc<SpiceKernel>`) can be
//! used concurrently from many threads without locking.
//!
//! Precedence follows SPICE: data from files loaded later take priority over earlier files,
//! and within one file, later segments take priority over earlier ones.

use std::collections::HashMap;
use std::hash::{BuildHasherDefault, Hasher};
use std::path::{Path, PathBuf};
use std::sync::Arc;

use super::bodies::{builtin_body_id, builtin_body_name, normalize_name, SpiceBody};
use super::daf::Daf;
use super::error::{Result, SpiceError};
use super::frames::{self, FrameInfo, StateTransform, J2000};
use super::pck::PckSegment;
use super::spk::SpkSegment;
use super::text::KernelPool;

/// Kilometers per astronomical unit (IAU 2012).
pub const AU_KM: f64 = 149_597_870.7;
/// Seconds per day.
pub const SECONDS_PER_DAY: f64 = 86_400.0;
/// Julian date of the J2000 epoch (TDB).
pub const J2000_JD: f64 = 2_451_545.0;

/// Convert a TDB Julian date to ephemeris time (TDB seconds past J2000).
/// 1962 Jan 20 00:00 UTC as ephemeris time (SPICE `str2et` with the NAIF leap-second kernel):
/// where sorcha, and so layup, switch from IAU_EARTH to ITRF93 for observatory positions.
pub const ITRF93_START_ET: f64 = -1197547158.8155186;

#[inline(always)]
pub fn et_from_jd(jd_tdb: f64) -> f64 {
    (jd_tdb - J2000_JD) * SECONDS_PER_DAY
}

/// Convert ephemeris time (TDB seconds past J2000) to a TDB Julian date.
#[inline(always)]
pub fn jd_from_et(et: f64) -> f64 {
    et / SECONDS_PER_DAY + J2000_JD
}

/// Convert a km, km/s state to AU, AU/day.
#[inline(always)]
pub fn km_to_au_state(s: &[f64; 6]) -> [f64; 6] {
    let v = SECONDS_PER_DAY / AU_KM;
    [s[0] / AU_KM, s[1] / AU_KM, s[2] / AU_KM, s[3] * v, s[4] * v, s[5] * v]
}

/// Canonical form of a kernel path, used to identify loaded files.
fn canonical(path: &Path) -> String {
    std::fs::canonicalize(path)
        .map(|p| p.display().to_string())
        .unwrap_or_else(|_| path.display().to_string())
}

/// Minimal, fast hasher for integer IDs.
#[derive(Default, Clone, Copy)]
pub(crate) struct IdHasher(u64);

impl Hasher for IdHasher {
    #[inline]
    fn finish(&self) -> u64 {
        self.0
    }
    #[inline]
    fn write(&mut self, bytes: &[u8]) {
        for &b in bytes {
            self.0 = (self.0.rotate_left(5) ^ b as u64).wrapping_mul(0x51_7c_c1_b7_27_22_0a_95);
        }
    }
    #[inline]
    fn write_i32(&mut self, i: i32) {
        self.0 = (i as u32 as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15);
    }
}

pub(crate) type IdMap<V> = HashMap<i32, V, BuildHasherDefault<IdHasher>>;

/// The kind of a loaded kernel file.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum KernelKind {
    Spk,
    Pck,
    Text,
    Meta,
}

impl std::fmt::Display for KernelKind {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let s = match self {
            KernelKind::Spk => "SPK",
            KernelKind::Pck => "PCK (binary)",
            KernelKind::Text => "text",
            KernelKind::Meta => "meta-kernel",
        };
        write!(f, "{}", s)
    }
}

/// A loaded SPK file.
#[derive(Debug)]
pub struct SpkFile {
    pub path: String,
    pub daf: Daf,
    pub segments: Vec<SpkSegment>,
    /// Constants listed in the comment area under "Initial conditions and constants used for
    /// integration" (JPL planetary ephemerides), as (name, value).
    pub constants: Vec<(String, f64)>,
}

/// Parse the integration constants JPL writes into planetary SPK comment areas.
fn integration_constants(comments: &str) -> Vec<(String, f64)> {
    let Some(start) = comments.find("Initial conditions and constants used for integration") else {
        return Vec::new();
    };
    comments[start..]
        .lines()
        .filter_map(|line| {
            let mut it = line.split_whitespace();
            let key = it.next()?;
            let value = it.next()?.replace(['D', 'd'], "E").parse::<f64>().ok()?;
            if it.next().is_some() || !key.chars().all(|c| c.is_ascii_alphanumeric() || c == '_') {
                return None;
            }
            Some((key.to_string(), value))
        })
        .collect()
}

/// A loaded binary PCK file.
#[derive(Debug)]
pub struct PckFile {
    pub path: String,
    pub daf: Daf,
    pub segments: Vec<PckSegment>,
}

/// A loaded text kernel (the raw text is kept so the pool can be rebuilt on unload).
#[derive(Debug)]
pub struct TextFile {
    pub path: String,
    pub text: String,
}

#[derive(Debug, Clone)]
enum Loaded {
    Spk(Arc<SpkFile>),
    Pck(Arc<PckFile>),
    Text(Arc<TextFile>),
    /// Meta-kernel: its own (non-loader) variables plus the files it caused to be loaded.
    Meta { text: Arc<TextFile>, children: Vec<String> },
}

impl Loaded {
    fn path(&self) -> &str {
        match self {
            Loaded::Spk(f) => &f.path,
            Loaded::Pck(f) => &f.path,
            Loaded::Text(f) => &f.path,
            Loaded::Meta { text, .. } => &text.path,
        }
    }

    fn kind(&self) -> KernelKind {
        match self {
            Loaded::Spk(_) => KernelKind::Spk,
            Loaded::Pck(_) => KernelKind::Pck,
            Loaded::Text(_) => KernelKind::Text,
            Loaded::Meta { .. } => KernelKind::Meta,
        }
    }
}

/// Information about one loaded kernel file.
#[derive(Debug, Clone, PartialEq)]
pub struct LoadedKernel {
    pub path: String,
    pub kind: KernelKind,
}

/// What one loaded kernel file contains, for inspection and display (see
/// [`SpiceKernel::summary`]).
#[derive(Debug, Clone, PartialEq)]
pub struct KernelSummary {
    pub path: String,
    pub kind: KernelKind,
    /// File size on disk, if the file can still be read.
    pub size_bytes: Option<u64>,
    /// SPK and binary PCK files: the segments grouped by what they describe, in file order.
    pub groups: Vec<SegmentGroup>,
    /// Text kernels and meta-kernels: the variables the file assigns, sorted.
    pub variables: Vec<String>,
    /// Meta-kernels: the files loaded through it.
    pub children: Vec<String>,
}

/// The segments of one file that share a body, center, frame and data type.
///
/// For SPK files `body` is the target and `center` the center of motion. For binary PCK
/// files `body` is the frame class ID (e.g. 3000 for ITRF93), `frame` is the inertial frame
/// the orientation is given relative to, and `center` is 0.
#[derive(Debug, Clone, PartialEq)]
pub struct SegmentGroup {
    pub body: i32,
    pub center: i32,
    pub frame: i32,
    pub data_type: i32,
    pub n_segments: usize,
    /// Merged ET coverage of the group's segments.
    pub coverage: Vec<Interval>,
}

fn group_segments(items: impl Iterator<Item = (i32, i32, i32, i32, f64, f64)>) -> Vec<SegmentGroup> {
    let mut order: Vec<(i32, i32, i32, i32)> = Vec::new();
    let mut acc: HashMap<(i32, i32, i32, i32), (usize, Vec<Interval>)> = HashMap::new();
    for (body, center, frame, ty, start, end) in items {
        let key = (body, center, frame, ty);
        let e = acc.entry(key).or_insert_with(|| {
            order.push(key);
            (0, Vec::new())
        });
        e.0 += 1;
        e.1.push(Interval { start, end });
    }
    order
        .into_iter()
        .map(|key| {
            let (n, mut iv) = acc.remove(&key).unwrap();
            SegmentGroup {
                body: key.0,
                center: key.1,
                frame: key.2,
                data_type: key.3,
                n_segments: n,
                coverage: merge_intervals(&mut iv),
            }
        })
        .collect()
}

fn pool_variable_names(text: &TextFile) -> Vec<String> {
    let mut scratch = KernelPool::new();
    if scratch.load_str(&text.text, &text.path).is_err() {
        return Vec::new();
    }
    let mut v: Vec<String> = scratch.names().cloned().collect();
    v.sort();
    v
}

/// A time interval (ET seconds) during which a body has ephemeris data.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Interval {
    pub start: f64,
    pub end: f64,
}

/// A collection of loaded SPICE kernels and everything derived from them.
#[derive(Clone, Default)]
pub struct SpiceKernel {
    loaded: Vec<Loaded>,
    // Derived state, rebuilt whenever the set of loaded files changes.
    pub(crate) spk_files: Vec<Arc<SpkFile>>,
    pub(crate) pck_files: Vec<Arc<PckFile>>,
    /// target -> [(spk file slot, segment index)], highest priority first.
    spk_index: IdMap<Vec<(u32, u32)>>,
    /// frame class ID -> [(pck file slot, segment index)], highest priority first.
    pck_index: IdMap<Vec<(u32, u32)>>,
    pub(crate) pool: KernelPool,
    names_to_codes: HashMap<String, i32>,
    codes_to_names: IdMap<String>,
    /// Frames defined in loaded text kernels (FRAME_<id>_NAME/CLASS/...).
    pub(crate) kernel_frames: IdMap<FrameInfo>,
    /// Frame names defined in loaded text kernels (FRAME_<name> = id).
    pub(crate) kernel_frame_names: HashMap<String, i32>,
}

impl std::fmt::Debug for SpiceKernel {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("SpiceKernel")
            .field("loaded", &self.loaded_kernels())
            .field("spk_segments", &self.spk_files.iter().map(|f| f.segments.len()).sum::<usize>())
            .field("pck_segments", &self.pck_files.iter().map(|f| f.segments.len()).sum::<usize>())
            .field("pool_variables", &self.pool.len())
            .finish()
    }
}

impl std::fmt::Display for SpiceKernel {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "SpiceKernel with {} loaded file(s):", self.loaded.len())?;
        for k in self.loaded_kernels() {
            writeln!(f, "  [{}] {}", k.kind, k.path)?;
        }
        Ok(())
    }
}

const MAX_CHAIN: usize = 32;

impl SpiceKernel {
    /// An empty kernel context.
    pub fn new() -> Self {
        Self::default()
    }

    // ---------------------------------------------------------------------------------------
    // Loading
    // ---------------------------------------------------------------------------------------

    /// Load a kernel, detecting its type from the file contents: SPK and binary PCK files
    /// (DAF), text kernels (LSK, text PCK, FK, ...) and meta-kernels.
    ///
    /// Loading a file that is already loaded moves it to the highest priority, as in SPICE.
    /// Files are identified by their canonical path.
    pub fn load<P: AsRef<Path>>(&mut self, path: P) -> Result<()> {
        let mut stack = Vec::new();
        let result = self.load_inner(path.as_ref(), &mut stack);
        // Rebuild even on failure so the derived state always matches `loaded`.
        self.rebuild();
        result
    }

    /// A kernel with the standard spacerocks kernel set loaded, downloading any missing
    /// files into the kernel cache (`~/.spacerocks/spice`, or `$SPACEROCKS_SPICE_DIR`).
    ///
    /// The set is: `latest_leapseconds.tls`, `de440s.bsp` (planets, 1849–2150), the newest
    /// `earth_1962_*_combined.bpc` (high-precision Earth orientation, history + prediction),
    /// `gm_de440.tpc` (GM values), and `sb441-n16.bsp` (the 16 asteroid perturbers used by
    /// the n-body integrator; about 650 MB on first download). See
    /// [`config::default_kernel_list`](super::config::default_kernel_list).
    pub fn defaults() -> Result<Self> {
        Self::from_config(&super::config::Config::default())
    }

    /// Like [`defaults`](Self::defaults), optionally without downloading (only files already
    /// in the cache or search paths are used; missing ones are an error).
    pub fn defaults_with_download(download: bool) -> Result<Self> {
        Self::from_config(&super::config::Config::default_with_download(download))
    }

    /// Build a kernel from a kernel-set configuration, downloading missing files if the
    /// configuration allows it.
    pub fn from_config(config: &super::config::Config) -> Result<Self> {
        let mut k = SpiceKernel::new();
        k.load_config(config)?;
        Ok(k)
    }

    /// Build a kernel from a TOML kernel-set file (see [`config`](super::config)).
    pub fn from_config_file<P: AsRef<Path>>(path: P) -> Result<Self> {
        Self::from_config(&super::config::Config::from_file(path)?)
    }

    /// Load every kernel of a kernel-set configuration into this kernel, in order.
    pub fn load_config(&mut self, config: &super::config::Config) -> Result<()> {
        for spec in &config.default_kernels {
            let (path, fresh) = config.resolve(spec)?;
            if let Err(e) = self.load(&path) {
                if fresh {
                    // A corrupt or non-kernel download (e.g. an HTML error page): don't keep it.
                    let _ = std::fs::remove_file(&path);
                }
                return Err(e);
            }
        }
        Ok(())
    }

    /// Load an SPK file.
    pub fn load_spk<P: AsRef<Path>>(&mut self, path: P) -> Result<()> {
        let key = canonical(path.as_ref());
        let daf = Daf::open(path)?;
        let entry = Self::spk_entry(daf, key)?;
        self.insert(entry);
        self.rebuild();
        Ok(())
    }

    /// Load a binary PCK file.
    pub fn load_pck<P: AsRef<Path>>(&mut self, path: P) -> Result<()> {
        let key = canonical(path.as_ref());
        let daf = Daf::open(path)?;
        let entry = Self::pck_entry(daf, key)?;
        self.insert(entry);
        self.rebuild();
        Ok(())
    }

    /// Alias for [`load_pck`](Self::load_pck) (binary PCK files usually end in `.bpc`).
    pub fn load_bpc<P: AsRef<Path>>(&mut self, path: P) -> Result<()> {
        self.load_pck(path)
    }

    /// Load a text kernel or meta-kernel.
    pub fn load_text<P: AsRef<Path>>(&mut self, path: P) -> Result<()> {
        let mut stack = Vec::new();
        let result = self.load_text_inner(path.as_ref(), &mut stack);
        self.rebuild();
        result
    }

    fn load_inner(&mut self, path: &Path, stack: &mut Vec<String>) -> Result<()> {
        let path_str = path.display().to_string();
        let mut head = [b' '; 8];
        {
            use std::io::Read;
            let mut f = std::fs::File::open(path).map_err(|e| SpiceError::io(&path_str, e))?;
            let mut n = 0;
            while n < 8 {
                let k = f.read(&mut head[n..]).map_err(|e| SpiceError::io(&path_str, e))?;
                if k == 0 {
                    break;
                }
                n += k;
            }
        }
        let id = String::from_utf8_lossy(&head).trim().to_string();
        let key = canonical(path);
        if id.starts_with("DAF/SPK") {
            let entry = Self::spk_entry(Daf::open(path)?, key)?;
            self.insert(entry);
            Ok(())
        } else if id.starts_with("DAF/PCK") {
            let entry = Self::pck_entry(Daf::open(path)?, key)?;
            self.insert(entry);
            Ok(())
        } else if id == "NAIF/DAF" {
            // Pre-1995 DAFs carry no type in the ID word; decide from the summary format.
            let daf = Daf::open(path)?;
            let entry = if daf.nd == 2 && daf.ni == 5 {
                Self::pck_entry(daf, key)?
            } else {
                Self::spk_entry(daf, key)?
            };
            self.insert(entry);
            Ok(())
        } else if id.starts_with("DAF/") {
            Err(SpiceError::invalid(&path_str, format!("unsupported DAF kernel type '{}'", id)))
        } else if id.starts_with("DAS/") {
            Err(SpiceError::invalid(&path_str, format!("DAS kernels ('{}') are not supported", id)))
        } else {
            self.load_text_inner(path, stack)
        }
    }

    fn spk_entry(daf: Daf, key: String) -> Result<Loaded> {
        if daf.nd != 2 || daf.ni != 6 {
            return Err(SpiceError::invalid(
                &daf.path,
                format!("SPK files have ND=2, NI=6; this file has ND={}, NI={}", daf.nd, daf.ni),
            ));
        }
        let words = daf.words();
        let segments = daf
            .summaries
            .iter()
            .map(|s| SpkSegment::from_summary(s, words))
            .collect::<Result<Vec<_>>>()?;
        let constants = integration_constants(&daf.comments());
        Ok(Loaded::Spk(Arc::new(SpkFile { path: key, daf, segments, constants })))
    }

    fn pck_entry(daf: Daf, key: String) -> Result<Loaded> {
        if daf.nd != 2 || daf.ni != 5 {
            return Err(SpiceError::invalid(
                &daf.path,
                format!("binary PCK files have ND=2, NI=5; this file has ND={}, NI={}", daf.nd, daf.ni),
            ));
        }
        let words = daf.words();
        let segments = daf
            .summaries
            .iter()
            .map(|s| PckSegment::from_summary(s, words))
            .collect::<Result<Vec<_>>>()?;
        Ok(Loaded::Pck(Arc::new(PckFile { path: key, daf, segments })))
    }

    /// Add a loaded file at the highest priority, replacing any previous copy of it.
    fn insert(&mut self, entry: Loaded) {
        let key = entry.path().to_string();
        self.claim(&key);
        self.loaded.push(entry);
    }

    /// Remove a previous copy of `key` and detach it from any meta-kernel that loaded it, so
    /// that unloading that meta-kernel later does not remove the new copy.
    fn claim(&mut self, key: &str) {
        self.remove_path(key);
        for l in self.loaded.iter_mut() {
            if let Loaded::Meta { children, .. } = l {
                children.retain(|c| c != key);
            }
        }
    }

    fn load_text_inner(&mut self, path: &Path, stack: &mut Vec<String>) -> Result<()> {
        let path_str = path.display().to_string();
        let bytes = std::fs::read(path).map_err(|e| SpiceError::io(&path_str, e))?;
        let text = String::from_utf8_lossy(&bytes).into_owned();
        if !text.contains("\\begindata") {
            return Err(SpiceError::invalid(
                &path_str,
                "not a recognized kernel (no DAF header and no \\begindata section)",
            ));
        }
        // Parse into a scratch pool to validate and to detect meta-kernels.
        let mut scratch = KernelPool::new();
        scratch.load_str(&text, &path_str)?;
        let key = canonical(path);
        if scratch.contains("KERNELS_TO_LOAD") {
            return self.load_meta(path, key, text, &scratch, stack);
        }
        self.insert(Loaded::Text(Arc::new(TextFile { path: key, text })));
        Ok(())
    }

    fn load_meta(&mut self, path: &Path, key: String, text: String, pool: &KernelPool, stack: &mut Vec<String>) -> Result<()> {
        const MAX_META_DEPTH: usize = 16;
        if stack.contains(&key) {
            return Err(SpiceError::invalid(&key, "meta-kernel loads itself (directly or through other meta-kernels)"));
        }
        if stack.len() >= MAX_META_DEPTH {
            return Err(SpiceError::invalid(&key, "meta-kernels nested too deeply"));
        }
        let files = pool.get_strs_joined("KERNELS_TO_LOAD").unwrap_or_default();
        let symbols = pool.get_strs("PATH_SYMBOLS").map(|v| v.to_vec()).unwrap_or_default();
        let values = pool.get_strs_joined("PATH_VALUES").unwrap_or_default();
        if symbols.len() != values.len() {
            return Err(SpiceError::PoolVariable {
                name: "PATH_SYMBOLS".into(),
                reason: format!("{} symbols but {} PATH_VALUES in {}", symbols.len(), values.len(), key),
            });
        }
        let base_dir = path.parent().map(|p| p.to_path_buf()).unwrap_or_else(|| PathBuf::from("."));

        // As in SPICE, the meta-kernel's own variables enter the pool before the kernels it
        // loads, so it is inserted ahead of its children.
        self.insert(Loaded::Meta {
            text: Arc::new(TextFile { path: key.clone(), text }),
            children: Vec::new(),
        });
        stack.push(key.clone());
        // Longest symbols first so $AB is not clobbered by $A.
        let mut order: Vec<usize> = (0..symbols.len()).collect();
        order.sort_by_key(|&i| std::cmp::Reverse(symbols[i].len()));
        let mut result = Ok(());
        for f in files {
            let mut resolved = f.clone();
            for &i in &order {
                resolved = resolved.replace(&format!("${}", symbols[i].trim()), values[i].trim());
            }
            let mut p = PathBuf::from(&resolved);
            if p.is_relative() && !p.exists() {
                // SPICE resolves relative paths against the working directory; as a convenience
                // we also try the meta-kernel's own directory.
                let alt = base_dir.join(&p);
                if alt.exists() {
                    p = alt;
                }
            }
            if let Err(e) = self.load_inner(&p, stack) {
                // Roll back: drop the partially loaded child and this meta-kernel (with the
                // children it already loaded), so a failed load leaves no partial state.
                self.remove_path(&canonical(&p));
                self.remove_path(&key);
                result = Err(e);
                break;
            }
            let child = canonical(&p);
            if let Some(Loaded::Meta { children, .. }) = self.loaded.iter_mut().find(|l| l.path() == key) {
                children.push(child);
            }
        }
        stack.pop();
        result
    }

    /// Unload a kernel (for a meta-kernel, also the files it loaded). Returns whether
    /// anything was unloaded.
    pub fn unload<P: AsRef<Path>>(&mut self, path: P) -> bool {
        let key = canonical(path.as_ref());
        let found = self.remove_path(&key);
        if found {
            self.rebuild();
        }
        found
    }

    /// Unload everything.
    pub fn clear(&mut self) {
        *self = SpiceKernel::new();
    }

    fn remove_path(&mut self, path: &str) -> bool {
        let mut removed = false;
        let mut children_to_remove = Vec::new();
        self.loaded.retain(|l| {
            if l.path() == path {
                if let Loaded::Meta { children, .. } = l {
                    children_to_remove.extend(children.iter().cloned());
                }
                removed = true;
                false
            } else {
                true
            }
        });
        for c in children_to_remove {
            self.remove_path(&c);
        }
        removed
    }

    /// Files currently loaded, in load order (lowest priority first).
    pub fn loaded_kernels(&self) -> Vec<LoadedKernel> {
        self.loaded
            .iter()
            .map(|l| LoadedKernel {
                path: l.path().to_string(),
                kind: l.kind(),
            })
            .collect()
    }

    /// What each loaded file contains, in load order (lowest priority first).
    pub fn summary(&self) -> Vec<KernelSummary> {
        self.loaded
            .iter()
            .map(|l| {
                let path = l.path().to_string();
                let size_bytes = std::fs::metadata(&path).ok().map(|m| m.len());
                let (groups, variables, children) = match l {
                    Loaded::Spk(f) => (
                        group_segments(f.segments.iter().map(|s| (s.target, s.center, s.frame, s.data_type, s.start_et, s.end_et))),
                        Vec::new(),
                        Vec::new(),
                    ),
                    Loaded::Pck(f) => (
                        group_segments(f.segments.iter().map(|s| (s.class_id, 0, s.reference_frame, s.data_type, s.start_et, s.end_et))),
                        Vec::new(),
                        Vec::new(),
                    ),
                    Loaded::Text(t) => (Vec::new(), pool_variable_names(t), Vec::new()),
                    Loaded::Meta { text, children } => (Vec::new(), pool_variable_names(text), children.clone()),
                };
                KernelSummary { path, kind: l.kind(), size_bytes, groups, variables, children }
            })
            .collect()
    }

    /// Rebuild indexes and the kernel pool from the loaded files.
    fn rebuild(&mut self) {
        self.spk_files.clear();
        self.pck_files.clear();
        self.spk_index.clear();
        self.pck_index.clear();
        self.pool = KernelPool::new();
        for l in &self.loaded {
            match l {
                Loaded::Spk(f) => self.spk_files.push(f.clone()),
                Loaded::Pck(f) => self.pck_files.push(f.clone()),
                Loaded::Text(t) => {
                    // Already validated at load time.
                    let _ = self.pool.load_str(&t.text, &t.path);
                }
                Loaded::Meta { text, .. } => {
                    let mut scratch = KernelPool::new();
                    if scratch.load_str(&text.text, &text.path).is_ok() {
                        let names: Vec<String> = scratch.names().cloned().collect();
                        for n in names {
                            if matches!(n.as_str(), "KERNELS_TO_LOAD" | "PATH_SYMBOLS" | "PATH_VALUES") {
                                continue;
                            }
                            match scratch.get(&n).cloned() {
                                Some(super::text::PoolValue::Numeric(v)) => self.pool.set_f64s(&n, v),
                                Some(super::text::PoolValue::Strings(v)) => self.pool.set_strs(&n, v),
                                None => {}
                            }
                        }
                    }
                }
            }
        }
        // Highest priority first: iterate files from last loaded, segments from last in file.
        for (slot, f) in self.spk_files.iter().enumerate().rev() {
            for (i, s) in f.segments.iter().enumerate().rev() {
                self.spk_index.entry(s.target).or_default().push((slot as u32, i as u32));
            }
        }
        for (slot, f) in self.pck_files.iter().enumerate().rev() {
            for (i, s) in f.segments.iter().enumerate().rev() {
                self.pck_index.entry(s.class_id).or_default().push((slot as u32, i as u32));
            }
        }
        // Body name mappings from the pool (NAIF_BODY_NAME / NAIF_BODY_CODE).
        self.names_to_codes.clear();
        self.codes_to_names.clear();
        if let (Some(names), Some(codes)) = (self.pool.get_strs("NAIF_BODY_NAME"), self.pool.get_f64s("NAIF_BODY_CODE")) {
            for (n, c) in names.iter().zip(codes.iter()) {
                let code = c.round() as i32;
                self.names_to_codes.insert(normalize_name(n), code);
                self.codes_to_names.insert(code, n.trim().to_string());
            }
        }
        self.rebuild_frames();
    }

    /// The kernel variable pool built from loaded text kernels.
    pub fn pool(&self) -> &KernelPool {
        &self.pool
    }

    // ---------------------------------------------------------------------------------------
    // Bodies
    // ---------------------------------------------------------------------------------------

    /// NAIF ID for a body name. Names defined in loaded text kernels take precedence over
    /// the built-in table; integer strings are accepted as IDs.
    pub fn body_id(&self, name: &str) -> Option<i32> {
        let n = normalize_name(name);
        if let Some(&c) = self.names_to_codes.get(&n) {
            return Some(c);
        }
        builtin_body_id(&n)
    }

    /// Preferred name for a NAIF ID.
    pub fn body_name(&self, code: i32) -> Option<String> {
        if let Some(n) = self.codes_to_names.get(&code) {
            return Some(n.clone());
        }
        builtin_body_name(code).map(|s| s.to_string())
    }

    /// Resolve a body name to a [`SpiceBody`], including its GM if a text PCK provides it.
    pub fn body(&self, name: &str) -> Result<SpiceBody> {
        let code = self.body_id(name).ok_or_else(|| SpiceError::UnknownBody(name.to_string()))?;
        Ok(SpiceBody::new(code, &normalize_name(name), self.gm(code).unwrap_or(0.0)))
    }

    /// GM (km^3/s^2) of a body from `BODYnnn_GM` in the kernel pool.
    pub fn gm(&self, code: i32) -> Option<f64> {
        self.pool.get_f64(&format!("BODY{}_GM", code))
    }

    // ---------------------------------------------------------------------------------------
    // Ephemeris
    // ---------------------------------------------------------------------------------------

    /// Highest-priority SPK segment for `body` covering `et`.
    #[inline]
    fn find_spk(&self, body: i32, et: f64) -> Option<(&SpkFile, &SpkSegment)> {
        let list = self.spk_index.get(&body)?;
        for &(slot, i) in list {
            let f = &self.spk_files[slot as usize];
            let s = &f.segments[i as usize];
            if s.covers(et) {
                return Some((f, s));
            }
        }
        None
    }

    /// Highest-priority PCK segment for frame class ID covering `et`: (file slot, segment).
    #[inline]
    pub(crate) fn find_pck(&self, class_id: i32, et: f64) -> Option<(usize, usize)> {
        let list = self.pck_index.get(&class_id)?;
        for &(slot, i) in list {
            let s = &self.pck_files[slot as usize].segments[i as usize];
            if s.covers(et) {
                return Some((slot as usize, i as usize));
            }
        }
        None
    }

    /// Evaluate one segment and express the result in J2000.
    #[inline]
    fn eval_segment_j2000(&self, f: &SpkFile, s: &SpkSegment, et: f64) -> Result<[f64; 6]> {
        let st = s.state(f.daf.words(), et)?;
        if s.frame == J2000 {
            return Ok(st);
        }
        let xf = self.frame_to_j2000(s.frame, et)?;
        Ok(xf.apply(&st))
    }

    /// Geometric state (km, km/s) of `target` relative to `observer` in J2000 at `et`
    /// (TDB seconds past J2000). Equivalent to SPICE `spkgeo(target, et, "J2000", observer)`.
    pub fn spkgeo(&self, target: i32, observer: i32, et: f64) -> Result<[f64; 6]> {
        if target == observer {
            return Ok([0.0; 6]);
        }
        // Chain from the target toward the root, with cumulative states.
        let mut nodes = [0i32; MAX_CHAIN];
        let mut cum = [[0.0f64; 6]; MAX_CHAIN];
        nodes[0] = target;
        let mut n = 1;
        let mut body = target;
        let mut acc = [0.0f64; 6];
        while let Some((f, s)) = self.find_spk(body, et) {
            let st = self.eval_segment_j2000(f, s, et)?;
            for k in 0..6 {
                acc[k] += st[k];
            }
            body = s.center;
            if body == observer {
                return Ok(acc);
            }
            if n == MAX_CHAIN {
                return Err(self.insufficient(target, observer, et, "ephemeris chain too long"));
            }
            nodes[n] = body;
            cum[n] = acc;
            n += 1;
        }
        // Walk up from the observer until we meet the target's chain.
        let mut oacc = [0.0f64; 6];
        let mut body = observer;
        for _ in 0..MAX_CHAIN {
            if let Some(k) = nodes[..n].iter().position(|&b| b == body) {
                let t = &cum[k];
                return Ok([
                    t[0] - oacc[0],
                    t[1] - oacc[1],
                    t[2] - oacc[2],
                    t[3] - oacc[3],
                    t[4] - oacc[4],
                    t[5] - oacc[5],
                ]);
            }
            match self.find_spk(body, et) {
                Some((f, s)) => {
                    let st = self.eval_segment_j2000(f, s, et)?;
                    for k in 0..6 {
                        oacc[k] += st[k];
                    }
                    body = s.center;
                }
                None => break,
            }
        }
        Err(self.insufficient(target, observer, et, ""))
    }

    fn insufficient(&self, target: i32, observer: i32, et: f64, extra: &str) -> SpiceError {
        let mut detail = String::new();
        for (id, label) in [(target, "target"), (observer, "observer")] {
            if id == 0 {
                continue; // the SSB is the root of the ephemeris tree; it never has data of its own
            }
            match self.spk_index.get(&id) {
                None => detail.push_str(&format!("; no SPK data loaded for {} {}", label, id)),
                Some(list) => {
                    if !list.iter().any(|&(f, i)| self.spk_files[f as usize].segments[i as usize].covers(et)) {
                        detail.push_str(&format!("; SPK data for {} {} does not cover this epoch", label, id));
                    }
                }
            }
        }
        if !extra.is_empty() {
            detail.push_str("; ");
            detail.push_str(extra);
        }
        SpiceError::InsufficientEphemerisData {
            target,
            observer,
            et,
            detail,
        }
    }

    /// Geometric state (km, km/s) of `target` relative to `observer` in `frame` at `et`.
    pub fn state(&self, target: i32, observer: i32, frame: &str, et: f64) -> Result<[f64; 6]> {
        let st = self.spkgeo(target, observer, et)?;
        let fid = self.frame(frame)?.id;
        if fid == J2000 {
            return Ok(st);
        }
        let xf = self.frame_to_j2000(fid, et)?.inverse();
        Ok(xf.apply(&st))
    }

    /// Like [`state`](Self::state) but with body names.
    pub fn state_by_name(&self, target: &str, observer: &str, frame: &str, et: f64) -> Result<[f64; 6]> {
        let t = self.body_id(target).ok_or_else(|| SpiceError::UnknownBody(target.to_string()))?;
        let o = self.body_id(observer).ok_or_else(|| SpiceError::UnknownBody(observer.to_string()))?;
        self.state(t, o, frame, et)
    }

    /// State of `target` relative to `observer` in J2000, in AU and AU/day, at a TDB Julian date.
    pub fn state_au(&self, target: i32, observer: i32, jd_tdb: f64) -> Result<[f64; 6]> {
        Ok(km_to_au_state(&self.spkgeo(target, observer, et_from_jd(jd_tdb))?))
    }

    /// Barycentric (SSB-relative) J2000 states of several bodies, in AU and AU/day, written into
    /// `out`. Shared intermediate nodes (e.g. the Earth-Moon barycenter for the Earth and the
    /// Moon) are evaluated once per call.
    pub fn barycentric_states_au(&self, ids: &[i32], jd_tdb: f64, out: &mut [[f64; 6]]) -> Result<()> {
        self.barycentric_states_au_et(ids, et_from_jd(jd_tdb), out)
    }

    /// [`SpiceKernel::barycentric_states_au`] at `jd_ref + dt` (TDB Julian date plus days),
    /// without rounding the sum to a Julian date first. A Julian date near 2.46e6 resolves only
    /// ~40 µs; ephemeris seconds resolve ~0.1 µs, which matters for fast-moving perturbers.
    pub fn barycentric_states_au_rel(&self, ids: &[i32], jd_ref: f64, dt: f64, out: &mut [[f64; 6]]) -> Result<()> {
        self.barycentric_states_au_et(ids, et_from_jd(jd_ref) + dt * 86400.0, out)
    }

    /// [`SpiceKernel::barycentric_states_au`] at ephemeris time `et` (TDB seconds past J2000).
    pub fn barycentric_states_au_et(&self, ids: &[i32], et: f64, out: &mut [[f64; 6]]) -> Result<()> {
        if out.len() < ids.len() {
            return Err(SpiceError::InvalidFile {
                path: String::new(),
                reason: format!("output buffer holds {} states but {} bodies were requested", out.len(), ids.len()),
            });
        }
        // Small memo of (node, state wrt SSB in km).
        let mut memo_ids = [0i32; 64];
        let mut memo_states = [[0.0f64; 6]; 64];
        let mut nmemo = 1usize; // node 0 (SSB) has zero state
        let mut chain_nodes = [0i32; MAX_CHAIN];
        let mut chain_states = [[0.0f64; 6]; MAX_CHAIN];

        for (idx, &id) in ids.iter().enumerate() {
            // Walk up until we reach a memoized node.
            let mut body = id;
            let mut depth = 0usize;
            let mut base: Option<[f64; 6]> = None;
            loop {
                if let Some(k) = memo_ids[..nmemo].iter().position(|&b| b == body) {
                    base = Some(memo_states[k]);
                    break;
                }
                if depth == MAX_CHAIN {
                    break;
                }
                match self.find_spk(body, et) {
                    Some((f, s)) => {
                        chain_nodes[depth] = body;
                        chain_states[depth] = self.eval_segment_j2000(f, s, et)?;
                        depth += 1;
                        body = s.center;
                    }
                    None => break,
                }
            }
            let result = match base {
                Some(b) => {
                    // Accumulate from the top of the chain down, memoizing intermediate nodes.
                    let mut acc = b;
                    for d in (0..depth).rev() {
                        for k in 0..6 {
                            acc[k] += chain_states[d][k];
                        }
                        if nmemo < memo_ids.len() {
                            memo_ids[nmemo] = chain_nodes[d];
                            memo_states[nmemo] = acc;
                            nmemo += 1;
                        }
                    }
                    acc
                }
                // Could not reach the SSB directly (e.g. data relative to a body whose own
                // ephemeris is loaded "below" it); fall back to the general algorithm.
                None => self.spkgeo(id, 0, et)?,
            };
            out[idx] = km_to_au_state(&result);
        }
        Ok(())
    }

    /// Comment areas of the loaded SPK files as (path, text), highest priority first.
    pub fn spk_comments(&self) -> Vec<(String, String)> {
        self.spk_files.iter().rev().map(|f| (f.path.clone(), f.daf.comments())).collect()
    }

    /// An integration constant (e.g. "AU", "CLIGHT", "J2E", "GMS") recorded in the comment
    /// area of a loaded JPL planetary SPK file, from the highest-priority file that has it.
    pub fn integration_constant(&self, name: &str) -> Option<f64> {
        self.spk_files
            .iter()
            .rev()
            .find_map(|f| f.constants.iter().find(|(k, _)| k == name).map(|(_, v)| *v))
    }

    /// All bodies with SPK data, sorted.
    pub fn spk_bodies(&self) -> Vec<i32> {
        let mut v: Vec<i32> = self.spk_index.keys().copied().collect();
        v.sort_unstable();
        v
    }

    /// All SPK segments for a body, highest priority first.
    pub fn spk_segments(&self, body: i32) -> Vec<&SpkSegment> {
        self.spk_index
            .get(&body)
            .map(|l| {
                l.iter()
                    .map(|&(f, i)| &self.spk_files[f as usize].segments[i as usize])
                    .collect()
            })
            .unwrap_or_default()
    }

    /// Merged ET coverage intervals of a body's SPK data.
    pub fn spk_coverage(&self, body: i32) -> Vec<Interval> {
        let mut iv: Vec<Interval> = self
            .spk_segments(body)
            .iter()
            .map(|s| Interval {
                start: s.start_et,
                end: s.end_et,
            })
            .collect();
        merge_intervals(&mut iv)
    }

    /// Merged ET coverage intervals of binary PCK data for a frame class ID.
    pub fn pck_coverage(&self, class_id: i32) -> Vec<Interval> {
        let mut iv: Vec<Interval> = self
            .pck_index
            .get(&class_id)
            .map(|l| {
                l.iter()
                    .map(|&(f, i)| {
                        let s = &self.pck_files[f as usize].segments[i as usize];
                        Interval {
                            start: s.start_et,
                            end: s.end_et,
                        }
                    })
                    .collect()
            })
            .unwrap_or_default();
        merge_intervals(&mut iv)
    }

    /// State transform from body-fixed ITRF93 (Earth) to J2000 at `et`; convenience for
    /// computing observatory positions. Requires a binary Earth PCK.
    pub fn itrf93_to_j2000(&self, et: f64) -> Result<StateTransform> {
        self.frame_to_j2000(frames::ITRF93, et)
    }

    /// State transform from the Earth-fixed frame to J2000 at `et`: ITRF93 (binary Earth PCK)
    /// from 1962 Jan 20 UTC on, where sorcha and layup switch, and before that, where the binary
    /// PCKs have no data, IAU_EARTH (the IAU rotation model, from a text PCK such as
    /// pck00010.tpc) corrected for the Earth's actual rotation ([`iau_earth_rotation_delay`]).
    pub fn earth_fixed_to_j2000(&self, et: f64) -> Result<StateTransform> {
        if et >= ITRF93_START_ET {
            return self.frame_to_j2000(frames::ITRF93, et);
        }
        self.frame_to_j2000(frames::IAU_EARTH, et - iau_earth_rotation_delay(et)).map_err(|e| match e {
            SpiceError::InsufficientOrientationData { frame, et, detail } => SpiceError::InsufficientOrientationData {
                frame,
                et,
                detail: format!("{}; before 1962 Jan 20 the Earth's orientation comes from the IAU rotation model: load pck00010.tpc (SpiceKernel.defaults() includes it)", detail),
            },
            other => other,
        })
    }
}

/// How far (seconds) the IAU_EARTH rotation model runs ahead of the Earth at `et`: evaluating
/// IAU_EARTH at `et` minus this gives the Earth's orientation at `et`.
///
/// IAU_EARTH turns uniformly with TDB, while the Earth turns with UT1, so the model's error in
/// rotation grows with ΔT = TT − UT1. The delay is ΔT (Stephenson, Morrison & Hohenkerk) less
/// `75.025 + 0.5620 (year − 2000)` s, the offset and rate that fit IAU_EARTH to ITRF93 over
/// 1962–2019. With it, IAU_EARTH agrees with ITRF93 to 0.32 km at the surface (rms 0.21 km, the
/// nutation IAU_EARTH omits) over that span; without it the difference reaches 9 km. Before 1962
/// the uncertainty of ΔT adds about 0.46 km per second of it.
pub fn iau_earth_rotation_delay(et: f64) -> f64 {
    let year = 2000.0 + et / (365.25 * 86400.0);
    crate::time::deltat::delta_t_year(year) - 75.025 - 0.5620 * (year - 2000.0)
}

fn merge_intervals(iv: &mut [Interval]) -> Vec<Interval> {
    iv.sort_by(|a, b| a.start.partial_cmp(&b.start).unwrap_or(std::cmp::Ordering::Equal));
    let mut out: Vec<Interval> = Vec::new();
    for i in iv.iter() {
        match out.last_mut() {
            Some(last) if i.start <= last.end => last.end = last.end.max(i.end),
            _ => out.push(*i),
        }
    }
    out
}
