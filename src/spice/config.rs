//! Kernel sets, configuration files, and automatic downloading of missing kernels.
//!
//! [`SpiceKernel::defaults`](super::SpiceKernel::defaults) loads a standard set of kernels
//! (leapseconds, DE440 planets, Earth orientation, GM values, and the 16 asteroid perturbers
//! used by the n-body code), downloading any that are missing into a local cache directory.
//!
//! The cache directory is `~/.spacerocks/spice`, or `$SPACEROCKS_SPICE_DIR` if that is set
//! (useful on clusters where the home directory has a small quota).
//!
//! A kernel set can also be described in a TOML file:
//!
//! ```toml
//! auto_download = true
//! check_for_updates = false
//! download_dir = "~/data/spice"
//! kernel_paths = ["/n/holylfs/shared/spice"]      # searched before download_dir
//!
//! default_kernels = [
//!     { name = "latest_leapseconds.tls", kernel_type = "lsk" },
//!     { name = "de440s.bsp", kernel_type = "spk/planets" },
//!     { name = "earth_1962_*_combined.bpc", kernel_type = "pck" },   # newest match
//!     { name = "sb441-n16.bsp", url = "https://ssd.jpl.nasa.gov/ftp/eph/small_bodies/asteroids_de441/sb441-n16.bsp" },
//!     { name = "my_local_kernel.bsp" },                                # must already exist locally
//! ]
//! ```
//!
//! `kernel_type` is the subdirectory of the NAIF generic kernels server
//! (<https://naif.jpl.nasa.gov/pub/naif/generic_kernels/>). A `*` in `name` selects the
//! lexicographically last matching file, which for NAIF's date-stamped files is the newest.

use std::fs;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};

use super::error::{Result, SpiceError};

/// Base URL of the NAIF generic kernels server.
pub const NAIF_GENERIC_KERNELS: &str = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels";

/// URL of the DE441-consistent ephemeris of the 16 most massive asteroids.
pub const SB441_N16_URL: &str = "https://ssd.jpl.nasa.gov/ftp/eph/small_bodies/asteroids_de441/sb441-n16.bsp";

/// Specification of one kernel file in a kernel set.
#[derive(Debug, Serialize, Deserialize, Clone, PartialEq)]
pub struct KernelSpec {
    /// File name. May contain `*` wildcards; the newest (lexicographically last) match is used.
    pub name: String,
    /// Subdirectory on the NAIF generic kernels server (e.g. "lsk", "spk/planets", "pck").
    #[serde(default)]
    pub kernel_type: String,
    /// Explicit download URL (of the file, or of its directory when `name` is a pattern).
    /// Overrides `kernel_type`.
    #[serde(default)]
    pub url: Option<String>,
}

impl KernelSpec {
    /// A kernel from the NAIF generic kernels server.
    pub fn naif(kernel_type: &str, name: &str) -> Self {
        KernelSpec {
            name: name.to_string(),
            kernel_type: kernel_type.to_string(),
            url: None,
        }
    }

    /// A kernel downloaded from an explicit URL.
    pub fn from_url(name: &str, url: &str) -> Self {
        KernelSpec {
            name: name.to_string(),
            kernel_type: String::new(),
            url: Some(url.to_string()),
        }
    }

    /// A kernel that must already exist in one of the search directories.
    pub fn local(name: &str) -> Self {
        KernelSpec {
            name: name.to_string(),
            kernel_type: String::new(),
            url: None,
        }
    }

    fn is_pattern(&self) -> bool {
        self.name.contains('*')
    }

    /// URL of the directory holding this kernel, if it can be downloaded.
    fn remote_dir(&self, naif_url: &str) -> Option<String> {
        match &self.url {
            Some(u) if self.is_pattern() => Some(u.trim_end_matches('/').to_string()),
            Some(u) => u.rsplit_once('/').map(|(dir, _)| dir.to_string()),
            None if !self.kernel_type.is_empty() => Some(format!(
                "{}/{}",
                naif_url.trim_end_matches('/'),
                self.kernel_type.trim_matches('/')
            )),
            None => None,
        }
    }

    /// URL of the concrete file `file_name` (a match of this spec).
    fn remote_file(&self, naif_url: &str, file_name: &str) -> Option<String> {
        match &self.url {
            Some(u) if !self.is_pattern() => Some(u.clone()),
            _ => self.remote_dir(naif_url).map(|d| format!("{}/{}", d, file_name)),
        }
    }
}

/// Configuration of a kernel set: which kernels, where to look for them, and whether to
/// download missing ones.
#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct Config {
    /// Kernels to load, in order (later kernels take priority).
    #[serde(default = "default_kernel_list")]
    pub default_kernels: Vec<KernelSpec>,
    /// Extra directories searched (in order) before `download_dir`.
    #[serde(default)]
    pub kernel_paths: Vec<PathBuf>,
    /// Download kernels that are not found locally.
    #[serde(default = "default_true")]
    pub auto_download: bool,
    /// For wildcard names, check the server for a newer version even if one exists locally.
    #[serde(default)]
    pub check_for_updates: bool,
    /// Where downloaded kernels are stored.
    #[serde(default = "default_download_dir")]
    pub download_dir: PathBuf,
    /// Base URL of the NAIF generic kernels tree (change to use a mirror).
    #[serde(default = "default_naif_url")]
    pub naif_url: String,
}

fn default_true() -> bool {
    true
}

fn default_naif_url() -> String {
    NAIF_GENERIC_KERNELS.to_string()
}

/// The standard kernel set loaded by [`SpiceKernel::defaults`](super::SpiceKernel::defaults).
pub fn default_kernel_list() -> Vec<KernelSpec> {
    vec![
        KernelSpec::naif("lsk", "latest_leapseconds.tls"),
        KernelSpec::naif("spk/planets", "de440s.bsp"),
        KernelSpec::naif("pck", "earth_1962_*_combined.bpc"),
        KernelSpec::naif("pck", "gm_de440.tpc"),
        KernelSpec::naif("pck", "pck00010.tpc"),
        KernelSpec::from_url("sb441-n16.bsp", SB441_N16_URL),
    ]
}

/// Default kernel cache directory: `$SPACEROCKS_SPICE_DIR`, else `~/.spacerocks/spice`.
pub fn default_download_dir() -> PathBuf {
    if let Some(d) = std::env::var_os("SPACEROCKS_SPICE_DIR") {
        if !d.is_empty() {
            return expand_tilde(Path::new(&d));
        }
    }
    dirs::home_dir()
        .unwrap_or_else(|| PathBuf::from("."))
        .join(".spacerocks")
        .join("spice")
}

fn expand_tilde(p: &Path) -> PathBuf {
    match p.to_str() {
        Some(s) if s == "~" || s.starts_with("~/") => match dirs::home_dir() {
            Some(h) => h.join(s.trim_start_matches('~').trim_start_matches('/')),
            None => p.to_path_buf(),
        },
        _ => p.to_path_buf(),
    }
}

impl Default for Config {
    fn default() -> Self {
        Config {
            default_kernels: default_kernel_list(),
            kernel_paths: Vec::new(),
            auto_download: true,
            check_for_updates: false,
            download_dir: default_download_dir(),
            naif_url: default_naif_url(),
        }
    }
}

impl Config {
    /// The default configuration with downloading enabled or disabled.
    pub fn default_with_download(download: bool) -> Self {
        Config {
            auto_download: download,
            ..Config::default()
        }
    }

    /// Read a configuration from a TOML file. Missing fields take their default values;
    /// `~` is expanded in paths.
    pub fn from_file<P: AsRef<Path>>(path: P) -> Result<Self> {
        let p = path.as_ref().display().to_string();
        let content = fs::read_to_string(path.as_ref()).map_err(|e| SpiceError::io(&p, e))?;
        Self::from_toml(&content).map_err(|e| match e {
            SpiceError::Config(msg) => SpiceError::Config(format!("{}: {}", p, msg)),
            other => other,
        })
    }

    /// Parse a configuration from TOML text.
    pub fn from_toml(text: &str) -> Result<Self> {
        let mut c: Config = toml::from_str(text).map_err(|e| SpiceError::Config(e.to_string()))?;
        c.download_dir = expand_tilde(&c.download_dir);
        c.kernel_paths = c.kernel_paths.iter().map(|p| expand_tilde(p)).collect();
        Ok(c)
    }

    /// Directories searched for kernels, in order.
    pub fn search_dirs(&self) -> Vec<PathBuf> {
        let mut v = self.kernel_paths.clone();
        v.push(self.download_dir.clone());
        v
    }

    /// Find `spec` locally, downloading it if allowed and necessary. Returns the local path and
    /// whether it was just downloaded.
    pub fn resolve(&self, spec: &KernelSpec) -> Result<(PathBuf, bool)> {
        let dirs = self.search_dirs();
        let local = find_local(&dirs, &spec.name);
        let can_download = self.auto_download && spec.remote_dir(&self.naif_url).is_some();

        if !spec.is_pattern() {
            if let Some(p) = local {
                return Ok((p, false));
            }
            if !can_download {
                return Err(not_found(spec, &dirs, self.auto_download));
            }
            let url = spec.remote_file(&self.naif_url, &spec.name).expect("checked above");
            let dest = self.download_dir.join(&spec.name);
            download(&url, &dest)?;
            return Ok((dest, true));
        }

        // Wildcard name: prefer a local match unless we were asked to look for updates.
        if let Some(p) = &local {
            if !(can_download && self.check_for_updates) {
                return Ok((p.clone(), false));
            }
        }
        if !can_download {
            return local.map(|p| (p, false)).ok_or_else(|| not_found(spec, &dirs, self.auto_download));
        }
        let dir_url = spec.remote_dir(&self.naif_url).expect("checked above");
        let remote_name = match latest_remote_match(&dir_url, &spec.name) {
            Ok(n) => n,
            // Offline or listing unavailable: fall back to what we have.
            Err(e) => return local.map(|p| (p, false)).ok_or(e),
        };
        if let Some(p) = &local {
            let local_name = p.file_name().map(|n| n.to_string_lossy().to_string()).unwrap_or_default();
            if local_name >= remote_name {
                return Ok((p.clone(), false));
            }
        }
        let dest = self.download_dir.join(&remote_name);
        if dest.exists() {
            return Ok((dest, false));
        }
        let url = spec.remote_file(&self.naif_url, &remote_name).expect("checked above");
        download(&url, &dest)?;
        Ok((dest, true))
    }
}

fn not_found(spec: &KernelSpec, dirs: &[PathBuf], download_enabled: bool) -> SpiceError {
    let searched = dirs.iter().map(|d| d.display().to_string()).collect::<Vec<_>>().join(", ");
    let why = if !download_enabled {
        "automatic download is disabled"
    } else {
        "no download source is configured for it"
    };
    SpiceError::KernelNotFound {
        name: spec.name.clone(),
        detail: format!("searched [{}]; {}", searched, why),
    }
}

/// Simple glob matching supporting `*` only.
pub fn glob_match(pattern: &str, name: &str) -> bool {
    let parts: Vec<&str> = pattern.split('*').collect();
    if parts.len() == 1 {
        return pattern == name;
    }
    let mut rest = name;
    for (i, part) in parts.iter().enumerate() {
        if i == 0 {
            match rest.strip_prefix(part) {
                Some(r) => rest = r,
                None => return false,
            }
        } else if i == parts.len() - 1 {
            return rest.len() >= part.len() && rest.ends_with(part);
        } else {
            match rest.find(part) {
                Some(k) => rest = &rest[k + part.len()..],
                None => return false,
            }
        }
    }
    true
}

/// Newest local file matching `name` (which may be a pattern), searching `dirs` in order.
fn find_local(dirs: &[PathBuf], name: &str) -> Option<PathBuf> {
    if !name.contains('*') {
        return dirs.iter().map(|d| d.join(name)).find(|p| p.is_file());
    }
    let mut best: Option<(String, PathBuf)> = None;
    for d in dirs {
        let Ok(entries) = fs::read_dir(d) else { continue };
        for e in entries.flatten() {
            let fname = e.file_name().to_string_lossy().to_string();
            if glob_match(name, &fname) && e.path().is_file() && best.as_ref().map_or(true, |(b, _)| fname > *b) {
                best = Some((fname, e.path()));
            }
        }
    }
    best.map(|(_, p)| p)
}

fn http_client() -> Result<reqwest::blocking::Client> {
    reqwest::blocking::Client::builder()
        .connect_timeout(std::time::Duration::from_secs(30))
        .timeout(None)
        .user_agent(concat!("spacerocks/", env!("CARGO_PKG_VERSION")))
        .build()
        .map_err(|e| SpiceError::Download {
            url: String::new(),
            reason: e.to_string(),
        })
}

/// Name of the newest file in a remote directory listing that matches `pattern`.
fn latest_remote_match(dir_url: &str, pattern: &str) -> Result<String> {
    let url = format!("{}/", dir_url.trim_end_matches('/'));
    let fail = |reason: String| SpiceError::Download { url: url.clone(), reason };
    let resp = http_client()?.get(&url).send().map_err(|e| fail(e.to_string()))?;
    if !resp.status().is_success() {
        return Err(fail(format!("HTTP {}", resp.status())));
    }
    let body = resp.text().map_err(|e| fail(e.to_string()))?;
    let mut best: Option<String> = None;
    for piece in body.split("href=\"").skip(1) {
        let Some(target) = piece.split('"').next() else { continue };
        let fname = target.rsplit('/').next().unwrap_or(target);
        if glob_match(pattern, fname) && best.as_deref().map_or(true, |b| fname > b) {
            best = Some(fname.to_string());
        }
    }
    best.ok_or_else(|| fail(format!("no file matching '{}' in the directory listing", pattern)))
}

/// Download `url` to `dest` atomically (via a `.part` file), reporting progress on stderr.
pub fn download(url: &str, dest: &Path) -> Result<()> {
    let fail = |reason: String| SpiceError::Download {
        url: url.to_string(),
        reason,
    };
    if let Some(parent) = dest.parent() {
        fs::create_dir_all(parent).map_err(|e| SpiceError::io(&parent.display().to_string(), e))?;
    }
    let mut resp = http_client()?.get(url).send().map_err(|e| fail(e.to_string()))?;
    if !resp.status().is_success() {
        return Err(fail(format!("HTTP {}", resp.status())));
    }
    let total = resp.content_length();
    let name = dest.file_name().map(|n| n.to_string_lossy().to_string()).unwrap_or_default();
    match total {
        Some(n) => eprintln!("spacerocks: downloading {} ({:.1} MB) from {}", name, n as f64 / 1e6, url),
        None => eprintln!("spacerocks: downloading {} from {}", name, url),
    }
    let part = dest.with_file_name(format!("{}.part", name));
    let result = (|| -> std::io::Result<u64> {
        let mut out = std::io::BufWriter::new(fs::File::create(&part)?);
        let mut buf = vec![0u8; 1 << 20];
        let mut done: u64 = 0;
        let mut next_report = 0.1f64;
        loop {
            let n = resp.read(&mut buf)?;
            if n == 0 {
                break;
            }
            out.write_all(&buf[..n])?;
            done += n as u64;
            if let Some(t) = total {
                let frac = done as f64 / t as f64;
                if t > 20_000_000 && frac >= next_report {
                    eprintln!("spacerocks:   {} {:3.0}%", name, frac * 100.0);
                    next_report += 0.1;
                }
            }
        }
        out.flush()?;
        Ok(done)
    })();
    match result {
        Ok(n) => {
            if let Some(t) = total {
                if n != t {
                    let _ = fs::remove_file(&part);
                    return Err(fail(format!("incomplete download ({} of {} bytes)", n, t)));
                }
            }
            fs::rename(&part, dest).map_err(|e| SpiceError::io(&dest.display().to_string(), e))?;
            Ok(())
        }
        Err(e) => {
            let _ = fs::remove_file(&part);
            Err(fail(e.to_string()))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn glob() {
        assert!(glob_match("earth_1962_*_combined.bpc", "earth_1962_260806_2126_combined.bpc"));
        assert!(!glob_match("earth_1962_*_combined.bpc", "earth_2026_260806_2126_predict.bpc"));
        assert!(glob_match("de440s.bsp", "de440s.bsp"));
        assert!(glob_match("*.tls", "naif0012.tls"));
        assert!(!glob_match("a*b*c", "abx"));
    }

    #[test]
    fn toml_config() {
        let c = Config::from_toml(
            r#"
            auto_download = false
            download_dir = "~/k"
            default_kernels = [
              { name = "de440s.bsp", kernel_type = "spk/planets" },
              { name = "x.bsp" },
            ]
        "#,
        )
        .unwrap();
        assert!(!c.auto_download);
        assert!(!c.download_dir.starts_with("~"));
        assert_eq!(c.default_kernels.len(), 2);
        assert_eq!(
            c.default_kernels[0].remote_file(NAIF_GENERIC_KERNELS, "de440s.bsp").unwrap(),
            "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440s.bsp"
        );
        assert!(c.default_kernels[1].remote_dir(NAIF_GENERIC_KERNELS).is_none());
        // Missing fields fall back to defaults.
        let d = Config::from_toml("").unwrap();
        assert_eq!(d.default_kernels, default_kernel_list());
    }
}
