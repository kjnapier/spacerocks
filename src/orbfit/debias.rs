//! Star-catalog debiasing of optical astrometry (Eggl, Farnocchia, Chamberlin & Chesley 2020,
//! Icarus 339, 113596), as layup applies it.
//!
//! The bias table gives, for each of 26 star catalogs and each HEALPix pixel (NESTED ordering),
//! the catalog's position offset (RA·cos Dec and Dec, arcsec) and proper-motion offset (mas/yr)
//! relative to Gaia, at J2000. JPL distributes it as `bias.dat` inside `debias_hires2018.tgz`
//! (nside 256). The first load parses the text file and writes a compact binary copy next to it,
//! `bias.bin`, which later loads memory-map.

use std::fs;
use std::io::{BufRead, BufReader, Read, Write};
use std::path::{Path, PathBuf};

use memmap2::Mmap;

/// Where JPL publishes the high-resolution (nside 256) table.
pub const DEBIAS_URL: &str = "https://ssd.jpl.nasa.gov/ftp/ssd/debias/debias_hires2018.tgz";

/// The catalogs of the table, in its column order: ADES name and MPC one-letter code.
pub const BIAS_CATALOGS: [(&str, &str); 26] = [
    ("USNOA1", "a"),
    ("USNOSA1", "b"),
    ("USNOA2", "c"),
    ("USNOSA2", "d"),
    ("UCAC1", "e"),
    ("Tyc2", "g"),
    ("GSC1.1", "i"),
    ("GSC1.2", "j"),
    ("ACT", "l"),
    ("GSCACT", "m"),
    ("SDSS8", "n"),
    ("USNOB1", "o"),
    ("PPM", "p"),
    ("UCAC4", "q"),
    ("UCAC2", "r"),
    ("PPMXL", "t"),
    ("UCAC3", "u"),
    ("NOMAD", "v"),
    ("CMC14", "w"),
    ("2MASS", "L"),
    ("SDSS7", "N"),
    ("CMC15", "Q"),
    ("SSTRC4", "R"),
    ("URAT1", "S"),
    ("Gaia1", "U"),
    ("Gaia3", "W"),
];

const NCOLS: usize = 4 * BIAS_CATALOGS.len();
const HEADER_LINES: usize = 23;
const MAGIC: &[u8; 8] = b"SRBIAS1\0";

type Result<T> = std::result::Result<T, Box<dyn std::error::Error + Send + Sync>>;

/// A star-catalog bias table, memory-mapped.
pub struct BiasTable {
    nside: u64,
    map: Mmap,
}

/// Column of a catalog (by ADES name or MPC code) in the table, if it has one.
fn catalog_index(catalog: &str) -> Option<usize> {
    BIAS_CATALOGS.iter().position(|(name, code)| *name == catalog || *code == catalog)
}

impl BiasTable {
    /// Load `bias.dat` (JPL's text table), or its binary copy `bias.bin` when that is present
    /// and newer. The binary copy is written on the first load.
    pub fn load<P: AsRef<Path>>(path: P) -> Result<BiasTable> {
        let path = path.as_ref();
        let bin = path.with_extension("bin");
        let fresh = |b: &Path| -> bool {
            let (Ok(mb), Ok(mt)) = (fs::metadata(b).and_then(|m| m.modified()), fs::metadata(path).and_then(|m| m.modified())) else {
                return fs::metadata(b).is_ok() && fs::metadata(path).is_err();
            };
            mb >= mt
        };
        if !fresh(&bin) {
            write_binary(path, &bin)?;
        }
        Self::open_binary(&bin)
    }

    /// The default table, in `$SPACEROCKS_DEBIAS_DIR` or `~/.spacerocks/debias`, downloaded from
    /// JPL ([`DEBIAS_URL`]) if missing and `download` is true.
    pub fn load_default(download: bool) -> Result<BiasTable> {
        let dir = default_dir();
        let dat = dir.join("bias.dat");
        if !dat.exists() && !dir.join("bias.bin").exists() {
            if !download {
                return Err(format!("{} not found and download disabled", dat.display()).into());
            }
            fetch(&dir)?;
        }
        Self::load(dat)
    }

    fn open_binary(bin: &Path) -> Result<BiasTable> {
        let file = fs::File::open(bin)?;
        let map = unsafe { Mmap::map(&file)? };
        if map.len() < 24 || &map[..8] != MAGIC {
            return Err(format!("{} is not a spacerocks bias table", bin.display()).into());
        }
        let nside = u64::from_le_bytes(map[8..16].try_into().unwrap());
        let ncols = u64::from_le_bytes(map[16..24].try_into().unwrap()) as usize;
        let npix = 12 * nside * nside;
        if ncols != NCOLS || map.len() != 24 + npix as usize * NCOLS * 4 {
            return Err(format!("{} has the wrong size", bin.display()).into());
        }
        Ok(BiasTable { nside, map })
    }

    /// HEALPix resolution of the table.
    pub fn nside(&self) -> u64 {
        self.nside
    }

    fn value(&self, pixel: u64, column: usize) -> f64 {
        let at = 24 + (pixel as usize * NCOLS + column) * 4;
        f32::from_le_bytes(self.map[at..at + 4].try_into().unwrap()) as f64
    }

    /// Debias one detection (layup's `debias`): RA and Dec in radians, the TDB Julian date of
    /// the detection, and its star catalog (ADES name or MPC code). A catalog the table doesn't
    /// cover (Gaia DR2, EDR3, UCAC-5, ...; or none) leaves the position unchanged.
    pub fn debias(&self, ra: f64, dec: f64, jd_tdb: f64, catalog: &str) -> (f64, f64) {
        let Some(c) = catalog_index(catalog) else { return (ra, dec) };
        // layup works in degrees; so does its HEALPix lookup.
        let (ra_deg, dec_deg) = (ra.to_degrees(), dec.to_degrees());
        let pixel = ang2pix_nest_lonlat(self.nside, ra_deg, dec_deg);
        let (ra_off, dec_off) = (self.value(pixel, 4 * c), self.value(pixel, 4 * c + 1));
        let (pm_ra, pm_dec) = (self.value(pixel, 4 * c + 2), self.value(pixel, 4 * c + 3));
        // Julian years since J2000, the table's reference epoch.
        let dt = (jd_tdb - 2451545.0) / 365.25;
        let ddec = dec_off + dt * pm_dec / 1000.0;
        let dec_deb = dec_deg - ddec / 3600.0;
        let dra = (ra_off + dt * pm_ra / 1000.0) / deg2rad(dec_deg).cos();
        let ra_deb = ra_deg - dra / 3600.0;
        // Through the unit vector, as layup does (it normalizes RA into [0, 360)).
        let (a, d) = (deg2rad(ra_deb), deg2rad(dec_deb));
        let cosd = d.cos();
        let (x, y, z) = (cosd * a.cos(), cosd * a.sin(), d.sin());
        let r = (x * x + y * y + z * z).sqrt();
        let (xu, yu, zu) = (x / r, y / r, z / r);
        let ra_out = (rad2deg(yu.atan2(xu)) + 360.0).rem_euclid(360.0);
        let dec_out = rad2deg(zu.asin());
        (deg2rad(ra_out), deg2rad(dec_out))
    }
}

// numpy's deg2rad / rad2deg: one multiplication by a rounded constant.
fn deg2rad(x: f64) -> f64 {
    x * (std::f64::consts::PI / 180.0)
}
fn rad2deg(x: f64) -> f64 {
    x * (180.0 / std::f64::consts::PI)
}

fn default_dir() -> PathBuf {
    if let Some(d) = std::env::var_os("SPACEROCKS_DEBIAS_DIR").filter(|d| !d.is_empty()) {
        return PathBuf::from(d);
    }
    dirs::home_dir().unwrap_or_else(|| PathBuf::from(".")).join(".spacerocks").join("debias")
}

/// Download JPL's archive and extract `bias.dat` into `dir`.
fn fetch(dir: &Path) -> Result<()> {
    fs::create_dir_all(dir)?;
    let tgz = dir.join("debias_hires2018.tgz");
    crate::spice::config::download(DEBIAS_URL, &tgz)?;
    let mut archive = tar::Archive::new(flate2::read::GzDecoder::new(fs::File::open(&tgz)?));
    for entry in archive.entries()? {
        let mut entry = entry?;
        if entry.path()?.file_name().map(|n| n == "bias.dat").unwrap_or(false) {
            let mut out = fs::File::create(dir.join("bias.dat"))?;
            std::io::copy(&mut entry, &mut out)?;
            return Ok(());
        }
    }
    Err(format!("no bias.dat in {}", tgz.display()).into())
}

/// Parse JPL's text table (23 header lines, then one row of 104 numbers per pixel) into the
/// binary layout: magic, nside, column count, then f32 values row by row.
fn write_binary(dat: &Path, bin: &Path) -> Result<()> {
    let mut reader = BufReader::with_capacity(1 << 20, fs::File::open(dat)?);
    let mut line = String::new();
    for _ in 0..HEADER_LINES {
        line.clear();
        reader.read_line(&mut line)?;
    }
    let mut values: Vec<f32> = Vec::new();
    let mut rows = 0usize;
    let mut rest = String::new();
    reader.read_to_string(&mut rest)?;
    for (k, l) in rest.lines().enumerate() {
        if l.trim().is_empty() {
            continue;
        }
        let before = values.len();
        for tok in l.split_whitespace() {
            values.push(tok.parse::<f64>().map_err(|e| format!("{}: line {}: {}", dat.display(), HEADER_LINES + k + 1, e))? as f32);
        }
        if values.len() - before != NCOLS {
            return Err(format!("{}: line {} has {} values, expected {}", dat.display(), HEADER_LINES + k + 1, values.len() - before, NCOLS).into());
        }
        rows += 1;
    }
    let nside = ((rows / 12) as f64).sqrt().round() as u64;
    if 12 * nside * nside != rows as u64 || !nside.is_power_of_two() {
        return Err(format!("{}: {} rows is not a HEALPix map", dat.display(), rows).into());
    }
    let tmp = bin.with_extension("bin.part");
    {
        let mut out = std::io::BufWriter::new(fs::File::create(&tmp)?);
        out.write_all(MAGIC)?;
        out.write_all(&nside.to_le_bytes())?;
        out.write_all(&(NCOLS as u64).to_le_bytes())?;
        for v in &values {
            out.write_all(&v.to_le_bytes())?;
        }
        out.flush()?;
    }
    fs::rename(&tmp, bin)?;
    Ok(())
}

/// HEALPix pixel (NESTED) of a longitude/latitude in degrees, as healpy's
/// `ang2pix(nside, lon, lat, nest=True, lonlat=True)`: theta = pi/2 - lat, then HEALPix's
/// `loc2pix` (with sin(theta) near the poles).
pub fn ang2pix_nest_lonlat(nside: u64, lon: f64, lat: f64) -> u64 {
    let theta = std::f64::consts::FRAC_PI_2 - deg2rad(lat);
    let phi = deg2rad(lon);
    let near_pole = theta < 0.01 || theta > 3.14159 - 0.01;
    let sth = if near_pole { theta.sin() } else { 0.0 };
    loc2pix_nest(nside, theta.cos(), phi, sth, near_pole)
}

fn fmodulo(v1: f64, v2: f64) -> f64 {
    if v1 >= 0.0 {
        if v1 < v2 {
            v1
        } else {
            v1 % v2
        }
    } else {
        let tmp = v1 % v2 + v2;
        if tmp == v2 {
            0.0
        } else {
            tmp
        }
    }
}

fn spread_bits(v: u64) -> u64 {
    let mut out = 0;
    for b in 0..32 {
        out |= ((v >> b) & 1) << (2 * b);
    }
    out
}

fn loc2pix_nest(nside: u64, z: f64, phi: f64, sth: f64, have_sth: bool) -> u64 {
    const INV_HALFPI: f64 = 0.6366197723675813430755350534900574;
    let order = nside.trailing_zeros();
    let ns = nside as i64;
    let za = z.abs();
    let tt = fmodulo(phi * INV_HALFPI, 4.0);
    let (ix, iy, face) = if za <= 2.0 / 3.0 {
        let temp1 = nside as f64 * (0.5 + tt);
        let temp2 = nside as f64 * (z * 0.75);
        let jp = (temp1 - temp2) as i64;
        let jm = (temp1 + temp2) as i64;
        let (ifp, ifm) = (jp >> order, jm >> order);
        let face = if ifp == ifm { ifp | 4 } else if ifp < ifm { ifp } else { ifm + 8 };
        (jm & (ns - 1), ns - (jp & (ns - 1)) - 1, face)
    } else {
        let ntt = (tt as i64).min(3);
        let tp = tt - ntt as f64;
        let tmp = if za < 0.99 || !have_sth { nside as f64 * (3.0 * (1.0 - za)).sqrt() } else { nside as f64 * sth / ((1.0 + za) / 3.0).sqrt() };
        let jp = ((tp * tmp) as i64).min(ns - 1);
        let jm = (((1.0 - tp) * tmp) as i64).min(ns - 1);
        if z >= 0.0 {
            (ns - jm - 1, ns - jp - 1, ntt)
        } else {
            (jp, jm, ntt + 8)
        }
    };
    ((face as u64) << (2 * order)) + spread_bits(ix as u64) + (spread_bits(iy as u64) << 1)
}
