use std::collections::HashMap;
use std::fs::File;
use std::io::{self, Read, Seek, SeekFrom};
use std::mem;
use std::path::Path;
use std::slice;


use memmap2::Mmap;

use crate::spice::error::SpiceError;

const RECORD_LENGTH: usize = 1024;

/// Convert SPK epoch (seconds since J2000.0) to Julian day.
fn jul(eph: f64) -> f64 {
    2451545.0 + eph / 86400.0
}

/// Representation of one summary record from the SPK file.
#[repr(C)]
#[derive(Debug, Clone, Copy)]
struct Sum {
    beg: f64, // begin epoch, seconds since J2000.0
    end: f64, // ending epoch
    tar: i32, // target code
    cen: i32, // centre code
    r#ref: i32, // reference frame (e.g., 1 = J2000.0)
    ver: i32, // type of ephemeris (e.g., 2 = Chebyshev)
    one: i32, // initial array address
    two: i32, // final array address
}

/// A summary record block. (Note: The original file layout defines 3 doubles then 25 summary entries.)
#[repr(C)]
#[derive(Copy, Clone)]
struct SummaryRecord {
    next: f64,       // record number of the next summary record (0 if final)
    prev: f64,       // record number of the previous summary record (0 if initial)
    nsum: f64,       // number of summaries in this record
    s: [Sum; 25],    // up to 25 summaries in one record
}

/// The file record (header) at the start of the file.
#[repr(C)]
#[derive(Copy, Clone)]
struct FileRecord {
    locidw: [u8; 8], // identification word, e.g. "DAF/SPK"
    nd: i32,         // number of double precision components per summary
    ni: i32,         // number of integer components per summary
    locifn: [u8; 60],// internal name/description of the file
    fward: i32,      // record number of the first summary record
    bward: i32,      // record number of the final summary record
}

/// A target with its associated summary data.
#[derive(Debug)]
pub struct SpkTarget {
    pub code: i32,
    pub cen: i32,
    pub beg: f64,
    pub res: f64,
    pub one: Vec<i32>,
    pub two: Vec<i32>,
    pub ind: usize,
    pub end: f64,
    pub mass: f64, // mass constant; if not set then 0.0 is used
}

/// The top-level SPK structure.
pub struct Spk {
    pub map: Mmap, // memory-mapped file contents
    pub len: usize,
    pub targets: Vec<SpkTarget>,
    // hash map of target codes to indices
    pub target_map: HashMap<i32, usize>,
}


impl Spk {
    /// Open and initialize an SPK file given by `path`.
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self, SpiceError> {
        let mut file = File::open(&path)?;

        // Read the file record (first 1024 bytes)
        let mut buf = [0u8; RECORD_LENGTH];
        file.read_exact(&mut buf)?;

        // Interpret the first record as FileRecord.
        let file_record: FileRecord = unsafe { ptr_from_bytes(&buf)? };

        // Check for a valid DAF/SPK header.
        if &file_record.locidw[..7] != b"DAF/SPK" {
            return Err(SpiceError::ParseError(
                "Error parsing DAF/SPK file. Incorrect header.".to_string(),
            ));
        }

        // Check that the summary record size matches our expectation.
        let nc = 8 * (file_record.nd as usize + ((file_record.ni as usize + 1) / 2));
        if nc != mem::size_of::<Sum>() {
            return Err(SpiceError::ParseError(
                "Error parsing DAF/SPK file. Wrong size of summary record.".to_string(),
            ));
        }

        // Seek to the first summary record.
        let first_summary_record = (file_record.fward as u64 - 1) * RECORD_LENGTH as u64;
        file.seek(SeekFrom::Start(first_summary_record))?;
        file.read_exact(&mut buf)?;

        // Validate that we have a summary block.
        if buf[8] != 0 {
            return Err(SpiceError::ParseError(
                "Error parsing DAF/SPK file. Cannot find summary block.".to_string(),
            ));
        }

        let mut targets: Vec<SpkTarget> = Vec::new();

        // Loop over summary records.
        loop {
            let summary_record: SummaryRecord = unsafe { ptr_from_bytes(&buf)? };
            let nsum = summary_record.nsum as usize;

            for i in 0..nsum {
                let sum = summary_record.s[i];
                // Check if we need to create a new target.
                if targets.is_empty() || targets.last().unwrap().code != sum.tar {
                    let new_target = SpkTarget {
                        code: sum.tar,
                        cen: sum.cen,
                        beg: jul(sum.beg),
                        res: jul(sum.end) - jul(sum.beg),
                        one: Vec::new(),
                        two: Vec::new(),
                        ind: 0,
                        end: jul(sum.end),
                        mass: 0.0,
                    };
                    targets.push(new_target);
                }
                // Append indices for the current target.
                if let Some(target) = targets.last_mut() {
                    target.one.push(sum.one);
                    target.two.push(sum.two);
                    target.end = jul(sum.end);
                    target.ind += 1;
                }
            }

            // Determine the next summary record.
            let next_record = summary_record.next as i64 - 1;
            if next_record < 0 {
                break;
            } else {
                file.seek(SeekFrom::Start(next_record as u64 * RECORD_LENGTH as u64))?;
                file.read_exact(&mut buf)?;
            }
        }

        // Get file size.
        let len = file.metadata()?.len() as usize;

        // Memory-map the file.
        let map = unsafe { Mmap::map(&file)? };

        // Create a hash map of target codes to indices.
        let target_map = targets
            .iter()
            .enumerate()
            .map(|(i, target)| (target.code, i))
            .collect();

        // Optionally, you could advise the OS for random access (using libc::madvise).
        // For simplicity, this version omits that step.

        Ok(Spk { map, len, targets, target_map })
    }

    /// Calculate position data given the Julian date `jde`, a time offset `rel`,
    /// and target index `m`. On success returns a tuple (GM, x, y, z).
    pub fn calc(&self, jde: f64, rel: f64, m: usize) -> Result<(f64, f64, f64, f64, f64, f64, i32), SpiceError> {
        if m >= self.targets.len() {
            return Err(SpiceError::Nast);
        }
        let target = &self.targets[m];

        if jde + rel < target.beg || jde + rel > target.end {
            return Err(SpiceError::Coverage);
        }

        let mut pos_u = [0.0f64; 3];
        let mut vel_u = [0.0f64; 3];
        // Note: The v component is computed in the C version but not used.

        // Find the appropriate summary block index.
        let n = ((jde + rel - target.beg) / target.res) as usize;

        // Interpret the mapped file as a slice of f64 values.
        let f64_slice = unsafe {
            slice::from_raw_parts(
                self.map.as_ptr() as *const f64,
                self.len / mem::size_of::<f64>(),
            )
        };

        // Calculate pointer to data using target.two.
        let two_entry = target.two.get(n).ok_or_else(|| {
            SpiceError::ParseError("Index out of range in target.two".to_string())
        })?;
        let index_val = (*two_entry - 1) as usize;
        if index_val >= f64_slice.len() {
            return Err(SpiceError::ParseError(
                "Index computed from target.two out of range".to_string(),
            ));
        }
        // In the C code, R is stored immediately before the data block.
        if index_val < 1 {
            return Err(SpiceError::ParseError("Invalid index for R".to_string()));
        }
        let r = f64_slice[index_val - 1] as i32;
        let p_count = (r - 2) / 3;
        if p_count < 0 || p_count as usize > 32 {
            return Err(SpiceError::ParseError(
                "Invalid number of coefficients".to_string(),
            ));
        }
        let p_count = p_count as usize;

        // Determine record timing information from data preceding the record.
        if index_val < 3 {
            return Err(SpiceError::ParseError(
                "Not enough data for record header".to_string(),
            ));
        }
        let record_start = f64_slice[index_val - 3];
        let record_duration = f64_slice[index_val - 2] / 86400.0;
        let b = (((jde - jul(record_start)) + rel) / record_duration) as i32;


        // Locate the start of the precise record.
        let one_entry = target.one.get(n).ok_or_else(|| {
            SpiceError::ParseError("Index out of range in target.one".to_string())
        })?;
        let base_index = (*one_entry - 1) as usize + (b as usize * r as usize);
        if base_index >= f64_slice.len() {
            return Err(SpiceError::ParseError(
                "Computed base index out of range".to_string(),
            ));
        }
        let record_ptr = &f64_slice[base_index..];

        // Scale to interpolation units.
        if record_ptr.len() < 2 {
            return Err(SpiceError::ParseError(
                "Not enough data in record for interpolation".to_string(),
            ));
        }
        let t0 = record_ptr[0];
        let t1 = record_ptr[1] / 86400.0;
        let z = ((jde - jul(t0)) + rel) / t1;

        // Set up Chebyshev polynomials.
        let mut T = [0.0f64; 32];
        let mut S = [0.0f64; 32];
        T[0] = 1.0;
        T[1] = z;
        S[0] = 0.0;
        S[1] = 1.0;
        for p in 2..p_count {
            T[p] = 2.0 * z * T[p - 1] - T[p - 2];
            S[p] = 2.0 * z * S[p - 1] + 2.0 * T[p - 1] - S[p - 2];
        }

        // Sum interpolation coefficients for each coordinate.
        for n_coord in 0..3 {
            let coeff_start = 2 + n_coord * p_count;
            let coeff_end = coeff_start + p_count;
            if coeff_end > record_ptr.len() {
                return Err(SpiceError::ParseError(
                    "Not enough coefficients in record".to_string(),
                ));
            }
            for p in 0..p_count {
                pos_u[n_coord] += record_ptr[coeff_start + p] * T[p];
                vel_u[n_coord] += record_ptr[coeff_start + p] * S[p];
            }
            pos_u[n_coord] /= 149597870.7;
            vel_u[n_coord] /= 149597870.7;
            vel_u[n_coord] /= t1;
            
        }

        Ok((pos_u[0], pos_u[1], pos_u[2], vel_u[0], vel_u[1], vel_u[2], target.cen))
    }

    pub fn state_at(&self, epoch: f64, target: usize) -> Result<(f64, f64, f64, f64, f64, f64, i32), SpiceError> {
        let rel = epoch - 2451545.0;
        self.calc(2451545.0, rel, target)
    }
}


/// Unsafe helper: Convert a byte buffer into a struct of type `T`.
unsafe fn ptr_from_bytes<T: Copy>(bytes: &[u8]) -> Result<T, SpiceError> {
    if bytes.len() < mem::size_of::<T>() {
        return Err(SpiceError::ParseError("Buffer too small".to_string()));
    }
    let ptr = bytes.as_ptr() as *const T;
    Ok(ptr.read_unaligned())
}