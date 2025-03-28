use std::collections::HashMap;
use std::fs::File;
use std::io::{self, Read, Seek, SeekFrom};
use std::mem;
use std::path::Path;
use std::slice;

use memmap2::Mmap;

use crate::spice::error::SpiceError;

const RECORD_LENGTH: usize = 1024;

/// Convert ephemeris seconds past J2000.0 to Julian Day.
fn jul(eph: f64) -> f64 {
    2451545.0 + eph / 86400.0
}

/// The DAF file record for a binary PCK file (1024 bytes).
#[repr(C)]
#[derive(Copy, Clone, Debug)]
struct FileRecord {
    locidw: [u8; 8],   // Identification word, e.g. "DAF/PCK"
    nd: i32,           // Number of double precision components per summary (should be 2)
    ni: i32,           // Number of integer components per summary (should be 5)
    locifn: [u8; 60],  // Internal file name/description
    fward: i32,        // Record number of the first summary record
    bward: i32,        // Record number of the final summary record
    free: i32,         // First free address in the file
    locfmt: [u8; 8],   // Numeric binary format string ("LTL-IEEE" or "BIG-IEEE")
    prenul: [u8; 603], // Padding null block
    ftpstr: [u8; 28],  // FTP validation string
    pstnul: [u8; 297], // Padding null block
}

/// The summary for a PCK segment. For an angles‑only (Type 2) segment the summary
/// contains two doubles (begin and end epochs) followed by five integers:
///   frame_class, inertial_code, representation, init_addr, final_addr,
/// plus one padding integer to fill out 40 bytes.
#[repr(C)]
#[derive(Debug, Clone, Copy)]
struct PckSummary {
    beg: f64,           // Begin epoch (ephemeris seconds past J2000.0)
    end: f64,           // End epoch
    frame_class: i32,   // Frame class ID
    inertial_code: i32, // NAIF code for the inertial reference frame
    representation: i32,// Representation type (should be 2 for angles‑only)
    init_addr: i32,     // Initial address of the segment data
    final_addr: i32,    // Final address of the segment data
    _pad: i32,          // Padding (unused)
}

/// A summary record block contains a header (three doubles) followed by up to 25 summaries.
#[repr(C)]
#[derive(Copy, Clone, Debug)]
struct PckSummaryRecord {
    next: f64,                   // Record number of the next summary record (0 if final)
    prev: f64,                   // Record number of the previous summary record (0 if initial)
    nsum: f64,                   // Number of summaries in this record
    s: [PckSummary; 25],         // Up to 25 summaries in one record
}

/// A PCK segment (extracted from a summary).
#[derive(Debug)]
pub struct PckSegment {
    pub beg: f64,           // Beginning epoch (in Julian Day)
    pub end: f64,           // Ending epoch (in Julian Day)
    pub frame_class: i32,
    pub inertial_code: i32,
    pub representation: i32,
    pub init_addr: i32,     // Starting address of the data block (1-indexed)
    pub final_addr: i32,    // Ending address of the data block (1-indexed)
}

/// The top-level binary PCK (BPC) structure.
#[derive(Debug)]
pub struct Bpc {
    pub map: Mmap,                   // Memory-mapped file contents
    pub len: usize,                  // Total file length in bytes
    pub segments: Vec<PckSegment>,   // Parsed segments
    // A hash map mapping an inertial frame code to segment indices.
    pub segment_map: HashMap<i32, Vec<usize>>,
}

impl Bpc {
    /// Open and initialize a binary PCK file.
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self, SpiceError> {
        let mut file = File::open(&path)?;

        // Read the file record (first 1024 bytes)
        let mut buf = [0u8; RECORD_LENGTH];
        file.read_exact(&mut buf)?;
        let file_record: FileRecord = unsafe { ptr_from_bytes(&buf)? };

        // println!("{:?}", file_record);

        if &file_record.locidw[..7] != b"DAF/PCK" {
            return Err(SpiceError::ParseError(
                "Error parsing DAF/PCK file. Incorrect header.".to_string(),
            ));
        }
        if file_record.nd != 2 || file_record.ni != 5 {
            return Err(SpiceError::ParseError(
                "Unexpected ND or NI in PCK file.".to_string(),
            ));
        }
        if mem::size_of::<PckSummary>() != 40 {
            return Err(SpiceError::ParseError(
                "Wrong size of PCK summary record.".to_string(),
            ));
        }

        // Check that the summary record size matches our expectation.
        let nc = 8 * (file_record.nd as usize + ((file_record.ni as usize + 1) / 2));
        if nc != mem::size_of::<PckSummary>() {
            return Err(SpiceError::ParseError(
                "Error parsing DAF/SPK file. Wrong size of summary record.".to_string(),
            ));
        }

        // Seek to the first summary record.
        let first_summary_record = (file_record.fward as u64 - 1) * RECORD_LENGTH as u64;
        file.seek(SeekFrom::Start(first_summary_record))?;
        file.read_exact(&mut buf)?;

        if buf[8] != 0 {
            return Err(SpiceError::ParseError(
                "Error parsing DAF/PCK file. Cannot find summary block.".to_string(),
            ));
        }

        let mut segments: Vec<PckSegment> = Vec::new();
        loop {
            let summary_record: PckSummaryRecord = unsafe { ptr_from_bytes(&buf)? };
            // println!("{:?}", summary_record);
            let nsum = summary_record.nsum as usize;

            for i in 0..nsum {
                let sum = summary_record.s[i];
                let segment = PckSegment {
                    beg: jul(sum.beg),
                    end: jul(sum.end),
                    frame_class: sum.frame_class,
                    inertial_code: sum.inertial_code,
                    representation: sum.representation,
                    init_addr: sum.init_addr,
                    final_addr: sum.final_addr,
                };
                segments.push(segment);
            }
            let next_record = summary_record.next as i64 - 1;
            if next_record < 0 {
                break;
            } else {
                file.seek(SeekFrom::Start(next_record as u64 * RECORD_LENGTH as u64))?;
                file.read_exact(&mut buf)?;
            }
        }

        let len = file.metadata()?.len() as usize;
        let map = unsafe { Mmap::map(&file)? };

        let mut segment_map: HashMap<i32, Vec<usize>> = HashMap::new();
        for (i, segment) in segments.iter().enumerate() {
            segment_map
                .entry(segment.inertial_code)
                .or_default()
                .push(i);
        }

        Ok(Bpc {
            map,
            len,
            segments,
            segment_map,
        })
    }

    /// Compute the rotation matrix from a Type 2 (angles-only Chebyshev)
    /// PCK segment at time `t` (in Julian Day).
    ///
    /// This function assumes the segment data consists of a single Chebyshev record
    /// with the following layout (all values are double precision):
    ///   [ t_mid, t_half, coeffs_phi..., coeffs_theta..., coeffs_psi... ]
    ///
    ///   - t_mid: the midpoint time of the record (ephemeris seconds past J2000.0)
    ///   - t_half: half the record duration (in seconds)
    ///   - The number of coefficients per angle, p, is determined by:
    ///         p = (n_doubles - 2) / 3
    ///   - The Chebyshev polynomials are evaluated at the normalized time:
    ///         τ = (t - t_mid_jd) / (t_half / 86400)
    ///     where t_mid_jd is t_mid converted to Julian Day.
    ///
    /// The Euler angles (φ, θ, ψ) are then used in a 3-1-3 rotation:
    ///   R = R_z(ψ) · R_x(θ) · R_z(φ)
    pub fn rotation_matrix_at(&self, segment: &PckSegment, t: f64) -> Result<[[f64; 3]; 3], SpiceError> {
        // Check that time t is within the segment interval.
        if t < segment.beg || t > segment.end {
            return Err(SpiceError::ParseError("Time out of segment bounds".to_string()));
        }

        // We assume that the segment data consists of a single Chebyshev record, Type 2
        // https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/req/pck.html#Type%202:%20Chebyshev%20Angles%20only
        // +---------------+
        // | Record 1      |
        // +---------------+
        // | Record 2      |
        // +---------------+
        //   .
        //   .
        //   .
        // +---------------+
        // | Record N      |
        // +---------------+
        // | INIT          |
        // +---------------+
        // | INTLEN        |
        // +---------------+
        // | RSIZE         |
        // +---------------+
        // | N             |
        // +---------------+
        //
        // 1. INIT is the initial epoch of the first record, given in ephemeris seconds past 2000 Jan 01 12:00:00, also known as J2000.
        // 2. INTLEN is the length of the interval covered by each record, in seconds.
        // 3. RSIZE is the total size of (number of array elements in) each record.
        // 4. N is the number of records contained in the segment.


        // println!("t: {}", t);
        // println!("segment: {:?}", segment);

        // Access the file's data as f64 slice.
        let f64_slice = unsafe {
            slice::from_raw_parts(
                self.map.as_ptr() as *const f64,
                self.len / mem::size_of::<f64>(),
            )
        };

        // Convert addresses (1-indexed) to slice indices.
        let start_idx = (segment.init_addr - 1) as usize;
        let end_idx = (segment.final_addr) as usize; // inclusive
        if end_idx > f64_slice.len() || start_idx >= end_idx {
            return Err(SpiceError::ParseError("Segment addresses out of range".to_string()));
        }
        let record_data = &f64_slice[start_idx..end_idx];

        let n_records = record_data[record_data.len() - 1] as usize;
        let rsize = record_data[record_data.len() - 2] as usize;
        let intlen = record_data[record_data.len() - 3];
        let init = record_data[record_data.len() - 4];

        let total_records_length = record_data.len() - 4;
        let init_jd = jul(init);
        let dt = (t - init_jd) * 86400.0;
        let int_record_idx = (dt / intlen) as usize;
        let start = int_record_idx * rsize;
        let end = start + rsize;
        let record_data = &record_data[start..end];

        if record_data.len() < 2 {
            return Err(SpiceError::ParseError("Not enough data in record".to_string()));
        }

        let t_mid = record_data[0];
        let t_mid = jul(t_mid);
        let t_half = record_data[1];
        let tau = (t - t_mid) / (t_half / 86400.0);

        // Determine the number of coefficients per angle.
        let n_doubles = record_data.len();  
        if (n_doubles - 2) % 3 != 0 {
            return Err(SpiceError::ParseError("Unexpected number of coefficients".to_string()));
        }

        let p = (n_doubles - 2) / 3;        
        let coeffs_phi   = &record_data[2 .. 2 + p];
        let coeffs_theta = &record_data[2 + p .. 2 + 2 * p];
        let coeffs_psi   = &record_data[2 + 2 * p .. 2 + 3 * p];

        let phi   = eval_chebyshev(tau, coeffs_phi);
        let theta = eval_chebyshev(tau, coeffs_theta);
        let psi   = eval_chebyshev(tau, coeffs_psi);

        // Compute rotation matrices for each Euler angle.
        let r_z_phi = rot_z(phi);
        let r_x_theta = rot_x(theta);
        let r_z_psi = rot_z(psi);


        // For a 3-1-3 rotation: R = R_z(ψ) · R_x(θ) · R_z(φ)
        let r_temp = mat_mult(&r_x_theta, &r_z_psi);
        let r = mat_mult(&r_z_phi, &r_temp);

        Ok(r)
    }
}

impl Bpc {
    // Existing methods (e.g. open, rotation_matrix_at) remain unchanged

    /// Compute the rotation matrix for a given epoch (in Julian Day).
    /// This function automatically finds the appropriate PCK segment covering the epoch.
    pub fn rotation_matrix_at_epoch(&self, t: f64) -> Result<[[f64; 3]; 3], SpiceError> {
        // Search for a segment that covers time t.
        let segment = self.segments
            .iter()
            .find(|s| t >= s.beg && t <= s.end)
            .ok_or_else(|| {
                SpiceError::ParseError("No segment covers the requested epoch".to_string())
            })?;
        
        // Compute the rotation matrix for the found segment.
        self.rotation_matrix_at(segment, t)
    }
}


/// Evaluate a Chebyshev polynomial at τ using Clenshaw's recurrence.
/// coeffs[0..p] are the Chebyshev coefficients.
fn eval_chebyshev(tau: f64, coeffs: &[f64]) -> f64 {
    let n = coeffs.len();
    let mut b_kplus1 = 0.0;
    let mut b_kplus2 = 0.0;
    for &a in coeffs.iter().rev() {
        let b_k = 2.0 * tau * b_kplus1 - b_kplus2 + a;
        b_kplus2 = b_kplus1;
        b_kplus1 = b_k;
    }
    b_kplus1 - tau * b_kplus2
}

/// Returns a 3x3 rotation matrix for a rotation about the Z-axis by `angle` radians.
fn rot_z(angle: f64) -> [[f64; 3]; 3] {
    let (s, c) = angle.sin_cos();
    [
        [ c, -s, 0.0 ],
        [ s,  c, 0.0 ],
        [0.0, 0.0, 1.0],
    ]
}

/// Returns a 3x3 rotation matrix for a rotation about the X-axis by `angle` radians.
fn rot_x(angle: f64) -> [[f64; 3]; 3] {
    let (s, c) = angle.sin_cos();
    [
        [1.0, 0.0,  0.0],
        [0.0,  c, -s],
        [0.0,  s,  c],
    ]
}

/// Multiply two 3x3 matrices: result = a · b.
fn mat_mult(a: &[[f64; 3]; 3], b: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut res = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                res[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    res
}

/// Unsafe helper: Convert a byte buffer into a struct of type `T`.
unsafe fn ptr_from_bytes<T: Copy>(bytes: &[u8]) -> Result<T, SpiceError> {
    if bytes.len() < mem::size_of::<T>() {
        return Err(SpiceError::ParseError("Buffer too small".to_string()));
    }
    let ptr = bytes.as_ptr() as *const T;
    Ok(ptr.read_unaligned())
}
