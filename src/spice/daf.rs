//! Generic reader for NAIF Double precision Array Files (DAF).
//!
//! SPK and binary PCK kernels are DAFs. A DAF is a sequence of 1024-byte records: a file
//! record, optional comment records, then a doubly-linked list of summary records (each
//! followed by a name record) that describe arrays of double precision numbers stored in the
//! remaining records.
//!
//! Files in the host's native byte order are memory mapped and read in place with zero copies.
//! Files in the other byte order are read once and byte-swapped into memory.

use std::fs::File;
use std::path::Path;

use memmap2::Mmap;

use super::error::{Result, SpiceError};

const RECORD_BYTES: usize = 1024;
const RECORD_WORDS: usize = 128;

/// Byte order of the numeric data in a DAF.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ByteOrder {
    Little,
    Big,
}

impl ByteOrder {
    fn native() -> ByteOrder {
        if cfg!(target_endian = "little") {
            ByteOrder::Little
        } else {
            ByteOrder::Big
        }
    }
}

/// One array summary: `nd` doubles followed by `ni` integers, plus its name.
#[derive(Debug, Clone)]
pub struct DafSummary {
    pub doubles: Vec<f64>,
    pub ints: Vec<i32>,
    pub name: String,
}

enum Storage {
    Mapped(Mmap),
    Owned(Vec<f64>),
}

/// An open DAF, with all array summaries read at open time.
pub struct Daf {
    pub path: String,
    /// The file's identification word, e.g. "DAF/SPK" or "DAF/PCK" (trimmed).
    pub id_word: String,
    pub internal_name: String,
    pub nd: usize,
    pub ni: usize,
    pub byte_order: ByteOrder,
    pub summaries: Vec<DafSummary>,
    storage: Storage,
    n_words: usize,
}

impl std::fmt::Debug for Daf {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Daf")
            .field("path", &self.path)
            .field("id_word", &self.id_word)
            .field("nd", &self.nd)
            .field("ni", &self.ni)
            .field("byte_order", &self.byte_order)
            .field("arrays", &self.summaries.len())
            .finish()
    }
}

fn read_i32(bytes: &[u8], order: ByteOrder) -> i32 {
    let b: [u8; 4] = bytes[..4].try_into().unwrap();
    match order {
        ByteOrder::Little => i32::from_le_bytes(b),
        ByteOrder::Big => i32::from_be_bytes(b),
    }
}

fn read_f64(bytes: &[u8], order: ByteOrder) -> f64 {
    let b: [u8; 8] = bytes[..8].try_into().unwrap();
    match order {
        ByteOrder::Little => f64::from_le_bytes(b),
        ByteOrder::Big => f64::from_be_bytes(b),
    }
}

fn ascii(bytes: &[u8]) -> String {
    bytes
        .iter()
        .map(|&b| if (32..127).contains(&b) { b as char } else { ' ' })
        .collect::<String>()
        .trim()
        .to_string()
}

impl Daf {
    /// Open and index a DAF file.
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Daf> {
        let path_str = path.as_ref().display().to_string();
        let file = File::open(path.as_ref()).map_err(|e| SpiceError::io(&path_str, e))?;
        // SAFETY: the file is opened read-only. As with any memory map, modifying the file on
        // disk while it is mapped is undefined behaviour; kernels are treated as immutable.
        let mmap = unsafe { Mmap::map(&file) }.map_err(|e| SpiceError::io(&path_str, e))?;
        let bytes: &[u8] = &mmap;

        if bytes.len() < RECORD_BYTES {
            return Err(SpiceError::invalid(&path_str, "file is shorter than one DAF record"));
        }

        let id_word = ascii(&bytes[0..8]);
        if !(id_word.starts_with("DAF/") || id_word == "NAIF/DAF") {
            return Err(SpiceError::invalid(
                &path_str,
                format!("not a binary DAF (id word '{}'); transfer-format files must be converted with tobin first", id_word),
            ));
        }

        // Determine byte order: from LOCFMT if present, otherwise by plausibility of ND/NI.
        let locfmt = ascii(&bytes[88..96]);
        let order = match locfmt.as_str() {
            "LTL-IEEE" => ByteOrder::Little,
            "BIG-IEEE" => ByteOrder::Big,
            _ => {
                let plausible = |o: ByteOrder| {
                    let nd = read_i32(&bytes[8..12], o);
                    let ni = read_i32(&bytes[12..16], o);
                    (0..=124).contains(&nd) && (2..=250).contains(&ni)
                };
                if plausible(ByteOrder::native()) {
                    ByteOrder::native()
                } else if plausible(ByteOrder::Big) {
                    ByteOrder::Big
                } else {
                    ByteOrder::Little
                }
            }
        };

        let nd = read_i32(&bytes[8..12], order);
        let ni = read_i32(&bytes[12..16], order);
        if !(0..=124).contains(&nd) || !(2..=250).contains(&ni) {
            return Err(SpiceError::invalid(&path_str, format!("implausible ND={} NI={}", nd, ni)));
        }
        let (nd, ni) = (nd as usize, ni as usize);
        let internal_name = ascii(&bytes[16..76]);
        let fward = read_i32(&bytes[76..80], order);

        // Integer components are packed two per double.
        let ss = nd + (ni + 1) / 2; // summary size in doubles
        let nc = 8 * ss; // name size in characters

        let mut summaries = Vec::new();
        let mut rec = fward;
        let mut visited = 0usize;
        let n_records = bytes.len() / RECORD_BYTES;
        while rec > 0 {
            visited += 1;
            if visited > n_records {
                return Err(SpiceError::invalid(&path_str, "cycle in summary record list"));
            }
            let rec_idx = rec as usize - 1;
            if rec_idx >= n_records {
                return Err(SpiceError::invalid(&path_str, "summary record out of range"));
            }
            let sbase = rec_idx * RECORD_BYTES;
            let srec = &bytes[sbase..sbase + RECORD_BYTES];
            let next_raw = read_f64(&srec[0..8], order);
            if !(0.0..=(n_records as f64)).contains(&next_raw) {
                return Err(SpiceError::invalid(&path_str, "invalid next-summary-record pointer"));
            }
            let next = next_raw as i32;
            let nsum_raw = read_f64(&srec[16..24], order);
            if !(0.0..=125.0).contains(&nsum_raw) {
                return Err(SpiceError::invalid(&path_str, "invalid summary count in summary record"));
            }
            let nsum = nsum_raw as usize;
            if 3 + nsum * ss > RECORD_WORDS {
                return Err(SpiceError::invalid(&path_str, "too many summaries in summary record"));
            }
            let name_rec = if sbase + 2 * RECORD_BYTES <= bytes.len() {
                Some(&bytes[sbase + RECORD_BYTES..sbase + 2 * RECORD_BYTES])
            } else {
                None
            };
            for i in 0..nsum {
                let off = 24 + i * ss * 8;
                let mut doubles = Vec::with_capacity(nd);
                for k in 0..nd {
                    doubles.push(read_f64(&srec[off + 8 * k..], order));
                }
                let ioff = off + nd * 8;
                let mut ints = Vec::with_capacity(ni);
                for k in 0..ni {
                    ints.push(read_i32(&srec[ioff + 4 * k..], order));
                }
                let name = name_rec
                    .map(|nr| {
                        let s = i * nc;
                        if s + nc <= nr.len() {
                            ascii(&nr[s..s + nc])
                        } else {
                            String::new()
                        }
                    })
                    .unwrap_or_default();
                summaries.push(DafSummary { doubles, ints, name });
            }
            rec = next;
        }

        let n_words = bytes.len() / 8;
        let aligned = (bytes.as_ptr() as usize).trailing_zeros() >= 3;
        let storage = if order == ByteOrder::native() && aligned {
            Storage::Mapped(mmap)
        } else {
            let mut words = Vec::with_capacity(n_words);
            for k in 0..n_words {
                words.push(read_f64(&bytes[8 * k..], order));
            }
            Storage::Owned(words)
        };

        Ok(Daf {
            path: path_str,
            id_word,
            internal_name,
            nd,
            ni,
            byte_order: order,
            summaries,
            storage,
            n_words,
        })
    }

    /// The whole file as double precision words. DAF address `a` (1-based) is `words()[a - 1]`.
    #[inline(always)]
    pub fn words(&self) -> &[f64] {
        match &self.storage {
            Storage::Mapped(m) => {
                // SAFETY: alignment was checked when the map was created, the length is a
                // whole number of f64s, and any bit pattern is a valid f64.
                unsafe { std::slice::from_raw_parts(m.as_ptr() as *const f64, self.n_words) }
            }
            Storage::Owned(v) => v,
        }
    }

    /// Words at DAF addresses `begin..=end` (1-based, inclusive), bounds checked.
    pub fn range(&self, begin: usize, end: usize) -> Option<&[f64]> {
        let w = self.words();
        if begin == 0 || end < begin || end > w.len() {
            return None;
        }
        Some(&w[begin - 1..end])
    }

    /// Read the comment area as text.
    pub fn comments(&self) -> String {
        // Comment records live between the file record and the first summary record.
        let owned;
        let bytes: &[u8] = match &self.storage {
            Storage::Mapped(m) => m,
            Storage::Owned(_) => {
                owned = match std::fs::read(&self.path) {
                    Ok(b) => b,
                    Err(_) => return String::new(),
                };
                &owned
            }
        };
        let fward = read_i32(&bytes[76..80], self.byte_order).max(2) as usize;
        let mut out = String::new();
        for rec in 1..(fward - 1) {
            let base = rec * RECORD_BYTES;
            if base + 1000 > bytes.len() {
                break;
            }
            for &b in &bytes[base..base + 1000] {
                match b {
                    0 => out.push('\n'),
                    4 => return out,
                    32..=126 => out.push(b as char),
                    _ => {}
                }
            }
        }
        out
    }
}
