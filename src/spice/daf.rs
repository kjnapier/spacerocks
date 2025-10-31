use std::fs::File;
use std::io::{self, Read};
use std::mem;
use std::path::Path;

use memmap2::Mmap;

/// Errors that can occur during DAF file reading and parsing.
#[derive(Debug)]
pub enum DAFError {
    IoError(io::Error),
    ParseError(String),
}

impl From<io::Error> for DAFError {
    fn from(err: io::Error) -> Self {
        DAFError::IoError(err)
    }
}
/// https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/daf.html
/// https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/req/spk.html

/// Each physical record in a DAF file is 1024 bytes long.
const RECORD_LENGTH: usize = 1024;

/// The DAF file record (header) as described in the standard.
///
/// This structure corresponds to the following fields:
/// 1. LOCIDW   (8 bytes): Identification word (e.g., "DAF/PCK")
/// 2. ND       (4 bytes): Number of double precision components per summary.
/// 3. NI       (4 bytes): Number of integer components per summary.
/// 4. LOCIFN   (60 bytes): Internal file name or description.
/// 5. FWARD    (4 bytes): Record number of the first summary record.
/// 6. BWARD    (4 bytes): Record number of the final summary record.
/// 7. FREE     (4 bytes): First free address in the file.
/// 8. LOCFMT   (8 bytes): Numeric binary format string ("LTL-IEEE" or "BIG-IEEE").
/// 9. PRENUL   (603 bytes): Padding null block.
/// 10. FTPSTR  (28 bytes): FTP validation string.
/// 11. PSTNUL  (297 bytes): Padding null block.
#[repr(C)]
#[derive(Copy, Clone, Debug)]
pub struct FileRecord {
    pub locidw: [u8; 8],
    pub nd: i32,
    pub ni: i32,
    pub locifn: [u8; 60],
    pub fward: i32,
    pub bward: i32,
    pub free: i32,
    pub locfmt: [u8; 8],
    pub prenul: [u8; 603],
    pub ftpstr: [u8; 28],
    pub pstnul: [u8; 297],
}

/// Unsafe helper function: Converts a byte slice into a struct of type `T`.
///
/// # Safety
/// This function performs an unaligned read and assumes the bytes exactly match the layout of `T`.
unsafe fn ptr_from_bytes<T: Copy>(bytes: &[u8]) -> Result<T, DAFError> {
    if bytes.len() < mem::size_of::<T>() {
        return Err(DAFError::ParseError("Buffer too small".to_string()));
    }
    let ptr = bytes.as_ptr() as *const T;
    Ok(ptr.read_unaligned())
}

/// A generic DAF reader that reads the file header and maps the file into memory.
pub struct DAFReader {
    pub file_record: FileRecord,
    pub mmap: Mmap,
}

impl DAFReader {
    /// Opens a DAF file from the given path, reads the header (first 1024 bytes),
    /// and parses it into a FileRecord.
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self, DAFError> {
        let mut file = File::open(&path)?;

        // Read the file record (the first 1024 bytes).
        let mut header_buf = [0u8; RECORD_LENGTH];
        file.read_exact(&mut header_buf)?;
        let file_record: FileRecord = unsafe { ptr_from_bytes(&header_buf)? };

        // Memory-map the entire file.
        let mmap = unsafe { Mmap::map(&file)? };

        Ok(DAFReader { file_record, mmap })
    }

    /// Returns the identification word (LOCIDW) as a UTF-8 string.
    pub fn id_word(&self) -> String {
        String::from_utf8_lossy(&self.file_record.locidw).to_string()
    }

    /// Returns the internal file name (LOCIFN) as a UTF-8 string.
    pub fn internal_file_name(&self) -> String {
        String::from_utf8_lossy(&self.file_record.locifn).to_string()
    }

    /// Returns the numeric binary format string (LOCFMT) as a UTF-8 string.
    pub fn binary_format(&self) -> String {
        String::from_utf8_lossy(&self.file_record.locfmt).to_string()
    }
}