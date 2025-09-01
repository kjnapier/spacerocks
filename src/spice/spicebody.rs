use crate::spice::spk::Spk;
use crate::spice::error::SpiceError;
use crate::spice::spicekernel::SpiceKernel;

use crate::constants::MASSES;

use std::collections::HashMap;
use std::cell::Cell;

use lazy_static::lazy_static;
use phf;

lazy_static! {
    static ref SPICE_BODY_MAP: phf::Map<&'static str, i32> = phf::phf_map! {
        "SSB" => 0,
        "SOLAR SYSTEM BARYCENTER" => 0,
        "SUN" => 10,
        "MERCURY" => 199,
        "VENUS" => 299,
        "EARTH" => 399,
        "MOON" => 301,
        "MERCURY BARYCENTER" => 1,
        "VENUS BARYCENTER" => 2,
        "EARTH BARYCENTER" => 3,
        "MARS BARYCENTER" => 4,
        "JUPITER BARYCENTER" => 5,
        "SATURN BARYCENTER" => 6,
        "URANUS BARYCENTER" => 7,
        "NEPTUNE BARYCENTER" => 8,
        "PLUTO BARYCENTER" => 9,
        "2000001" => 2000001,
        "CERES" => 2000001,
        "2000002" => 2000002,
        "PALLAS" => 2000002,
        "2000003" => 2000003,
        "2000004" => 2000004,
        "2000007" => 2000007,
        "2000010" => 2000010,
        "2000015" => 2000015,
        "2000016" => 2000016,
        "2000031" => 2000031,
        "2000052" => 2000052,
        "2000065" => 2000065,
        "2000087" => 2000087,
        "2000088" => 2000088,
        "2000107" => 2000107,
        "2000511" => 2000511,
        "2000704" => 2000704,
        "HST" => -48
    };
}

// #[derive(Debug, Clone)]
// pub struct SpiceBody {
//     pub code: i32,
//     pub name: String,
//     pub gm: f64,
//     pub mass: f64,
//     // Cached index for the target.
//     cached_index: Cell<Option<usize>>,
//     cached_spk_index: Cell<Option<usize>>,
// }

use std::sync::Mutex;

#[derive(Debug, Clone)]
pub struct SpiceBody {
    pub code: i32,
    pub name: String,
    pub gm: f64,
    pub mass: f64,
    // Cached index for the target.
    cached_index: Cell<Option<usize>>,
    cached_spk_index: Cell<Option<usize>>,
}


impl SpiceBody {
    pub fn new(code: i32, name: &str, gm: f64) -> Self {
        let mass = MASSES.get(&name.to_string().to_lowercase()).cloned().unwrap_or(0.0);
        SpiceBody {
            code,
            name: name.to_string(),
            gm,
            mass: mass, // TODO: Get the mass value from the SPK file.
            cached_index: Cell::new(None),
            cached_spk_index: Cell::new(None),
        }
    }

    pub fn from_name(name: &str) -> Result<Self, Box<dyn std::error::Error>> {

        // first get name to all upper case
        let name = name.to_uppercase();

        let spiceid = SPICE_BODY_MAP.get(&name).ok_or_else(|| SpiceError::BodyNotFound(name.to_string()))?;
        let gm = 0.0; // TODO: Get the GM value from the SPK file.
        Ok(SpiceBody::new(*spiceid, &name, gm))

    }

    // Get the target index with caching.
    pub fn get_index(&self, spk: &Spk) -> Option<usize> {

        if let Some(index) = self.cached_index.get() {
            Some(index)
        } else {
            // Lookup in the target_map.
            let index = spk.target_map.get(&self.code).cloned();
            if let Some(i) = index {
                // Cache the value so that future lookups are faster.
                self.cached_index.set(Some(i));
            }
            index
        }

    }

}

unsafe impl Sync for SpiceBody {}