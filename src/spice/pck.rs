//! Binary PCK (orientation) segments: types 2, 3 and 20.
//!
//! Each segment gives, for a body-fixed frame (identified by its "class ID"), the rotation
//! from an inertial reference frame to the body-fixed frame and its time derivative.

use super::daf::DafSummary;
use super::error::{Result, SpiceError};
use super::generic::GenericSegment;
use super::math::{self, M3};
use super::spk::{cheb20_state, file_count, file_count_trunc};

#[derive(Debug, Clone)]
pub(crate) enum PckData {
    Chebyshev2 {
        init: f64,
        intlen: f64,
        rsize: usize,
        n: usize,
    },
    Generic3 {
        gs: GenericSegment,
        ncoeff: usize,
    },
    Cheb20 {
        dscale: f64,
        tscale: f64,
        init: f64,
        initjd: f64,
        initfr: f64,
        intlen_days: f64,
        rsize: usize,
        n: usize,
    },
}

/// One binary PCK segment.
#[derive(Debug, Clone)]
pub struct PckSegment {
    /// Frame class ID (e.g. 3000 for ITRF93).
    pub class_id: i32,
    /// Inertial reference frame the rotation is relative to (e.g. 17 = ECLIPJ2000).
    pub reference_frame: i32,
    pub data_type: i32,
    pub start_et: f64,
    pub end_et: f64,
    pub name: String,
    pub(crate) begin: usize,
    pub(crate) data: PckData,
}

impl PckSegment {
    pub(crate) fn from_summary(s: &DafSummary, words: &[f64]) -> Result<PckSegment> {
        let start_et = s.doubles[0];
        let end_et = s.doubles[1];
        let class_id = s.ints[0];
        let reference_frame = s.ints[1];
        let data_type = s.ints[2];
        let begin = s.ints[3].max(0) as usize;
        let end = s.ints[4].max(0) as usize;
        let malformed = |reason: &str| SpiceError::MalformedSegment {
            kind: "PCK",
            data_type,
            body: class_id,
            reason: reason.to_string(),
        };
        if begin == 0 || end < begin || end > words.len() {
            return Err(malformed("addresses out of range"));
        }
        let len = end - begin + 1;
        let w = |a: usize| words[a - 1];
        let data = match data_type {
            2 => {
                if len < 4 {
                    return Err(malformed("segment too short"));
                }
                let init = w(end - 3);
                let intlen = w(end - 2);
                let rsize = file_count_trunc(w(end - 1)).ok_or_else(|| malformed("invalid record size"))?;
                let n = file_count_trunc(w(end)).ok_or_else(|| malformed("invalid record count"))?;
                let fits = n.checked_mul(rsize).and_then(|v| v.checked_add(4)).is_some_and(|v| v <= len);
                if n == 0 || rsize < 5 || (rsize - 2) % 3 != 0 || !fits || !(intlen > 0.0) || (rsize - 2) / 3 > math::MAX_CHEB {
                    return Err(malformed("inconsistent record layout"));
                }
                PckData::Chebyshev2 { init, intlen, rsize, n }
            }
            3 => {
                let gs = GenericSegment::parse(words, begin, end, "PCK", data_type, class_id)?;
                let ncoeff = gs
                    .constant(words, 1)
                    .and_then(file_count)
                    .ok_or_else(|| malformed("missing or invalid coefficient count"))?;
                if ncoeff == 0 || ncoeff > math::MAX_CHEB {
                    return Err(malformed("bad coefficient count"));
                }
                PckData::Generic3 { gs, ncoeff }
            }
            20 => {
                if len < 7 {
                    return Err(malformed("segment too short"));
                }
                let dscale = w(end - 6);
                let tscale = w(end - 5);
                let initjd = w(end - 4);
                let initfr = w(end - 3);
                let intlen_days = w(end - 2);
                let rsize = file_count_trunc(w(end - 1)).ok_or_else(|| malformed("invalid record size"))?;
                let n = file_count_trunc(w(end)).ok_or_else(|| malformed("invalid record count"))?;
                let fits = n.checked_mul(rsize).and_then(|v| v.checked_add(7)).is_some_and(|v| v <= len);
                if n == 0 || rsize < 6 || rsize % 3 != 0 || !fits || rsize / 3 > math::MAX_CHEB + 1 || !(intlen_days > 0.0) {
                    return Err(malformed("inconsistent record layout"));
                }
                PckData::Cheb20 {
                    dscale,
                    tscale,
                    init: (initjd - 2_451_545.0 + initfr) * 86_400.0,
                    initjd,
                    initfr,
                    intlen_days,
                    rsize,
                    n,
                }
            }
            other => {
                return Err(SpiceError::UnsupportedSegmentType {
                    kind: "PCK",
                    data_type: other,
                    body: class_id,
                })
            }
        };
        Ok(PckSegment {
            class_id,
            reference_frame,
            data_type,
            start_et,
            end_et,
            name: s.name.clone(),
            begin,
            data,
        })
    }

    #[inline(always)]
    pub fn covers(&self, et: f64) -> bool {
        et >= self.start_et && et <= self.end_et
    }

    /// Rotation from the reference frame to the body-fixed frame, and its time derivative.
    pub(crate) fn rotation(&self, words: &[f64], et: f64) -> Result<(M3, M3)> {
        match &self.data {
            PckData::Chebyshev2 { init, intlen, rsize, n } => {
                let idx = (((et - init) / intlen) as i64).clamp(0, *n as i64 - 1) as usize;
                let base = self.begin - 1 + idx * rsize;
                let rec = &words[base..base + rsize];
                let (mid, radius) = (rec[0], rec[1]);
                let nc = (rsize - 2) / 3;
                let (ang, rate) = math::cheb3_val_der(&rec[2..2 + 3 * nc], nc, mid, radius, et);
                let w = math::fmod(ang[2], std::f64::consts::TAU);
                Ok(math::euler_with_rate([w, ang[1], ang[0]], [rate[2], rate[1], rate[0]], [3, 1, 3]))
            }
            PckData::Cheb20 {
                dscale,
                tscale,
                init,
                initjd,
                initfr,
                intlen_days,
                rsize,
                n,
            } => {
                let s = cheb20_state(
                    words, self.begin, *dscale, *tscale, *init, *initjd, *initfr, *intlen_days, *rsize, *n, et,
                );
                let w = math::fmod(s[2], std::f64::consts::TAU);
                Ok(math::euler_with_rate([w, s[1], s[0]], [s[5], s[4], s[3]], [3, 1, 3]))
            }
            PckData::Generic3 { gs, ncoeff } => {
                let bad = |r: &str| SpiceError::MalformedSegment {
                    kind: "PCK",
                    data_type: 3,
                    body: self.class_id,
                    reason: r.to_string(),
                };
                let idx = gs.find_packet(words, et).ok_or_else(|| bad("no packet for epoch"))?;
                let pkt = gs.packet(words, idx).ok_or_else(|| bad("packet out of range"))?;
                if pkt.len() < 2 + 6 * ncoeff {
                    return Err(bad("short packet"));
                }
                let (mid, radius) = (pkt[0], pkt[1]);
                let rpd = std::f64::consts::PI / 180.0;
                let mut e = [0.0; 6];
                for c in 0..6 {
                    e[c] = rpd * math::cheb_val(&pkt[2 + c * ncoeff..2 + (c + 1) * ncoeff], mid, radius, et);
                }
                let ra = std::f64::consts::FRAC_PI_2 + e[0];
                let dec = std::f64::consts::FRAC_PI_2 - e[1];
                let rot = math::mxm(&math::rotate(e[2], 3), &math::mxm(&math::rotate(dec, 1), &math::rotate(ra, 3)));
                let mav = [-e[3], -e[4], -e[5]];
                let mut drot = [[0.0; 3]; 3];
                for j in 0..3 {
                    let col = [rot[0][j], rot[1][j], rot[2][j]];
                    let c = math::cross(&mav, &col);
                    for i in 0..3 {
                        drot[i][j] = c[i];
                    }
                }
                Ok((rot, drot))
            }
        }
    }
}
