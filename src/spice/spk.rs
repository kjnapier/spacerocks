//! SPK (ephemeris) segments: parsing and evaluation for all commonly used data types.
//!
//! Supported types: 1, 2, 3, 5, 8, 9, 12, 13, 14, 15, 17, 18, 19, 20, 21.
//! (Types 10 — two-line elements — and 16 are not supported.)
//!
//! Evaluation returns the state of the segment's target relative to its center, in the
//! segment's reference frame, in km and km/s. Epochs are TDB seconds past J2000 (ET).
//! Algorithms follow the corresponding SPICELIB `SPKRnn`/`SPKEnn` routines.

use super::daf::DafSummary;
use super::error::{Result, SpiceError};
use super::generic::GenericSegment;
use super::math::{self, MAX_WINDOW};

/// Per-type metadata extracted once when the segment is loaded.
#[derive(Debug, Clone)]
pub(crate) enum SpkData {
    /// Types 2 and 3: fixed-length Chebyshev records.
    Chebyshev {
        init: f64,
        intlen: f64,
        rsize: usize,
        n: usize,
        with_velocity: bool,
    },
    /// Types 1 and 21: modified difference arrays.
    Mda {
        maxdim: usize,
        n: usize,
        epochs: usize, // 0-based word index of first epoch
    },
    /// Type 5: discrete states with two-body propagation.
    TwoBody { n: usize, gm: f64, epochs: usize },
    /// Types 8 and 12: equally spaced discrete states.
    EqualStates {
        start: f64,
        step: f64,
        degree: usize,
        n: usize,
        hermite: bool,
    },
    /// Types 9 and 13: unequally spaced discrete states.
    UnequalStates {
        degree: usize,
        n: usize,
        epochs: usize,
        hermite: bool,
    },
    /// Type 14: Chebyshev, unequal time steps (generic segment).
    Generic14 { gs: GenericSegment, ncoeff: usize },
    /// Type 15: precessing conic.
    Precessing,
    /// Type 17: equinoctial elements.
    Equinoctial,
    /// Type 18: a single ESOC/DDID mini-segment.
    Type18 { subtype: u8, window: usize, n: usize },
    /// Type 19: list of type-18-like mini-segments.
    Type19 { n: usize, select_last: bool },
    /// Type 20: Chebyshev velocity only.
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

/// One SPK segment.
#[derive(Debug, Clone)]
pub struct SpkSegment {
    pub target: i32,
    pub center: i32,
    pub frame: i32,
    pub data_type: i32,
    pub start_et: f64,
    pub end_et: f64,
    /// Segment identifier (the DAF array name).
    pub name: String,
    pub(crate) begin: usize,
    pub(crate) end: usize,
    pub(crate) data: SpkData,
}

#[inline(always)]
fn w(words: &[f64], addr: usize) -> f64 {
    words[addr - 1]
}

/// A non-negative integer count/size read from a file, rounded to nearest; rejects negative,
/// non-finite or absurdly large values so later index arithmetic cannot overflow.
#[inline]
pub(crate) fn file_count(x: f64) -> Option<usize> {
    if x.is_finite() && (0.0..=1.0e12).contains(&x) {
        Some(x.round() as usize)
    } else {
        None
    }
}

/// Like [`file_count`] but truncating (SPICE reads some sizes with `INT` rather than `NINT`).
#[inline]
pub(crate) fn file_count_trunc(x: f64) -> Option<usize> {
    if x.is_finite() && (0.0..=1.0e12).contains(&x) {
        Some(x as usize)
    } else {
        None
    }
}

/// Number of values in `sorted` strictly less than `x`.
#[inline(always)]
fn count_lt(sorted: &[f64], x: f64) -> usize {
    sorted.partition_point(|&t| t < x)
}

impl SpkSegment {
    pub(crate) fn from_summary(s: &DafSummary, words: &[f64]) -> Result<SpkSegment> {
        let start_et = s.doubles[0];
        let end_et = s.doubles[1];
        let target = s.ints[0];
        let center = s.ints[1];
        let frame = s.ints[2];
        let data_type = s.ints[3];
        let begin = s.ints[4].max(0) as usize;
        let end = s.ints[5].max(0) as usize;

        let malformed = |reason: String| SpiceError::MalformedSegment {
            kind: "SPK",
            data_type,
            body: target,
            reason,
        };
        if begin == 0 || end < begin || end > words.len() {
            return Err(malformed(format!("addresses {}..{} out of range", begin, end)));
        }
        let len = end - begin + 1;
        let need = |k: usize| -> Result<()> {
            if len < k {
                Err(malformed(format!("segment has {} words, need at least {}", len, k)))
            } else {
                Ok(())
            }
        };
        let rint = |x: f64| x.round() as i64;
        let count = |x: f64, what: &str| file_count(x).ok_or_else(|| malformed(format!("invalid {} ({})", what, x)));
        let count_t = |x: f64, what: &str| file_count_trunc(x).ok_or_else(|| malformed(format!("invalid {} ({})", what, x)));
        // `a * b + c <= len`, without overflow.
        let fits = |a: usize, b: usize, c: usize| a.checked_mul(b).and_then(|v| v.checked_add(c)).is_some_and(|v| v <= len);

        let data = match data_type {
            2 | 3 => {
                need(4)?;
                let init = w(words, end - 3);
                let intlen = w(words, end - 2);
                let rsize = count_t(w(words, end - 1), "record size")?;
                let n = count_t(w(words, end), "record count")?;
                let per = if data_type == 2 { 3 } else { 6 };
                if n == 0 || rsize < 2 + per || (rsize - 2) % per != 0 || !fits(n, rsize, 4) || !(intlen > 0.0) || !init.is_finite() {
                    return Err(malformed("inconsistent Chebyshev record layout".into()));
                }
                if (rsize - 2) / per > math::MAX_CHEB {
                    return Err(malformed("too many Chebyshev coefficients".into()));
                }
                SpkData::Chebyshev {
                    init,
                    intlen,
                    rsize,
                    n,
                    with_velocity: data_type == 3,
                }
            }
            1 | 21 => {
                let (maxdim, n, epochs_end) = if data_type == 1 {
                    need(1)?;
                    let n = count_t(w(words, end), "record count")?;
                    let ndir = n / 100;
                    (15usize, n, end.checked_sub(ndir + 1)) // last epoch address
                } else {
                    need(2)?;
                    let maxdim = count(w(words, end - 1), "MAXDIM")?;
                    let n = count(w(words, end), "record count")?;
                    let ndir = n / 100;
                    (maxdim, n, end.checked_sub(ndir + 2))
                };
                if maxdim == 0 || maxdim > 25 {
                    return Err(malformed(format!("MAXDIM {} out of range", maxdim)));
                }
                let dflsiz = 4 * maxdim + 11;
                let epochs_end = epochs_end.unwrap_or(0);
                if n == 0 || !fits(n, dflsiz + 1, 0) || epochs_end < n || epochs_end + 1 - n < begin {
                    return Err(malformed("inconsistent MDA layout".into()));
                }
                let first_epoch_addr = epochs_end + 1 - n;
                SpkData::Mda {
                    maxdim,
                    n,
                    epochs: first_epoch_addr - 1,
                }
            }
            5 => {
                need(2)?;
                let gm = w(words, end - 1);
                let n = count(w(words, end), "state count")?;
                if n == 0 || !fits(n, 7, 2) {
                    return Err(malformed("inconsistent type 5 layout".into()));
                }
                SpkData::TwoBody {
                    n,
                    gm,
                    epochs: begin + 6 * n - 1,
                }
            }
            8 | 12 => {
                need(4)?;
                let start = w(words, end - 3);
                let step = w(words, end - 2);
                let degree = count(w(words, end - 1), "degree")?;
                let n = count(w(words, end), "state count")?;
                if n == 0 || degree + 1 > MAX_WINDOW || degree + 1 > n || !fits(n, 6, 4) || !(step > 0.0) {
                    return Err(malformed("inconsistent type 8/12 layout".into()));
                }
                SpkData::EqualStates {
                    start,
                    step,
                    degree,
                    n,
                    hermite: data_type == 12,
                }
            }
            9 | 13 => {
                need(2)?;
                let degree = count(w(words, end - 1), "degree")?;
                let n = count(w(words, end), "state count")?;
                if n == 0 || degree + 1 > MAX_WINDOW || degree + 1 > n || !fits(n, 7, 2) {
                    return Err(malformed("inconsistent type 9/13 layout".into()));
                }
                SpkData::UnequalStates {
                    degree,
                    n,
                    epochs: begin + 6 * n - 1,
                    hermite: data_type == 13,
                }
            }
            14 => {
                let gs = GenericSegment::parse(words, begin, end, "SPK", data_type, target)?;
                let ncoeff = count(
                    gs.constant(words, 1).ok_or_else(|| malformed("missing coefficient count".into()))?,
                    "coefficient count",
                )?;
                if ncoeff == 0 || ncoeff > math::MAX_CHEB {
                    return Err(malformed("bad coefficient count".into()));
                }
                SpkData::Generic14 { gs, ncoeff }
            }
            15 => {
                if len != 16 {
                    return Err(malformed(format!("type 15 segment should have 16 words, has {}", len)));
                }
                SpkData::Precessing
            }
            17 => {
                if len != 12 {
                    return Err(malformed(format!("type 17 segment should have 12 words, has {}", len)));
                }
                SpkData::Equinoctial
            }
            18 => {
                need(3)?;
                let subtype = rint(w(words, end - 2));
                let window = count(w(words, end - 1), "window size")?;
                let n = count(w(words, end), "packet count")?;
                if !(0..=1).contains(&subtype) {
                    return Err(malformed(format!("unknown type 18 subtype {}", subtype)));
                }
                let packsz = if subtype == 0 { 12 } else { 6 };
                if n < 2 || !(2..=MAX_WINDOW).contains(&window) || !fits(n, packsz + 1, 3 + (n - 1) / 100) {
                    return Err(malformed("bad type 18 window/size".into()));
                }
                SpkData::Type18 {
                    subtype: subtype as u8,
                    window,
                    n,
                }
            }
            19 => {
                need(2)?;
                let isel = rint(w(words, end - 1));
                let n = count(w(words, end), "interval count")?;
                if n == 0 || !fits(n + 1, 2, n / 100 + 2) {
                    return Err(malformed("bad type 19 interval count".into()));
                }
                SpkData::Type19 {
                    n,
                    select_last: isel == 1,
                }
            }
            20 => {
                need(7)?;
                let dscale = w(words, end - 6);
                let tscale = w(words, end - 5);
                let initjd = w(words, end - 4);
                let initfr = w(words, end - 3);
                let intlen_days = w(words, end - 2);
                let rsize = count_t(w(words, end - 1), "record size")?;
                let n = count_t(w(words, end), "record count")?;
                if n == 0 || rsize < 6 || rsize % 3 != 0 || !fits(n, rsize, 7) || rsize / 3 > math::MAX_CHEB + 1 || !(intlen_days > 0.0) {
                    return Err(malformed("inconsistent type 20 layout".into()));
                }
                let init = (initjd - 2_451_545.0 + initfr) * 86_400.0;
                SpkData::Cheb20 {
                    dscale,
                    tscale,
                    init,
                    initjd,
                    initfr,
                    intlen_days,
                    rsize,
                    n,
                }
            }
            other => {
                return Err(SpiceError::UnsupportedSegmentType {
                    kind: "SPK",
                    data_type: other,
                    body: target,
                })
            }
        };

        Ok(SpkSegment {
            target,
            center,
            frame,
            data_type,
            start_et,
            end_et,
            name: s.name.clone(),
            begin,
            end,
            data,
        })
    }

    /// Does this segment cover `et` (inclusive at both ends, as in SPICE)?
    #[inline(always)]
    pub fn covers(&self, et: f64) -> bool {
        et >= self.start_et && et <= self.end_et
    }

    fn malformed(&self, reason: &str) -> SpiceError {
        SpiceError::MalformedSegment {
            kind: "SPK",
            data_type: self.data_type,
            body: self.target,
            reason: reason.to_string(),
        }
    }

    /// State (km, km/s) of `target` relative to `center` in the segment frame at `et`.
    pub(crate) fn state(&self, words: &[f64], et: f64) -> Result<[f64; 6]> {
        match &self.data {
            SpkData::Chebyshev {
                init,
                intlen,
                rsize,
                n,
                with_velocity,
            } => {
                let mut idx = ((et - init) / intlen) as i64;
                idx = idx.clamp(0, *n as i64 - 1);
                let base = self.begin - 1 + idx as usize * rsize;
                let rec = &words[base..base + rsize];
                let mid = rec[0];
                let radius = rec[1];
                let mut out = [0.0; 6];
                if !with_velocity {
                    let nc = (rsize - 2) / 3;
                    let (p, v) = math::cheb3_val_der(&rec[2..2 + 3 * nc], nc, mid, radius, et);
                    out[..3].copy_from_slice(&p);
                    out[3..].copy_from_slice(&v);
                } else {
                    let nc = (rsize - 2) / 6;
                    let p = math::cheb3_val(&rec[2..2 + 3 * nc], nc, mid, radius, et);
                    let v = math::cheb3_val(&rec[2 + 3 * nc..2 + 6 * nc], nc, mid, radius, et);
                    out[..3].copy_from_slice(&p);
                    out[3..].copy_from_slice(&v);
                }
                Ok(out)
            }

            SpkData::Mda { maxdim, n, epochs } => {
                let ep = &words[*epochs..*epochs + n];
                let recno = count_lt(ep, et).min(n - 1); // 0-based: first epoch >= et
                let dflsiz = 4 * maxdim + 11;
                let base = self.begin - 1 + recno * dflsiz;
                mda_state(&words[base..base + dflsiz], *maxdim, et)
                    .ok_or_else(|| self.malformed("invalid difference-line record"))
            }

            SpkData::TwoBody { n, gm, epochs } => {
                let ep = &words[*epochs..*epochs + n];
                let i = count_lt(ep, et); // number of epochs strictly before et
                let st = |k: usize| -> [f64; 6] {
                    let b = self.begin - 1 + 6 * k;
                    [words[b], words[b + 1], words[b + 2], words[b + 3], words[b + 4], words[b + 5]]
                };
                let (k1, k2) = if i == 0 {
                    (0, 0)
                } else if i >= *n {
                    (n - 1, n - 1)
                } else {
                    (i - 1, i)
                };
                let t1 = ep[k1];
                let t2 = ep[k2];
                let fail = || self.malformed("two-body propagation failed");
                if k1 == k2 {
                    return math::prop2b(*gm, &st(k1), et - t1).ok_or_else(fail);
                }
                let s1 = math::prop2b(*gm, &st(k1), et - t1).ok_or_else(fail)?;
                let s2 = math::prop2b(*gm, &st(k2), et - t2).ok_or_else(fail)?;
                let denom = t2 - t1;
                let arg = (et - t1) * std::f64::consts::PI / denom;
                let dargdt = std::f64::consts::PI / denom;
                let wgt = 0.5 * arg.cos() + 0.5;
                let dwdt = -0.5 * arg.sin() * dargdt;
                let mut out = [0.0; 6];
                for c in 0..6 {
                    out[c] = wgt * s1[c] + (1.0 - wgt) * s2[c];
                }
                for c in 0..3 {
                    out[c + 3] += dwdt * s1[c] - dwdt * s2[c];
                }
                Ok(out)
            }

            SpkData::EqualStates {
                start,
                step,
                degree,
                n,
                hermite,
            } => {
                let size = degree + 1;
                let first = if size % 2 == 1 {
                    let near = ((et - start) / step).round() as i64 + 1;
                    (near - (*degree / 2) as i64).max(1).min((*n - degree) as i64)
                } else {
                    let low = ((et - start) / step) as i64 + 1;
                    (low - (*degree / 2) as i64).max(1).min((*n - degree) as i64)
                } as usize;
                let mut xs = [0.0f64; MAX_WINDOW];
                for k in 0..size {
                    xs[k] = start + ((first - 1 + k) as f64) * step;
                }
                let base = self.begin - 1 + (first - 1) * 6;
                Ok(interp_states(&words[base..base + 6 * size], &xs[..size], et, *hermite))
            }

            SpkData::UnequalStates {
                degree,
                n,
                epochs,
                hermite,
            } => {
                let ep = &words[*epochs..*epochs + n];
                let size = degree + 1;
                let lt = count_lt(ep, et);
                let low = lt.max(1); // 1-based index of last epoch < et (or 1)
                let high = low + 1;
                let first = if size % 2 == 1 {
                    let near = if lt == 0 || high > *n || (et - ep[low - 1]).abs() < (et - ep[high - 1]).abs() {
                        low
                    } else {
                        high
                    };
                    (near as i64 - (*degree / 2) as i64).max(1).min((*n - degree) as i64)
                } else {
                    (low as i64 - (*degree / 2) as i64).max(1).min((*n - degree) as i64)
                } as usize;
                let base = self.begin - 1 + (first - 1) * 6;
                let xs = &ep[first - 1..first - 1 + size];
                Ok(interp_states(&words[base..base + 6 * size], xs, et, *hermite))
            }

            SpkData::Generic14 { gs, ncoeff } => {
                let idx = gs
                    .find_packet(words, et)
                    .ok_or_else(|| self.malformed("no packet for epoch"))?;
                let pkt = gs
                    .packet(words, idx)
                    .ok_or_else(|| self.malformed("packet out of range"))?;
                if pkt.len() < 2 + 6 * ncoeff {
                    return Err(self.malformed("short packet"));
                }
                let mid = pkt[0];
                let radius = pkt[1];
                let mut out = [0.0; 6];
                for c in 0..6 {
                    out[c] = math::cheb_val(&pkt[2 + c * ncoeff..2 + (c + 1) * ncoeff], mid, radius, et);
                }
                Ok(out)
            }

            SpkData::Precessing => {
                let rec = &words[self.begin - 1..self.end];
                type15_state(rec, et).ok_or_else(|| self.malformed("invalid precessing conic elements"))
            }

            SpkData::Equinoctial => {
                let rec = &words[self.begin - 1..self.end];
                type17_state(rec, et).ok_or_else(|| self.malformed("invalid equinoctial elements"))
            }

            SpkData::Type18 { subtype, window, n } => {
                let packsz = if *subtype == 0 { 12 } else { 6 };
                mini_segment_state(words, self.begin, *n, packsz, *subtype, *window, et)
                    .ok_or_else(|| self.malformed("could not evaluate type 18 data"))
            }

            SpkData::Type19 { n, select_last } => {
                // Layout from the end: [N][ISEL][pointers N+1][interval dir N/100][bounds N+1]
                let ndir = n / 100;
                let ptr_base = self.end - 2 - (n + 1); // address before first pointer
                let dir_base = ptr_base - ndir;
                let iv_base = dir_base - (n + 1);
                let bounds = &words[iv_base..iv_base + n + 1];
                let mini = if *select_last {
                    // interval i with b_i <= et < b_{i+1}
                    let k = bounds.partition_point(|&b| b <= et); // count <= et
                    k.clamp(1, *n)
                } else {
                    // interval i with b_i < et <= b_{i+1}
                    let k = bounds.partition_point(|&b| b < et);
                    k.clamp(1, *n)
                };
                let bad = || self.malformed("bad mini-segment pointers");
                let p1 = file_count(words[ptr_base + mini - 1]).ok_or_else(bad)?;
                let p2 = file_count(words[ptr_base + mini]).ok_or_else(bad)?;
                if p1 == 0 || p2 < p1 + 4 {
                    return Err(bad());
                }
                let minib = p1 + self.begin - 1;
                let minie = p2 + self.begin - 2;
                if minie > self.end {
                    return Err(bad());
                }
                let subtype = words[minie - 3].round() as i64;
                let window = file_count(words[minie - 2]).ok_or_else(bad)?;
                let npkt = file_count(words[minie - 1]).ok_or_else(bad)?;
                let packsz = match subtype {
                    0 => 12,
                    1 | 2 => 6,
                    _ => return Err(self.malformed("unknown type 19 subtype")),
                };
                let needed = npkt.checked_mul(packsz + 1).map(|v| v + (npkt.saturating_sub(1)) / 100 + 3);
                if !(2..=MAX_WINDOW).contains(&window) || npkt < 1 || needed.map_or(true, |v| v > minie + 1 - minib) {
                    return Err(self.malformed("bad type 19 mini-segment window/size"));
                }
                mini_segment_state(words, minib, npkt, packsz, subtype as u8, window, et)
                    .ok_or_else(|| self.malformed("could not evaluate type 19 data"))
            }

            SpkData::Cheb20 {
                dscale,
                tscale,
                init,
                initjd,
                initfr,
                intlen_days,
                rsize,
                n,
            } => Ok(cheb20_state(
                words,
                self.begin,
                *dscale,
                *tscale,
                *init,
                *initjd,
                *initfr,
                *intlen_days,
                *rsize,
                *n,
                et,
            )),
        }
    }
}

/// Evaluate a type 20 (or PCK type 20) record: velocity Chebyshev, position by integration.
#[allow(clippy::too_many_arguments)]
pub(crate) fn cheb20_state(
    words: &[f64],
    begin: usize,
    dscale: f64,
    tscale: f64,
    init: f64,
    initjd: f64,
    initfr: f64,
    intlen_days: f64,
    rsize: usize,
    n: usize,
    et: f64,
) -> [f64; 6] {
    let intrvl = intlen_days * 86_400.0;
    let recno = (((et - init) / intrvl) as i64 + 1).clamp(1, n as i64) as usize;
    let recbeg = (initjd - 2_451_545.0 + (recno as f64 - 1.0) * intlen_days) * 86_400.0;
    let radius = intrvl / 2.0;
    let mid = recbeg + initfr * 86_400.0 + radius;
    let base = begin - 1 + (recno - 1) * rsize;
    let rec = &words[base..base + rsize];
    let nterms = rsize / 3;
    let ncof = nterms - 1;
    let vscale = dscale / tscale;
    let mut out = [0.0; 6];
    let mut buf = [0.0f64; math::MAX_CHEB + 1];
    for c in 0..3 {
        let comp = &rec[c * nterms..(c + 1) * nterms];
        for (k, v) in comp[..ncof].iter().enumerate() {
            buf[k] = v * vscale;
        }
        let pos_mid = comp[ncof] * dscale;
        let (vel, integral) = math::cheb_val_integral(&buf[..ncof], mid, radius, et);
        out[c] = pos_mid + integral;
        out[c + 3] = vel;
    }
    out
}

/// Interpolate a window of 6-component states (types 8, 9, 12, 13).
fn interp_states(states: &[f64], xs: &[f64], et: f64, hermite: bool) -> [f64; 6] {
    let m = xs.len();
    let mut out = [0.0; 6];
    let mut a = [0.0f64; MAX_WINDOW];
    let mut b = [0.0f64; MAX_WINDOW];
    if hermite {
        for c in 0..3 {
            for k in 0..m {
                a[k] = states[6 * k + c];
                b[k] = states[6 * k + c + 3];
            }
            let (v, d) = math::hermite(xs, &a[..m], &b[..m], et);
            out[c] = v;
            out[c + 3] = d;
        }
    } else {
        for c in 0..6 {
            for k in 0..m {
                a[k] = states[6 * k + c];
            }
            out[c] = math::lagrange(xs, &a[..m], et);
        }
    }
    out
}

/// Evaluate a type 18 segment or type 19 mini-segment whose packets start at DAF address `minib`.
///
/// Layout: [packets npkt*packsz][epochs npkt][directory (npkt-1)/100][subtype, window, npkt].
fn mini_segment_state(
    words: &[f64],
    minib: usize,
    npkt: usize,
    packsz: usize,
    subtype: u8,
    window: usize,
    et: f64,
) -> Option<[f64; 6]> {
    let ep_start = minib - 1 + npkt * packsz; // 0-based index of first epoch
    let ep = words.get(ep_start..ep_start + npkt)?;
    let lt = count_lt(ep, et);
    let low = lt.max(1);
    let high = low + 1;
    let lsize = (window / 2).min(low);
    let rsize = (window / 2).min((npkt + 1).saturating_sub(high));
    let count = lsize + rsize;
    if count == 0 {
        return None;
    }
    let first = low - lsize + 1; // 1-based
    let xs = &ep[first - 1..first - 1 + count];
    let pk = words.get(minib - 1 + (first - 1) * packsz..minib - 1 + (first - 1 + count) * packsz)?;

    let mut out = [0.0; 6];
    let mut a = [0.0f64; MAX_WINDOW];
    let mut b = [0.0f64; MAX_WINDOW];
    match subtype {
        0 => {
            for c in 0..3 {
                for k in 0..count {
                    a[k] = pk[packsz * k + c];
                    b[k] = pk[packsz * k + c + 3];
                }
                out[c] = math::hermite(xs, &a[..count], &b[..count], et).0;
                for k in 0..count {
                    a[k] = pk[packsz * k + c + 6];
                    b[k] = pk[packsz * k + c + 9];
                }
                out[c + 3] = math::hermite(xs, &a[..count], &b[..count], et).0;
            }
        }
        1 => {
            for c in 0..6 {
                for k in 0..count {
                    a[k] = pk[packsz * k + c];
                }
                out[c] = math::lagrange(xs, &a[..count], et);
            }
        }
        _ => {
            for c in 0..3 {
                for k in 0..count {
                    a[k] = pk[packsz * k + c];
                    b[k] = pk[packsz * k + c + 3];
                }
                let (v, d) = math::hermite(xs, &a[..count], &b[..count], et);
                out[c] = v;
                out[c + 3] = d;
            }
        }
    }
    Some(out)
}

/// Modified difference array evaluation (SPKE01 / SPKE21). `r` is one record of
/// `4*maxdim + 11` words.
fn mda_state(r: &[f64], maxdim: usize, et: f64) -> Option<[f64; 6]> {
    let tl = r[0];
    let g = &r[1..1 + maxdim];
    let refpos = [r[maxdim + 1], r[maxdim + 3], r[maxdim + 5]];
    let refvel = [r[maxdim + 2], r[maxdim + 4], r[maxdim + 6]];
    let dt = |c: usize, j: usize| r[(c + 1) * maxdim + 7 + j];
    let kqmax1 = file_count_trunc(r[4 * maxdim + 7])?;
    let kq = [
        file_count_trunc(r[4 * maxdim + 8])?,
        file_count_trunc(r[4 * maxdim + 9])?,
        file_count_trunc(r[4 * maxdim + 10])?,
    ];
    if kqmax1 < 2 || kqmax1 > maxdim + 1 || kq.iter().any(|&k| k > maxdim || k + 1 > kqmax1 + 1) {
        return None;
    }

    let mut fc = [0.0f64; 32];
    let mut wc = [0.0f64; 32];
    let mut wv = [0.0f64; 34];
    fc[0] = 1.0;

    let delta = et - tl;
    let mut tp = delta;
    let mq2 = kqmax1.saturating_sub(2);
    let mut ks = kqmax1 as i64 - 1;
    for j in 1..=mq2 {
        fc[j] = tp / g[j - 1];
        wc[j - 1] = delta / g[j - 1];
        tp = delta + g[j - 1];
    }
    for j in 1..=kqmax1 {
        wv[j - 1] = 1.0 / j as f64;
    }
    let mut jx: usize = 0;
    let mut ks1 = ks - 1;
    while ks >= 2 {
        jx += 1;
        for j in 1..=jx {
            let a = (j as i64 + ks - 1) as usize;
            let b = (j as i64 + ks1 - 1) as usize;
            wv[a] = fc[j] * wv[b] - wc[j - 1] * wv[a];
        }
        ks = ks1;
        ks1 -= 1;
    }
    let mut out = [0.0; 6];
    for c in 0..3 {
        let mut sum = 0.0;
        for j in (1..=kq[c]).rev() {
            sum += dt(c, j - 1) * wv[(j as i64 + ks - 1) as usize];
        }
        out[c] = refpos[c] + delta * (refvel[c] + delta * sum);
    }
    for j in 1..=jx {
        let a = (j as i64 + ks - 1) as usize;
        let b = (j as i64 + ks1 - 1) as usize;
        wv[a] = fc[j] * wv[b] - wc[j - 1] * wv[a];
    }
    ks -= 1;
    for c in 0..3 {
        let mut sum = 0.0;
        for j in (1..=kq[c]).rev() {
            sum += dt(c, j - 1) * wv[(j as i64 + ks - 1) as usize];
        }
        out[c + 3] = refvel[c] + delta * sum;
    }
    Some(out)
}

/// Precessing conic propagation (SPKE15).
fn type15_state(rec: &[f64], et: f64) -> Option<[f64; 6]> {
    use math::{cross, dot, norm, vrotv, vsep};
    let epoch = rec[0];
    let mut tp = [rec[1], rec[2], rec[3]];
    let mut pa = [rec[4], rec[5], rec[6]];
    let p = rec[7];
    let ecc = rec[8];
    let j2flg = rec[9] as i64;
    let mut pv = [rec[10], rec[11], rec[12]];
    let gm = rec[13];
    let oj2 = rec[14];
    let rpl = rec[15];
    if p <= 0.0 || ecc < 0.0 || gm <= 0.0 || rpl < 0.0 {
        return None;
    }
    for v in [&mut pa, &mut tp, &mut pv] {
        let n = norm(v);
        if n == 0.0 {
            return None;
        }
        v.iter_mut().for_each(|x| *x /= n);
    }
    if dot(&pa, &tp).abs() > 1e-5 {
        return None;
    }
    let near = p / (1.0 + ecc);
    let speed = (gm / p).sqrt() * (1.0 + ecc);
    let vdir = cross(&tp, &pa);
    let state0 = [
        near * pa[0],
        near * pa[1],
        near * pa[2],
        speed * vdir[0],
        speed * vdir[1],
        speed * vdir[2],
    ];
    let dt = et - epoch;
    let mut state = math::prop2b(gm, &state0, dt)?;
    if j2flg != 3 && oj2 != 0.0 && ecc < 1.0 && near > rpl {
        let oneme2 = 1.0 - ecc * ecc;
        let dmdt = oneme2 / p * (gm * oneme2 / p).sqrt();
        let manom = dmdt * dt;
        let two_pi = std::f64::consts::TAU;
        let mut theta = math::fmod(manom, two_pi);
        if theta.abs() > std::f64::consts::PI {
            theta -= two_pi.copysign(theta);
        }
        let k2pi = manom - theta;
        let pos = [state[0], state[1], state[2]];
        let mut ta = vsep(&pa, &pos);
        ta = ta.abs().copysign(theta);
        ta += k2pi;
        let cosinc = dot(&pv, &tp);
        let z = ta * 1.5 * oj2 * (rpl / p) * (rpl / p);
        let dnode = -z * cosinc;
        let dperi = z * (2.5 * cosinc * cosinc - 0.5);
        if j2flg != 1 {
            let r = vrotv(&[state[0], state[1], state[2]], &tp, dperi);
            let v = vrotv(&[state[3], state[4], state[5]], &tp, dperi);
            state = [r[0], r[1], r[2], v[0], v[1], v[2]];
        }
        if j2flg != 2 {
            let r = vrotv(&[state[0], state[1], state[2]], &pv, dnode);
            let v = vrotv(&[state[3], state[4], state[5]], &pv, dnode);
            state = [r[0], r[1], r[2], v[0], v[1], v[2]];
        }
    }
    Some(state)
}

/// Equinoctial elements propagation (SPKE17 / EQNCPV).
fn type17_state(rec: &[f64], et: f64) -> Option<[f64; 6]> {
    let epoch = rec[0];
    let eqel = &rec[1..10];
    let rapol = rec[10];
    let decpol = rec[11];
    let a = eqel[0];
    if a <= 0.0 || (eqel[1] * eqel[1] + eqel[2] * eqel[2]).sqrt() > 0.9 {
        return None;
    }
    let (sa, ca) = rapol.sin_cos();
    let (sd, cd) = decpol.sin_cos();
    // Column-major TRANS from EQNCPV, written here row-major.
    let trans = [[-sa, -ca * sd, ca * cd], [ca, -sa * sd, sa * cd], [0.0, cd, sd]];

    let dt = et - epoch;
    let dlpdt = eqel[6];
    let dlp = dt * dlpdt;
    let (san, can) = dlp.sin_cos();
    let h = eqel[1] * can + eqel[2] * san;
    let k = eqel[2] * can - eqel[1] * san;
    let l = eqel[3];
    let nodedt = eqel[8];
    let node = dt * nodedt;
    let (sn, cn) = node.sin_cos();
    let p = eqel[4] * cn + eqel[5] * sn;
    let q = eqel[5] * cn - eqel[4] * sn;
    let mldt = eqel[7];
    let prate = dlpdt - nodedt;
    let mut b = (1.0 - h * h - k * k).sqrt();
    b = 1.0 / (1.0 + b);
    let di = 1.0 / (1.0 + p * p + q * q);
    let vf = [(1.0 - p * p + q * q) * di, 2.0 * p * q * di, -2.0 * p * di];
    let vg = [2.0 * p * q * di, (1.0 + p * p - q * q) * di, 2.0 * q * di];

    let ml = l + math::fmod(mldt * dt, std::f64::consts::TAU);
    let eecan = math::kepleq(ml, h, k);
    let (sf, cf) = eecan.sin_cos();
    let x1 = a * ((1.0 - b * h * h) * cf + (h * k * b * sf - k));
    let y1 = a * ((1.0 - b * k * k) * sf + (h * k * b * cf - h));
    let rb = h * sf + k * cf;
    let r = a * (1.0 - rb);
    let ra = mldt * a * a / r;
    let dx1 = ra * (-sf + h * b * rb);
    let dy1 = ra * (cf - k * b * rb);
    let nfac = 1.0 - dlpdt / mldt;
    let dx = nfac * dx1 - prate * y1;
    let dy = nfac * dy1 + prate * x1;

    let pos = [x1 * vf[0] + y1 * vg[0], x1 * vf[1] + y1 * vg[1], x1 * vf[2] + y1 * vg[2]];
    let temp = [-nodedt * pos[1], nodedt * pos[0], 0.0];
    let vel = [
        temp[0] + dx * vf[0] + dy * vg[0],
        temp[1] + dx * vf[1] + dy * vg[1],
        temp[2] + dx * vf[2] + dy * vg[2],
    ];
    let p_out = math::mxv(&trans, &pos);
    let v_out = math::mxv(&trans, &vel);
    Some([p_out[0], p_out[1], p_out[2], v_out[0], v_out[1], v_out[2]])
}
