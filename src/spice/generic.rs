//! NAIF "generic segment" layout, used by SPK type 14 and binary PCK type 3.

use super::error::{Result, SpiceError};

/// Parsed generic-segment metadata. All `*_base` values are absolute DAF addresses of the
/// word *before* the first element (element `k`, 1-based, is at `base + k`).
#[derive(Debug, Clone)]
pub struct GenericSegment {
    pub const_base: usize,
    pub n_const: usize,
    pub ref_base: usize,
    pub n_ref: usize,
    pub ref_type: i32,
    pub pkt_dir_type: i32,
    pub pkt_dir_base: usize,
    pub pkt_base: usize,
    pub n_pkt: usize,
    pub pkt_size: usize,
    pub pkt_offset: usize,
    /// Last DAF address of the segment.
    pub end: usize,
}

impl GenericSegment {
    pub fn parse(words: &[f64], begin: usize, end: usize, kind: &'static str, data_type: i32, body: i32) -> Result<Self> {
        let bad = |reason: &str| SpiceError::MalformedSegment {
            kind,
            data_type,
            body,
            reason: reason.to_string(),
        };
        if end == 0 || end > words.len() {
            return Err(bad("segment end address out of range"));
        }
        let raw = words[end - 1];
        if !(raw.is_finite() && (0.0..=1.0e6).contains(&raw)) {
            return Err(bad("invalid metadata size"));
        }
        let metasz_raw = raw.round() as i64;
        if metasz_raw < 15 {
            return Err(bad("generic segment has fewer than 15 metadata items"));
        }
        // SGMETA's handling of old (15-item) and new (17-item) metadata layouts.
        let (metasz, ametas) = if metasz_raw == 15 {
            (16usize, 16usize)
        } else if metasz_raw > 17 {
            (17usize, metasz_raw as usize)
        } else {
            (metasz_raw as usize, metasz_raw as usize)
        };
        if ametas > end {
            return Err(bad("metadata extends past segment start"));
        }
        let begm = end - ametas + 1;
        let mut meta = [0i64; 17];
        for i in 0..metasz {
            let v = words[begm - 1 + i];
            if !(v.is_finite() && v.abs() <= 1.0e12) {
                return Err(bad("invalid metadata value"));
            }
            meta[i] = v.round() as i64;
        }
        for item in meta.iter_mut().take(16).skip(metasz.saturating_sub(1)) {
            // SGMETA zeroes items metasz..16 (only relevant for 15-item metadata: PKTOFF = 0)
            *item = 0;
        }
        let b = begin as i64 - 1;
        let to_addr = |v: i64| -> Result<usize> {
            let a = v + b;
            if a < 0 || a as usize > end {
                Err(bad("address in metadata outside the segment"))
            } else {
                Ok(a as usize)
            }
        };
        let gs = GenericSegment {
            const_base: to_addr(meta[0])?,
            n_const: meta[1].max(0) as usize,
            ref_base: to_addr(meta[5])?,
            n_ref: meta[6].max(0) as usize,
            ref_type: meta[4] as i32,
            pkt_dir_base: to_addr(meta[7])?,
            pkt_dir_type: meta[9] as i32,
            pkt_base: to_addr(meta[10])?,
            n_pkt: meta[11].max(0) as usize,
            pkt_size: meta[14].max(0) as usize,
            pkt_offset: meta[15].max(0) as usize,
            end,
        };
        // Constants, references and fixed-size packets must lie inside the segment.
        let within = |base: usize, n: usize| base.checked_add(n).is_some_and(|v| v <= end);
        if !within(gs.const_base, gs.n_const)
            || !within(gs.ref_base, gs.n_ref.max(2))
            || (gs.pkt_dir_type == 0
                && !(gs.pkt_size + gs.pkt_offset)
                    .checked_mul(gs.n_pkt)
                    .is_some_and(|v| within(gs.pkt_base, v)))
        {
            return Err(bad("metadata describes data outside the segment"));
        }
        Ok(gs)
    }

    /// Constant `k` (1-based).
    pub fn constant(&self, words: &[f64], k: usize) -> Option<f64> {
        if k == 0 || k > self.n_const {
            return None;
        }
        words.get(self.const_base + k - 1).copied()
    }

    /// Index (1-based) of the packet whose reference value applies to `x`, following SGFRVI.
    pub fn find_packet(&self, words: &[f64], x: f64) -> Option<usize> {
        match self.ref_type {
            0 | 1 => {
                let begin = *words.get(self.ref_base)?;
                let step = *words.get(self.ref_base + 1)?;
                let endref = begin + (self.n_pkt as f64 - 1.0) * step;
                if x < begin {
                    return if self.ref_type == 1 { Some(1) } else { None };
                }
                if x > endref {
                    return Some(self.n_pkt);
                }
                if self.n_pkt <= 1 {
                    return Some(1);
                }
                let t = (x - begin) / step;
                let idx = if self.ref_type == 0 { t.trunc() } else { t.round() } as usize + 1;
                Some(idx.min(self.n_pkt).max(1))
            }
            2..=4 => {
                let refs = words.get(self.ref_base..self.ref_base + self.n_ref)?;
                match self.ref_type {
                    2 => {
                        let idx = refs.partition_point(|&r| r < x);
                        if idx == 0 {
                            None
                        } else {
                            Some(idx)
                        }
                    }
                    3 => {
                        let idx = refs.partition_point(|&r| r <= x);
                        if idx == 0 {
                            None
                        } else {
                            Some(idx)
                        }
                    }
                    _ => {
                        let n = refs.len();
                        let idx = refs.partition_point(|&r| r <= x);
                        if idx == 0 {
                            Some(1)
                        } else if idx >= n {
                            Some(n)
                        } else if refs[idx] - x <= x - refs[idx - 1] {
                            Some(idx + 1)
                        } else {
                            Some(idx)
                        }
                    }
                }
            }
            _ => None,
        }
    }

    /// Contents of packet `i` (1-based).
    pub fn packet<'a>(&self, words: &'a [f64], i: usize) -> Option<&'a [f64]> {
        if i == 0 || i > self.n_pkt {
            return None;
        }
        if self.pkt_dir_type == 0 {
            let size = self.pkt_size + self.pkt_offset;
            let start = self.pkt_base + (i - 1) * size + self.pkt_offset; // 0-based index of first word
            words.get(start..start + self.pkt_size)
        } else {
            let b = self.pkt_dir_base + i - 1;
            let begin1 = super::spk::file_count_trunc(*words.get(b)?)?;
            let begin2 = super::spk::file_count_trunc(*words.get(b + 1)?)?;
            let size = begin2.checked_sub(begin1.checked_add(self.pkt_offset)?)?;
            let start = self.pkt_base.checked_add(begin1)?.checked_sub(1)?;
            if start.checked_add(size)? > self.end {
                return None;
            }
            words.get(start..start + size)
        }
    }
}
