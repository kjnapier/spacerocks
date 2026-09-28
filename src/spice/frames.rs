//! Reference frames: built-in inertial frames, PCK body-fixed frames (binary and text PCK),
//! and text-kernel fixed-offset (TK) frames.
//!
//! All transforms are "state transforms" represented as a rotation `R` and its time
//! derivative `dR` (per second): for a state (r, v) in frame A, the state in frame B is
//! `(R r, dR r + R v)`.

use std::sync::OnceLock;

use nalgebra::Matrix3;

use super::builtin_frames::BUILTIN_FRAMES;
use super::error::{Result, SpiceError};
use super::kernel::SpiceKernel;
use super::math::{self, M3, IDENT, ZERO3};

/// Frame class numbers used by SPICE.
pub mod class {
    pub const INERTIAL: i32 = 1;
    pub const PCK: i32 = 2;
    pub const CK: i32 = 3;
    pub const TK: i32 = 4;
    pub const DYNAMIC: i32 = 5;
    pub const SWITCH: i32 = 6;
}

/// J2000 frame ID.
pub const J2000: i32 = 1;
/// ECLIPJ2000 frame ID.
pub const ECLIPJ2000: i32 = 17;
/// ITRF93 (high-precision Earth body-fixed) frame ID.
pub const ITRF93: i32 = 13000;
/// IAU_EARTH (the IAU rotation model of the Earth, from a text PCK such as pck00010.tpc) frame ID.
pub const IAU_EARTH: i32 = 10013;

/// Description of a reference frame.
#[derive(Debug, Clone, PartialEq)]
pub struct FrameInfo {
    pub id: i32,
    pub name: String,
    pub class: i32,
    pub class_id: i32,
    pub center: i32,
}

/// A state transformation: rotation and its derivative (per second).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StateTransform {
    pub rotation: M3,
    pub rate: M3,
}

impl StateTransform {
    pub const IDENTITY: StateTransform = StateTransform { rotation: IDENT, rate: ZERO3 };

    /// `self ∘ other`: first apply `other`, then `self`.
    #[inline]
    pub fn compose(&self, other: &StateTransform) -> StateTransform {
        StateTransform {
            rotation: math::mxm(&self.rotation, &other.rotation),
            rate: math::madd(
                &math::mxm(&self.rate, &other.rotation),
                &math::mxm(&self.rotation, &other.rate),
            ),
        }
    }

    #[inline]
    pub fn inverse(&self) -> StateTransform {
        StateTransform {
            rotation: math::transpose(&self.rotation),
            rate: math::transpose(&self.rate),
        }
    }

    #[inline]
    pub fn is_identity(&self) -> bool {
        self.rotation == IDENT && self.rate == ZERO3
    }

    /// Apply to a 6-vector state.
    #[inline]
    pub fn apply(&self, s: &[f64; 6]) -> [f64; 6] {
        let r = [s[0], s[1], s[2]];
        let v = [s[3], s[4], s[5]];
        let rr = math::mxv(&self.rotation, &r);
        let a = math::mxv(&self.rate, &r);
        let b = math::mxv(&self.rotation, &v);
        [rr[0], rr[1], rr[2], a[0] + b[0], a[1] + b[1], a[2] + b[2]]
    }

    pub fn rotation_matrix(&self) -> Matrix3<f64> {
        to_nalgebra(&self.rotation)
    }

    pub fn rate_matrix(&self) -> Matrix3<f64> {
        to_nalgebra(&self.rate)
    }
}

pub(crate) fn to_nalgebra(m: &M3) -> Matrix3<f64> {
    Matrix3::new(
        m[0][0], m[0][1], m[0][2], m[1][0], m[1][1], m[1][2], m[2][0], m[2][1], m[2][2],
    )
}

// ---------------------------------------------------------------------------------------------
// Built-in inertial frames (SPICELIB CHGIRF).
// ---------------------------------------------------------------------------------------------

const INERTIAL_NAMES: [&str; 21] = [
    "J2000", "B1950", "FK4", "DE-118", "DE-96", "DE-102", "DE-108", "DE-111", "DE-114", "DE-122", "DE-125",
    "DE-130", "GALACTIC", "DE-200", "DE-202", "MARSIAU", "ECLIPJ2000", "ECLIPB1950", "DE-140", "DE-142", "DE-143",
];

/// (base frame index (1-based), rotation sequence (arcsec, axis) applied as [a1]_x1 [a2]_x2 [a3]_x3)
const INERTIAL_DEFS: [(usize, &[(f64, u8)]); 21] = [
    (1, &[(0.0, 1)]),
    (1, &[(1152.84248596724, 3), (-1002.26108439117, 2), (1153.04066200330, 3)]),
    (2, &[(0.525, 3)]),
    (2, &[(0.53155, 3)]),
    (2, &[(0.4107, 3)]),
    (2, &[(0.1359, 3)]),
    (2, &[(0.4775, 3)]),
    (2, &[(0.5880, 3)]),
    (2, &[(0.5529, 3)]),
    (2, &[(0.5316, 3)]),
    (2, &[(0.5754, 3)]),
    (2, &[(0.5247, 3)]),
    (3, &[(1177200.0, 3), (225360.0, 1), (1016100.0, 3)]),
    (1, &[(0.0, 3)]),
    (1, &[(0.0, 3)]),
    (1, &[(324000.0, 3), (133610.4, 2), (-152348.4, 3)]),
    (1, &[(84381.448, 1)]),
    (2, &[(84404.836, 1)]),
    (1, &[(1152.71013777252, 3), (-1002.25042010533, 2), (1153.75719544491, 3)]),
    (1, &[(1152.72061453864, 3), (-1002.25052830351, 2), (1153.74663857521, 3)]),
    (1, &[(1153.03919093833, 3), (-1002.24822382286, 2), (1153.42900222357, 3)]),
];

/// Rotation matrices from J2000 to each built-in inertial frame (index = frame ID - 1).
fn inertial_from_j2000() -> &'static [M3; 21] {
    static TABLE: OnceLock<[M3; 21]> = OnceLock::new();
    TABLE.get_or_init(|| {
        let arcsec = std::f64::consts::PI / (180.0 * 3600.0);
        let mut out = [IDENT; 21];
        for i in 0..21 {
            let (base, seq) = INERTIAL_DEFS[i];
            let mut m = IDENT;
            for &(angle, axis) in seq.iter().rev() {
                m = math::mxm(&math::rotate(angle * arcsec, axis), &m);
            }
            // base index < i+1 always, so it is already computed
            let b = if base == i + 1 { IDENT } else { out[base - 1] };
            out[i] = math::mxm(&m, &b);
        }
        out
    })
}

/// Rotation taking vectors in built-in inertial frame `id` to J2000.
pub fn inertial_to_j2000(id: i32) -> Option<M3> {
    if (1..=21).contains(&id) {
        Some(math::transpose(&inertial_from_j2000()[(id - 1) as usize]))
    } else {
        None
    }
}

fn builtin_frames_by_id() -> &'static std::collections::HashMap<i32, FrameInfo> {
    static MAP: OnceLock<std::collections::HashMap<i32, FrameInfo>> = OnceLock::new();
    MAP.get_or_init(|| {
        let mut m = std::collections::HashMap::new();
        for (i, name) in INERTIAL_NAMES.iter().enumerate() {
            let id = i as i32 + 1;
            m.insert(
                id,
                FrameInfo {
                    id,
                    name: name.to_string(),
                    class: class::INERTIAL,
                    class_id: id,
                    center: 0,
                },
            );
        }
        for f in BUILTIN_FRAMES {
            m.insert(
                f.1,
                FrameInfo {
                    id: f.1,
                    name: f.0.to_string(),
                    class: f.3,
                    class_id: f.4,
                    center: f.2,
                },
            );
        }
        m
    })
}

fn builtin_frames_by_name() -> &'static std::collections::HashMap<String, i32> {
    static MAP: OnceLock<std::collections::HashMap<String, i32>> = OnceLock::new();
    MAP.get_or_init(|| builtin_frames_by_id().values().map(|f| (f.name.clone(), f.id)).collect())
}

fn builtin_frame_by_name(name: &str) -> Option<FrameInfo> {
    builtin_frames_by_name().get(name).and_then(|id| builtin_frame_by_id(*id))
}

fn builtin_frame_by_id(id: i32) -> Option<FrameInfo> {
    builtin_frames_by_id().get(&id).cloned()
}

const MAX_FRAME_DEPTH: usize = 32;

impl SpiceKernel {
    /// Collect frame definitions from the kernel pool (called whenever kernels change).
    pub(crate) fn rebuild_frames(&mut self) {
        self.kernel_frames.clear();
        self.kernel_frame_names.clear();
        let names: Vec<String> = self.pool.names().filter(|n| n.starts_with("FRAME_")).cloned().collect();
        for var in names {
            let rest = &var[6..];
            // FRAME_<id>_NAME defines a frame; FRAME_<name> = <id> maps a name.
            if let Some(id_str) = rest.strip_suffix("_NAME") {
                if let Ok(id) = id_str.parse::<i32>() {
                    let p = &self.pool;
                    let (Some(name), Some(cls)) = (p.get_str(&var), p.get_i32(&format!("FRAME_{}_CLASS", id))) else {
                        continue;
                    };
                    let class_id = p.get_i32(&format!("FRAME_{}_CLASS_ID", id)).unwrap_or(id);
                    let center = match p.get(&format!("FRAME_{}_CENTER", id)) {
                        Some(super::text::PoolValue::Numeric(v)) => v.first().map(|x| x.round() as i32).unwrap_or(0),
                        Some(super::text::PoolValue::Strings(v)) => v.first().and_then(|s| self.body_id(s)).unwrap_or(0),
                        None => 0,
                    };
                    self.kernel_frames.insert(
                        id,
                        FrameInfo {
                            id,
                            name: name.trim().to_uppercase(),
                            class: cls,
                            class_id,
                            center,
                        },
                    );
                    continue;
                }
            }
            if let Some(id) = self.pool.get_i32(&var) {
                let keyed_by_id = rest.split_once('_').map(|(a, _)| a.parse::<i32>().is_ok()).unwrap_or(false);
                if !keyed_by_id {
                    self.kernel_frame_names.insert(rest.to_string(), id);
                }
            }
        }
    }

    /// Frame ID for a frame name (kernel-defined frames take precedence over built-ins).
    pub fn frame_id(&self, name: &str) -> Option<i32> {
        let up = name.trim().to_uppercase();
        if let Some(&id) = self.kernel_frame_names.get(&up) {
            return Some(id);
        }
        if let Some(f) = builtin_frame_by_name(&up) {
            return Some(f.id);
        }
        up.parse::<i32>().ok().filter(|id| self.frame_info(*id).is_some())
    }

    /// Full description of a frame.
    pub fn frame_info(&self, id: i32) -> Option<FrameInfo> {
        if let Some(f) = self.kernel_frames.get(&id) {
            return Some(f.clone());
        }
        builtin_frame_by_id(id)
    }

    #[inline]
    fn with_frame_info<R>(&self, id: i32, f: impl FnOnce(&FrameInfo) -> R) -> Option<R> {
        if let Some(info) = self.kernel_frames.get(&id) {
            return Some(f(info));
        }
        builtin_frames_by_id().get(&id).map(f)
    }

    /// Resolve a frame given by name or ID string.
    pub fn frame(&self, name: &str) -> Result<FrameInfo> {
        self.frame_id(name)
            .and_then(|id| self.frame_info(id))
            .ok_or_else(|| SpiceError::UnknownFrame(name.to_string()))
    }

    /// Transform from frame `id` to its parent frame at `et`. Returns (parent, transform).
    fn frame_to_parent(&self, info: &FrameInfo, et: f64) -> Result<(i32, StateTransform)> {
        match info.class {
            class::INERTIAL => {
                let r = inertial_to_j2000(info.class_id).ok_or_else(|| SpiceError::UnsupportedFrameClass {
                    frame: info.name.clone(),
                    class: info.class,
                })?;
                Ok((J2000, StateTransform { rotation: r, rate: ZERO3 }))
            }
            class::PCK => {
                // Binary PCK data takes precedence over text PCK constants.
                if let Some((file, seg)) = self.find_pck(info.class_id, et) {
                    let f = &self.pck_files[file];
                    let s = &f.segments[seg];
                    let (m, dm) = s.rotation(f.daf.words(), et)?;
                    // m: reference -> body. We need body -> reference.
                    return Ok((
                        s.reference_frame,
                        StateTransform {
                            rotation: math::transpose(&m),
                            rate: math::transpose(&dm),
                        },
                    ));
                }
                self.text_pck_to_parent(info, et)
            }
            class::TK => self.tk_to_parent(info),
            other => Err(SpiceError::UnsupportedFrameClass {
                frame: info.name.clone(),
                class: other,
            }),
        }
    }

    /// State transform from frame `id` to J2000 at `et`.
    pub fn frame_to_j2000(&self, id: i32, et: f64) -> Result<StateTransform> {
        if id == J2000 {
            return Ok(StateTransform::IDENTITY);
        }
        if !self.kernel_frames.contains_key(&id) {
            if let Some(r) = inertial_to_j2000(id) {
                return Ok(StateTransform { rotation: r, rate: ZERO3 });
            }
        }
        let mut xf = StateTransform::IDENTITY;
        let mut cur = id;
        for _ in 0..MAX_FRAME_DEPTH {
            if cur == J2000 {
                return Ok(xf);
            }
            let (parent, step) = self
                .with_frame_info(cur, |info| self.frame_to_parent(info, et))
                .ok_or_else(|| SpiceError::UnknownFrame(cur.to_string()))??;
            xf = if xf.is_identity() { step } else { step.compose(&xf) };
            cur = parent;
        }
        Err(SpiceError::FrameChainTooDeep(id.to_string()))
    }

    /// State transformation matrix from frame `from` to frame `to` at `et`
    /// (TDB seconds past J2000). Equivalent to SPICE `sxform`.
    pub fn sxform(&self, from: &str, to: &str, et: f64) -> Result<StateTransform> {
        let f = self.frame(from)?.id;
        let t = self.frame(to)?.id;
        self.sxform_ids(f, t, et)
    }

    pub fn sxform_ids(&self, from: i32, to: i32, et: f64) -> Result<StateTransform> {
        if from == to {
            return Ok(StateTransform::IDENTITY);
        }
        let a = self.frame_to_j2000(from, et)?;
        let b = self.frame_to_j2000(to, et)?;
        Ok(b.inverse().compose(&a))
    }

    /// Rotation matrix from frame `from` to frame `to` at `et` (SPICE `pxform`).
    pub fn pxform(&self, from: &str, to: &str, et: f64) -> Result<Matrix3<f64>> {
        Ok(self.sxform(from, to, et)?.rotation_matrix())
    }

    // ---- TK frames -------------------------------------------------------------------------

    fn tk_var(&self, info: &FrameInfo, suffix: &str) -> Option<String> {
        let by_id = format!("TKFRAME_{}_{}", info.class_id, suffix);
        if self.pool.contains(&by_id) {
            return Some(by_id);
        }
        let by_name = format!("TKFRAME_{}_{}", info.name, suffix);
        if self.pool.contains(&by_name) {
            return Some(by_name);
        }
        None
    }

    fn tk_to_parent(&self, info: &FrameInfo) -> Result<(i32, StateTransform)> {
        let missing = |what: &str| SpiceError::PoolVariable {
            name: format!("TKFRAME_{}_{}", info.class_id, what),
            reason: format!("required to define TK frame {}", info.name),
        };
        let rel_var = self.tk_var(info, "RELATIVE").ok_or_else(|| missing("RELATIVE"))?;
        let rel_name = self.pool.get_str(&rel_var).ok_or_else(|| missing("RELATIVE"))?;
        let parent = self
            .frame_id(rel_name)
            .ok_or_else(|| SpiceError::UnknownFrame(rel_name.to_string()))?;
        let spec_var = self.tk_var(info, "SPEC").ok_or_else(|| missing("SPEC"))?;
        let spec = self.pool.get_str(&spec_var).ok_or_else(|| missing("SPEC"))?.trim().to_uppercase();
        let m: M3 = match spec.as_str() {
            "MATRIX" => {
                let v = self
                    .tk_var(info, "MATRIX")
                    .and_then(|n| self.pool.get_f64s(&n))
                    .filter(|v| v.len() == 9)
                    .ok_or_else(|| missing("MATRIX"))?;
                // Column-major in the kernel.
                let raw = [[v[0], v[3], v[6]], [v[1], v[4], v[7]], [v[2], v[5], v[8]]];
                sharpen(&raw)
            }
            "ANGLES" => {
                let a = self
                    .tk_var(info, "ANGLES")
                    .and_then(|n| self.pool.get_f64s(&n))
                    .filter(|v| v.len() == 3)
                    .ok_or_else(|| missing("ANGLES"))?;
                let ax = self
                    .tk_var(info, "AXES")
                    .and_then(|n| self.pool.get_f64s(&n))
                    .filter(|v| v.len() == 3)
                    .ok_or_else(|| missing("AXES"))?;
                let units = self
                    .tk_var(info, "UNITS")
                    .and_then(|n| self.pool.get_str(&n))
                    .ok_or_else(|| missing("UNITS"))?
                    .trim()
                    .to_uppercase();
                let f = angle_unit(&units).ok_or_else(|| SpiceError::PoolVariable {
                    name: format!("TKFRAME_{}_UNITS", info.class_id),
                    reason: format!("unrecognized angle unit '{}'", units),
                })?;
                if ax.iter().any(|&x| ![1.0, 2.0, 3.0].contains(&x.round())) {
                    return Err(SpiceError::PoolVariable {
                        name: format!("TKFRAME_{}_AXES", info.class_id),
                        reason: "axes must be 1, 2 or 3".into(),
                    });
                }
                let axis = |x: f64| x.round() as u8;
                math::mxm(
                    &math::rotate(a[0] * f, axis(ax[0])),
                    &math::mxm(&math::rotate(a[1] * f, axis(ax[1])), &math::rotate(a[2] * f, axis(ax[2]))),
                )
            }
            "QUATERNION" => {
                let q = self
                    .tk_var(info, "Q")
                    .and_then(|n| self.pool.get_f64s(&n))
                    .filter(|v| v.len() == 4)
                    .ok_or_else(|| missing("Q"))?;
                q2m(q)
            }
            other => {
                return Err(SpiceError::PoolVariable {
                    name: spec_var.clone(),
                    reason: format!("unrecognized TK frame specification '{}'", other),
                })
            }
        };
        Ok((parent, StateTransform { rotation: m, rate: ZERO3 }))
    }

    // ---- Text PCK (IAU rotation models) -----------------------------------------------------

    fn text_pck_to_parent(&self, info: &FrameInfo, et: f64) -> Result<(i32, StateTransform)> {
        let body = info.class_id;
        let p = &self.pool;
        let no_data = |detail: String| SpiceError::InsufficientOrientationData {
            frame: info.name.clone(),
            et,
            detail,
        };
        let ra = p
            .get_f64s(&format!("BODY{}_POLE_RA", body))
            .ok_or_else(|| no_data(format!("no binary PCK data covers this epoch and BODY{}_POLE_RA is not in the kernel pool", body)))?;
        let dec = p
            .get_f64s(&format!("BODY{}_POLE_DEC", body))
            .ok_or_else(|| no_data(format!("BODY{}_POLE_DEC not in kernel pool", body)))?;
        let pm = p
            .get_f64s(&format!("BODY{}_PM", body))
            .ok_or_else(|| no_data(format!("BODY{}_PM not in kernel pool", body)))?;

        // Reference frame and epoch of the constants are keyed on the body's barycenter
        // (SPICELIB ZZBODBRY): planets and satellites use their system barycenter.
        let refid = if body > 100 && body < 1000 { body / 100 } else { body };
        let parent = p
            .get_i32(&format!("BODY{}_CONSTANTS_REF_FRAME", refid))
            .or_else(|| p.get_i32(&format!("BODY{}_CONSTS_REF_FRAME", refid)))
            .unwrap_or(J2000);
        let jed = p
            .get_f64(&format!("BODY{}_CONSTANTS_JED_EPOCH", refid))
            .or_else(|| p.get_f64(&format!("BODY{}_CONSTS_JED_EPOCH", refid)))
            .unwrap_or(2_451_545.0);

        const D: f64 = 86_400.0;
        const T: f64 = 86_400.0 * 36_525.0;
        let epoch = et - D * (jed - 2_451_545.0);
        let td = epoch / D;
        let tc = epoch / T;
        let rpd = std::f64::consts::PI / 180.0;
        let c3 = |v: &[f64]| [v.first().copied().unwrap_or(0.0), v.get(1).copied().unwrap_or(0.0), v.get(2).copied().unwrap_or(0.0)];
        let (rc, dc, wc) = (c3(ra), c3(dec), c3(pm));

        // Degrees and degrees/second, evaluated exactly as TISBOD does.
        let mut ra_v = rc[0] + tc * (rc[1] + tc * rc[2]);
        let mut dec_v = dc[0] + tc * (dc[1] + tc * dc[2]);
        let mut w_v = wc[0] + td * (wc[1] + td * wc[2]);
        let mut ra_d = (rc[1] + 2.0 * tc * rc[2]) / T;
        let mut dec_d = (dc[1] + 2.0 * tc * dc[2]) / T;
        let mut w_d = (wc[1] + 2.0 * td * wc[2]) / D;

        let nut_ra = p.get_f64s(&format!("BODY{}_NUT_PREC_RA", body));
        let nut_dec = p.get_f64s(&format!("BODY{}_NUT_PREC_DEC", body));
        let nut_pm = p.get_f64s(&format!("BODY{}_NUT_PREC_PM", body));
        let angles_var = format!("BODY{}_NUT_PREC_ANGLES", refid);
        if let Some(angles) = p.get_f64s(&angles_var) {
            let nphsco = match p.get_f64(&format!("BODY{}_MAX_PHASE_DEGREE", refid)) {
                Some(deg) => {
                    let deg = deg.round() as usize;
                    if !(1..=3).contains(&deg) {
                        return Err(no_data(format!("BODY{}_MAX_PHASE_DEGREE must be 1, 2 or 3", refid)));
                    }
                    deg + 1
                }
                None => 2,
            };
            let nphase = angles.len() / nphsco;
            let longest = [nut_ra, nut_dec, nut_pm].iter().map(|v| v.map(|x| x.len()).unwrap_or(0)).max().unwrap_or(0);
            if longest > nphase {
                return Err(no_data(format!(
                    "{} has only {} phase angles but {} nutation/precession coefficients are given",
                    angles_var, nphase, longest
                )));
            }
            let (mut sra, mut sdec, mut sw, mut sdra, mut sddec, mut sdw) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
            for i in 0..nphase {
                let c = &angles[i * nphsco..(i + 1) * nphsco];
                let (theta, dtheta) = if nphsco == 2 {
                    ((c[0] + tc * c[1]) * rpd, c[1] / T * rpd)
                } else {
                    let mut th = 0.0;
                    for (j, cj) in c.iter().enumerate() {
                        th += tc.powi(j as i32) * cj;
                    }
                    // Matches TISBOD term by term (including its T^(l-1) divisor).
                    let mut dth = c[1] / T;
                    for (l, cl) in c.iter().enumerate().skip(2) {
                        dth += l as f64 * tc.powi(l as i32 - 1) * cl / T.powi(l as i32 - 1);
                    }
                    (th * rpd, dth * rpd)
                };
                let (s, co) = theta.sin_cos();
                let (ds, dc) = (co * dtheta, -s * dtheta);
                if let Some(v) = nut_ra.and_then(|v| v.get(i)) {
                    sra += v * s;
                    sdra += v * ds;
                }
                if let Some(v) = nut_dec.and_then(|v| v.get(i)) {
                    sdec += v * co;
                    sddec += v * dc;
                }
                if let Some(v) = nut_pm.and_then(|v| v.get(i)) {
                    sw += v * s;
                    sdw += v * ds;
                }
            }
            ra_v += sra;
            dec_v += sdec;
            w_v += sw;
            ra_d += sdra;
            dec_d += sddec;
            w_d += sdw;
        } else if nut_ra.is_some() || nut_dec.is_some() || nut_pm.is_some() {
            return Err(no_data(format!("{} not in kernel pool", angles_var)));
        }
        let ra_r = ra_v * rpd;
        let dec_r = dec_v * rpd;
        let w_r = math::fmod(w_v * rpd, std::f64::consts::TAU);
        let (m, dm) = math::euler_with_rate(
            [w_r, std::f64::consts::FRAC_PI_2 - dec_r, ra_r + std::f64::consts::FRAC_PI_2],
            [w_d * rpd, -dec_d * rpd, ra_d * rpd],
            [3, 1, 3],
        );
        // m: inertial -> body; we return body -> inertial.
        Ok((
            parent,
            StateTransform {
                rotation: math::transpose(&m),
                rate: math::transpose(&dm),
            },
        ))
    }
}

fn angle_unit(units: &str) -> Option<f64> {
    let deg = std::f64::consts::PI / 180.0;
    Some(match units {
        "RADIANS" => 1.0,
        "DEGREES" => deg,
        "ARCMINUTES" => deg / 60.0,
        "ARCSECONDS" => deg / 3600.0,
        "HOURANGLE" => 15.0 * deg,
        "MINUTEANGLE" => 15.0 * deg / 60.0,
        "SECONDANGLE" => 15.0 * deg / 3600.0,
        _ => return None,
    })
}

/// SPICE `Q2M` for a (possibly non-unit) SPICE-style quaternion (scalar first).
fn q2m(q: &[f64]) -> M3 {
    let n = (q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]).sqrt();
    let (q0, q1, q2, q3) = (q[0] / n, q[1] / n, q[2] / n, q[3] / n);
    let q01 = q0 * q1;
    let q02 = q0 * q2;
    let q03 = q0 * q3;
    let q11 = q1 * q1;
    let q12 = q1 * q2;
    let q13 = q1 * q3;
    let q22 = q2 * q2;
    let q23 = q2 * q3;
    let q33 = q3 * q3;
    [
        [1.0 - 2.0 * (q22 + q33), 2.0 * (q12 - q03), 2.0 * (q13 + q02)],
        [2.0 * (q12 + q03), 1.0 - 2.0 * (q11 + q33), 2.0 * (q23 - q01)],
        [2.0 * (q13 - q02), 2.0 * (q23 + q01), 1.0 - 2.0 * (q11 + q22)],
    ]
}

/// SPICE `SHARPR` (orthonormalize columns 1,2 then 3), with the sign fix-up TKFRAM applies.
fn sharpen(m: &M3) -> M3 {
    let col = |m: &M3, j: usize| [m[0][j], m[1][j], m[2][j]];
    let unit = |v: [f64; 3]| {
        let n = math::norm(&v);
        [v[0] / n, v[1] / n, v[2] / n]
    };
    let c1 = unit(col(m, 0));
    let mut c3 = unit(math::cross(&c1, &col(m, 1)));
    let mut c2 = unit(math::cross(&c3, &c1));
    if math::dot(&c2, &col(m, 1)) < 0.0 {
        c2 = [-c2[0], -c2[1], -c2[2]];
    }
    if math::dot(&c3, &col(m, 2)) < 0.0 {
        c3 = [-c3[0], -c3[1], -c3[2]];
    }
    [[c1[0], c2[0], c3[0]], [c1[1], c2[1], c3[1]], [c1[2], c2[2], c3[2]]]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn eclipj2000_matches_obliquity() {
        let m = inertial_from_j2000()[16];
        let eps = 84381.448f64 / 3600.0 * std::f64::consts::PI / 180.0;
        assert!((m[1][1] - eps.cos()).abs() < 1e-15);
        assert!((m[1][2] - eps.sin()).abs() < 1e-15);
    }
}
