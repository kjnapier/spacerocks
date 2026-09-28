//! Numerical building blocks shared by the SPK/PCK evaluators.
//!
//! Everything here is allocation-free; interpolation routines use fixed-size stack
//! buffers sized for the largest windows allowed by the SPICE file formats.

/// Largest number of points in any interpolation window we accept.
/// (Types 9/13/18/19 allow up to 27 points; we allow a little headroom.)
pub const MAX_WINDOW: usize = 32;

/// Largest Chebyshev expansion we accept (degree + 1).
pub const MAX_CHEB: usize = 64;

/// Three Chebyshev expansions of equal length (stored back to back in `coeffs`) evaluated
/// together, with derivatives. Per component the arithmetic is identical to
/// [`cheb_val_der`]; interleaving the three recurrences lets the CPU overlap them.
#[inline(always)]
pub fn cheb3_val_der(coeffs: &[f64], nc: usize, mid: f64, radius: f64, x: f64) -> ([f64; 3], [f64; 3]) {
    let (cx, rest) = coeffs.split_at(nc);
    let (cy, cz) = rest.split_at(nc);
    let cz = &cz[..nc];
    if nc == 0 {
        return ([0.0; 3], [0.0; 3]);
    }
    let s = (x - mid) / radius;
    let s2 = s * 2.0;
    let mut w0 = [0.0f64; 3];
    let mut w1 = [0.0f64; 3];
    let mut w2;
    let mut d0 = [0.0f64; 3];
    let mut d1 = [0.0f64; 3];
    let mut d2;
    let mut j = nc;
    while j > 1 {
        let c = [cx[j - 1], cy[j - 1], cz[j - 1]];
        w2 = w1;
        w1 = w0;
        d2 = d1;
        d1 = d0;
        for k in 0..3 {
            w0[k] = c[k] + (s2 * w1[k] - w2[k]);
            d0[k] = w1[k] * 2.0 + d1[k] * s2 - d2[k];
        }
        j -= 1;
    }
    let c0 = [cx[0], cy[0], cz[0]];
    let mut v = [0.0; 3];
    let mut dv = [0.0; 3];
    for k in 0..3 {
        v[k] = c0[k] + (s * w0[k] - w1[k]);
        dv[k] = (w0[k] + s * d0[k] - d1[k]) / radius;
    }
    (v, dv)
}

/// Three Chebyshev expansions evaluated together (values only; arithmetic as [`cheb_val`]).
#[inline(always)]
pub fn cheb3_val(coeffs: &[f64], nc: usize, mid: f64, radius: f64, x: f64) -> [f64; 3] {
    let (cx, rest) = coeffs.split_at(nc);
    let (cy, cz) = rest.split_at(nc);
    let cz = &cz[..nc];
    if nc == 0 {
        return [0.0; 3];
    }
    let s = (x - mid) / radius;
    let s2 = s * 2.0;
    let mut w0 = [0.0f64; 3];
    let mut w1 = [0.0f64; 3];
    let mut w2;
    let mut j = nc;
    while j > 1 {
        let c = [cx[j - 1], cy[j - 1], cz[j - 1]];
        w2 = w1;
        w1 = w0;
        for k in 0..3 {
            w0[k] = c[k] + (s2 * w1[k] - w2[k]);
        }
        j -= 1;
    }
    [
        s * w0[0] - w1[0] + cx[0],
        s * w0[1] - w1[1] + cy[0],
        s * w0[2] - w1[2] + cz[0],
    ]
}

/// Evaluate a Chebyshev expansion (value only; same operation order as SPICELIB `CHBVAL`).
#[inline(always)]
pub fn cheb_val(coeffs: &[f64], mid: f64, radius: f64, x: f64) -> f64 {
    let n = coeffs.len();
    if n == 0 {
        return 0.0;
    }
    let s = (x - mid) / radius;
    let s2 = s * 2.0;
    let (mut w0, mut w1, mut w2) = (0.0f64, 0.0f64, 0.0f64);
    let mut j = n;
    while j > 1 {
        w2 = w1;
        w1 = w0;
        w0 = coeffs[j - 1] + (s2 * w1 - w2);
        j -= 1;
    }
    let _ = w2;
    s * w0 - w1 + coeffs[0]
}

/// Value of a Chebyshev expansion and its integral from `mid` to `x`
/// (the integral is in units of `x`). Direct port of SPICELIB `CHBIGR`.
pub fn cheb_val_integral(cp: &[f64], mid: f64, radius: f64, x: f64) -> (f64, f64) {
    let nterms = cp.len();
    let degp = nterms as i64 - 1;
    let s = (x - mid) / radius;
    let s2 = 2.0 * s;

    let a2 = if nterms >= 3 { cp[0] - 0.5 * cp[2] } else { cp[0] };
    let mut adegp1 = 0.0;
    let mut adegp2 = 0.0;
    if degp >= 2 {
        adegp1 = 0.5 * cp[(degp - 1) as usize] / degp as f64;
    }
    if degp >= 1 {
        adegp2 = 0.5 * cp[degp as usize] / (degp + 1) as f64;
    }
    let mut f = [0.0f64; 3];
    f[0] = if degp == 0 { a2 } else { adegp2 };
    let mut z = [f[0], 0.0, 0.0];
    let mut w = [0.0f64; 3];

    let mut i = nterms;
    while i > 1 {
        let ai = if i == 2 {
            a2
        } else if i < nterms {
            0.5 * (cp[i - 2] - cp[i]) / (i - 1) as f64
        } else {
            adegp1
        };
        f[2] = f[1];
        f[1] = f[0];
        f[0] = ai + (s2 * f[1] - f[2]);
        z[2] = z[1];
        z[1] = z[0];
        z[0] = ai - z[2];
        w[2] = w[1];
        w[1] = w[0];
        w[0] = cp[i - 1] + (s2 * w[1] - w[2]);
        i -= 1;
    }
    let c0 = z[1];
    let integral = radius * (c0 + s * f[0] - f[1]);
    let value = cp[0] + (s * w[0] - w[1]);
    (value, integral)
}

/// Fortran `MOD` for doubles as implemented by f2c (`x - y*trunc(x/y)`), used where SPICE
/// uses it so results round identically.
#[inline(always)]
pub fn fmod(x: f64, y: f64) -> f64 {
    let q = x / y;
    let q = if q >= 0.0 { q.floor() } else { -(-q).floor() };
    x - y * q
}

/// Lagrange interpolation (Neville's algorithm) of `ys` given at abscissas `xs`, evaluated at `x`.
#[inline]
pub fn lagrange(xs: &[f64], ys: &[f64], x: f64) -> f64 {
    let n = xs.len();
    debug_assert!(n <= MAX_WINDOW && ys.len() == n);
    let mut p = [0.0f64; MAX_WINDOW];
    p[..n].copy_from_slice(&ys[..n]);
    for j in 1..n {
        for i in 0..(n - j) {
            let denom = xs[i] - xs[i + j];
            p[i] = ((x - xs[i + j]) * p[i] + (xs[i] - x) * p[i + 1]) / denom;
        }
    }
    p[0]
}

/// Hermite interpolation using function values `f` and derivatives `df` at abscissas `xs`.
/// Returns `(value, derivative)` at `x`.
///
/// This is a direct port of SPICELIB `HRMINT` (a Neville-style scheme), so results agree with
/// CSPICE to the last bit or two even for poorly spaced abscissas.
#[inline]
#[allow(clippy::assign_op_pattern, clippy::manual_div_ceil)]
pub fn hermite(xs: &[f64], f: &[f64], df: &[f64], x: f64) -> (f64, f64) {
    let n = xs.len();
    debug_assert!((1..=MAX_WINDOW).contains(&n));
    // 1-based arrays of length 2n as in the Fortran; index 0 unused.
    let mut w1 = [0.0f64; 2 * MAX_WINDOW + 1];
    let mut w2 = [0.0f64; 2 * MAX_WINDOW + 1];
    for i in 0..n {
        w1[2 * i + 1] = f[i];
        w1[2 * i + 2] = df[i];
    }
    let xv = |k: usize| xs[k - 1];
    for i in 1..n {
        let c1 = xv(i + 1) - x;
        let c2 = x - xv(i);
        let denom = xv(i + 1) - xv(i);
        let prev = 2 * i - 1;
        let this = prev + 1;
        let next = this + 1;
        w2[prev] = w1[this];
        w2[this] = (w1[next] - w1[prev]) / denom;
        let temp = w1[this] * (x - xv(i)) + w1[prev];
        w1[this] = (c1 * w1[prev] + c2 * w1[next]) / denom;
        w1[prev] = temp;
    }
    w2[2 * n - 1] = w1[2 * n];
    w1[2 * n - 1] = w1[2 * n] * (x - xv(n)) + w1[2 * n - 1];
    for j in 2..=(2 * n - 1) {
        for i in 1..=(2 * n - j) {
            let xi = (i + 1) / 2;
            let xij = (i + j + 1) / 2;
            let c1 = xv(xij) - x;
            let c2 = x - xv(xi);
            let denom = xv(xij) - xv(xi);
            w2[i] = (c1 * w2[i] + c2 * w2[i + 1] + (w1[i + 1] - w1[i])) / denom;
            w1[i] = (c1 * w1[i] + c2 * w1[i + 1]) / denom;
        }
    }
    (w1[1], w2[1])
}

/// Stumpff functions c0..c3 as computed by SPICELIB `STMP03`.
pub fn stumpff(x: f64) -> (f64, f64, f64, f64) {
    // lbound = -(ln 2 + ln dpmax)^2
    const Y: f64 = std::f64::consts::LN_2 + 709.782_712_893_384; // ln 2 + ln(DPMAX)
    const LBOUND: f64 = -Y * Y;
    if x <= LBOUND {
        return (f64::NAN, f64::NAN, f64::NAN, f64::NAN);
    }
    if x < -1.0 {
        let z = (-x).sqrt();
        let c0 = z.cosh();
        let c1 = z.sinh() / z;
        return (c0, c1, (1.0 - c0) / x, (1.0 - c1) / x);
    }
    if x > 1.0 {
        let z = x.sqrt();
        let c0 = z.cos();
        let c1 = z.sin() / z;
        return (c0, c1, (1.0 - c0) / x, (1.0 - c1) / x);
    }
    let pairs = |i: usize| 1.0 / ((i as f64) * ((i + 1) as f64));
    let mut c3 = 1.0;
    let mut i = 20;
    while i >= 4 {
        c3 = 1.0 - x * pairs(i) * c3;
        i -= 2;
    }
    c3 *= pairs(2);
    let mut c2 = 1.0;
    let mut i = 19;
    while i >= 3 {
        c2 = 1.0 - x * pairs(i) * c2;
        i -= 2;
    }
    c2 *= pairs(1);
    let c1 = 1.0 - x * c3;
    let c0 = 1.0 - x * c2;
    (c0, c1, c2, c3)
}

/// Two-body propagation (port of SPICELIB `PROP2B`) of state `pv` by `dt` under gravitational
/// parameter `gm`. Units must be consistent (km, km/s, s for SPK use).
pub fn prop2b(gm: f64, pv: &[f64; 6], dt: f64) -> Option<[f64; 6]> {
    let pos = [pv[0], pv[1], pv[2]];
    let vel = [pv[3], pv[4], pv[5]];
    if gm <= 0.0 || pos == [0.0; 3] || vel == [0.0; 3] {
        return None;
    }
    if dt == 0.0 {
        return Some(*pv);
    }
    let r0 = norm(&pos);
    let rv = dot(&pos, &vel);
    let hvec = cross(&pos, &vel);
    let h2 = dot(&hvec, &hvec);
    if h2 == 0.0 {
        return None;
    }
    let tmp = cross(&vel, &hvec);
    let eqvec = [
        tmp[0] / gm - pos[0] / r0,
        tmp[1] / gm - pos[1] / r0,
        tmp[2] / gm - pos[2] / r0,
    ];
    let e = norm(&eqvec);
    let q = h2 / (gm * (1.0 + e));
    let f = 1.0 - e;
    let b = (q / gm).sqrt();
    let br0 = b * r0;
    let b2rv = b * b * rv;
    let bq = b * q;
    let qovr0 = q / r0;
    let maxc = 1.0f64
        .max(br0.abs())
        .max(b2rv.abs())
        .max(bq.abs())
        .max((qovr0 / bq).abs());
    let bound = if f < 0.0 {
        let logmxc = maxc.ln();
        let logdpm = (f64::MAX / 2.0).ln();
        let fixed = logdpm - logmxc;
        let rootf = (-f).sqrt();
        let logf = (-f).ln();
        (fixed / rootf).min((fixed + 1.5 * logf) / rootf)
    } else {
        let logbnd = (1.5f64.ln() + f64::MAX.ln() - maxc.ln()) / 3.0;
        logbnd.exp()
    };

    // PROP2B evaluates Kepler's equation with two slightly different operation orders in
    // different places; both are mirrored so the bisection follows the same path.
    if !(bound > 0.0) {
        return None;
    }
    let kfun = |x: f64| {
        let (_, c1, c2, c3) = stumpff(f * x * x);
        x * (br0 * c1 + x * (b2rv * c2 + x * (bq * c3)))
    };
    let kfun2 = |x: f64| {
        let (_, c1, c2, c3) = stumpff(f * x * x);
        x * (br0 * c1 + x * (b2rv * c2 + x * bq * c3))
    };

    let mut x = (dt / bq).clamp(-bound, bound);
    let mut k = kfun(x);
    let (mut lower, mut upper);
    if dt < 0.0 {
        upper = 0.0;
        lower = x;
        while k > dt {
            upper = lower;
            lower *= 2.0;
            let oldx = x;
            x = lower.clamp(-bound, bound);
            if x == oldx {
                return None;
            }
            k = kfun(x);
        }
    } else {
        lower = 0.0;
        upper = x;
        while k < dt {
            lower = upper;
            upper *= 2.0;
            let oldx = x;
            x = upper.clamp(-bound, bound);
            if x == oldx {
                return None;
            }
            k = kfun2(x);
        }
    }

    x = upper.min(lower.max((lower + upper) / 2.0));
    let mut lcount = 0;
    let mut mostc = 1000;
    while x > lower && x < upper && lcount < mostc {
        let k = kfun2(x);
        if k > dt {
            upper = x;
        } else if k < dt {
            lower = x;
        } else {
            upper = x;
            lower = x;
        }
        if mostc > 64 && upper != 0.0 && lower != 0.0 {
            mostc = 64;
            lcount = 0;
        }
        x = upper.min(lower.max((lower + upper) / 2.0));
        lcount += 1;
    }

    let (c0, c1, c2, c3) = stumpff(f * x * x);
    let x2 = x * x;
    let x3 = x2 * x;
    let br = br0 * c0 + x * (b2rv * c1 + x * (bq * c2));
    let pc = 1.0 - qovr0 * x2 * c2;
    let vc = dt - bq * x3 * c3;
    let pcdot = -(qovr0 / br) * x * c1;
    let vcdot = 1.0 - bq / br * x2 * c2;
    Some([
        pc * pos[0] + vc * vel[0],
        pc * pos[1] + vc * vel[1],
        pc * pos[2] + vc * vel[2],
        pcdot * pos[0] + vcdot * vel[0],
        pcdot * pos[1] + vcdot * vel[1],
        pcdot * pos[2] + vcdot * vel[2],
    ])
}

/// Solve the equinoctial form of Kepler's equation `ml = f + h cos f - k sin f` for `f`.
pub fn kepleq(ml: f64, h: f64, k: f64) -> f64 {
    let two_pi = std::f64::consts::TAU;
    let ml = ml.rem_euclid(two_pi);
    let mut f = ml;
    for _ in 0..100 {
        let (s, c) = f.sin_cos();
        let func = f + h * c - k * s - ml;
        let dfunc = 1.0 - h * s - k * c;
        let delta = func / dfunc;
        f -= delta;
        if delta.abs() <= 1e-15 * (1.0 + f.abs()) {
            break;
        }
    }
    f
}

// ---------------------------------------------------------------------------------------------
// Small 3-vector / 3x3 matrix helpers (row-major [[f64;3];3]).
// ---------------------------------------------------------------------------------------------

pub type M3 = [[f64; 3]; 3];

pub const IDENT: M3 = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
pub const ZERO3: M3 = [[0.0; 3]; 3];

#[inline(always)]
pub fn dot(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline(always)]
pub fn norm(a: &[f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

#[inline(always)]
pub fn cross(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[inline(always)]
pub fn mxv(m: &M3, v: &[f64; 3]) -> [f64; 3] {
    [
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
    ]
}

#[inline(always)]
pub fn mxm(a: &M3, b: &M3) -> M3 {
    let mut out = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            out[i][j] = a[i][0] * b[0][j] + a[i][1] * b[1][j] + a[i][2] * b[2][j];
        }
    }
    out
}

#[inline(always)]
pub fn madd(a: &M3, b: &M3) -> M3 {
    let mut out = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            out[i][j] = a[i][j] + b[i][j];
        }
    }
    out
}

#[inline(always)]
pub fn transpose(a: &M3) -> M3 {
    [
        [a[0][0], a[1][0], a[2][0]],
        [a[0][1], a[1][1], a[2][1]],
        [a[0][2], a[1][2], a[2][2]],
    ]
}

/// SPICE `ROTATE`: matrix that rotates the *frame* by `angle` about `axis` (1, 2 or 3).
#[inline]
pub fn rotate(angle: f64, axis: u8) -> M3 {
    let (s, c) = angle.sin_cos();
    match axis {
        1 => [[1.0, 0.0, 0.0], [0.0, c, s], [0.0, -s, c]],
        2 => [[c, 0.0, -s], [0.0, 1.0, 0.0], [s, 0.0, c]],
        _ => [[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]],
    }
}

/// Derivative of `rotate(angle, axis)` with respect to `angle`.
#[inline]
pub fn drotate(angle: f64, axis: u8) -> M3 {
    let (s, c) = angle.sin_cos();
    match axis {
        1 => [[0.0, 0.0, 0.0], [0.0, -s, c], [0.0, -c, -s]],
        2 => [[-s, 0.0, -c], [0.0, 0.0, 0.0], [c, 0.0, -s]],
        _ => [[-s, c, 0.0], [-c, -s, 0.0], [0.0, 0.0, 0.0]],
    }
}

/// Rotation `[a1]_ax1 [a2]_ax2 [a3]_ax3` and its time derivative, given angle rates.
pub fn euler_with_rate(
    angles: [f64; 3],
    rates: [f64; 3],
    axes: [u8; 3],
) -> (M3, M3) {
    let r1 = rotate(angles[0], axes[0]);
    let r2 = rotate(angles[1], axes[1]);
    let r3 = rotate(angles[2], axes[2]);
    let d1 = scale(&drotate(angles[0], axes[0]), rates[0]);
    let d2 = scale(&drotate(angles[1], axes[1]), rates[1]);
    let d3 = scale(&drotate(angles[2], axes[2]), rates[2]);
    let r23 = mxm(&r2, &r3);
    let r = mxm(&r1, &r23);
    let dr = madd(
        &madd(&mxm(&d1, &r23), &mxm(&r1, &mxm(&d2, &r3))),
        &mxm(&r1, &mxm(&r2, &d3)),
    );
    (r, dr)
}

#[inline(always)]
pub fn scale(a: &M3, s: f64) -> M3 {
    let mut out = *a;
    for row in out.iter_mut() {
        for v in row.iter_mut() {
            *v *= s;
        }
    }
    out
}

/// Rotate vector `v` about axis `axis` by `theta` (SPICE `VROTV`).
pub fn vrotv(v: &[f64; 3], axis: &[f64; 3], theta: f64) -> [f64; 3] {
    let n = norm(axis);
    if n == 0.0 {
        return *v;
    }
    let x = [axis[0] / n, axis[1] / n, axis[2] / n];
    let p = dot(v, &x);
    let proj = [p * x[0], p * x[1], p * x[2]];
    let v1 = [v[0] - proj[0], v[1] - proj[1], v[2] - proj[2]];
    let v2 = cross(&x, &v1);
    let (s, c) = theta.sin_cos();
    [
        proj[0] + c * v1[0] + s * v2[0],
        proj[1] + c * v1[1] + s * v2[1],
        proj[2] + c * v1[2] + s * v2[2],
    ]
}

/// Angular separation of two vectors (SPICE `VSEP`), numerically robust.
pub fn vsep(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    let na = norm(a);
    let nb = norm(b);
    if na == 0.0 || nb == 0.0 {
        return 0.0;
    }
    let ua = [a[0] / na, a[1] / na, a[2] / na];
    let ub = [b[0] / nb, b[1] / nb, b[2] / nb];
    if dot(&ua, &ub) > 0.0 {
        let d = [ua[0] - ub[0], ua[1] - ub[1], ua[2] - ub[2]];
        2.0 * (0.5 * norm(&d)).asin()
    } else if dot(&ua, &ub) < 0.0 {
        let s = [ua[0] + ub[0], ua[1] + ub[1], ua[2] + ub[2]];
        std::f64::consts::PI - 2.0 * (0.5 * norm(&s)).asin()
    } else {
        std::f64::consts::FRAC_PI_2
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn chebyshev_matches_direct() {
        let c = [1.0, -0.5, 0.25, 0.125, -0.3];
        let x: f64 = 0.37;
        let t = [1.0, x, 2.0 * x * x - 1.0, 4.0 * x * x * x - 3.0 * x, 8.0 * x.powi(4) - 8.0 * x * x + 1.0];
        let dt = [0.0, 1.0, 4.0 * x, 12.0 * x * x - 3.0, 32.0 * x.powi(3) - 16.0 * x];
        let v: f64 = c.iter().zip(t.iter()).map(|(a, b)| a * b).sum();
        let d: f64 = c.iter().zip(dt.iter()).map(|(a, b)| a * b).sum();
        let mut three = Vec::new();
        for _ in 0..3 {
            three.extend_from_slice(&c);
        }
        let (vv, dd) = cheb3_val_der(&three, c.len(), 0.0, 1.0, x);
        assert!((v - vv[2]).abs() < 1e-15 && (d - dd[2]).abs() < 1e-14);
        assert!((cheb_val(&c, 0.0, 1.0, x) - v).abs() < 1e-15);
    }

    #[test]
    fn hermite_reproduces_cubic() {
        let p = |x: f64| 2.0 * x.powi(3) - x + 3.0;
        let dp = |x: f64| 6.0 * x * x - 1.0;
        let xs = [0.0, 1.5];
        let (v, d) = hermite(&xs, &[p(0.0), p(1.5)], &[dp(0.0), dp(1.5)], 0.7);
        assert!((v - p(0.7)).abs() < 1e-13 && (d - dp(0.7)).abs() < 1e-13);
    }

    #[test]
    fn lagrange_reproduces_quadratic() {
        let xs = [0.0, 1.0, 3.0];
        let ys: Vec<f64> = xs.iter().map(|x| x * x - 2.0 * x).collect();
        assert!((lagrange(&xs, &ys, 2.0) - 0.0).abs() < 1e-14);
    }

    #[test]
    fn prop2b_full_orbit() {
        let gm = 398600.4418;
        let r: f64 = 7000.0;
        let v = (gm / r).sqrt();
        let pv = [r, 0.0, 0.0, 0.0, v, 0.0];
        let period = 2.0 * std::f64::consts::PI * (r.powi(3) / gm).sqrt();
        let out = prop2b(gm, &pv, period).unwrap();
        assert!((out[0] - r).abs() < 1e-6 && out[1].abs() < 1e-6);
        let back = prop2b(gm, &pv, -period / 4.0).unwrap();
        assert!((back[1] + r).abs() < 1e-6);
    }
}
