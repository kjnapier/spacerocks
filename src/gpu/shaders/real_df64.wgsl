// Real numbers as double-single ("df64") pairs: the unevaluated sum x + y of two f32s with
// |y| <= ulp(x) / 2, about 48 bits of mantissa with f32's exponent range. The error-free
// transformations follow Dekker (1971), Knuth and the QD library (Hida, Li & Bailey 2001). The
// exact product uses fma where the driver fuses it (SPIR-V lets drivers round in between, and
// Mesa's software driver does) and Dekker's splitting otherwise: the host sets FMA_EXACT.
//
// The error-free transformations only work if each operation is rounded as written. Shader
// compilers don't promise that: SPIR-V lets drivers fuse a * b + c into an fma and Mesa's
// optimizer rewrites inexact float expressions, Metal compiles with fast math. So every
// intermediate they depend on goes through keep(), a bitwise AND with a mask that is all ones at
// run time but unknown to the compiler: an integer operation it can't see through, which pins
// the rounded value. The host checks these primitives against f64 on each adapter before it uses
// them (see `GpuContext::check`), and the Kepler test checks them inside a large kernel.
alias R = vec2<f32>;

fn keep(x: f32) -> f32 {
    return bitcast<f32>(bitcast<u32>(x) & MASK());
}

fn two_sum(a: f32, b: f32) -> vec2<f32> {
    let s = keep(a + b);
    let bb = keep(s - a);
    let e = (a - (s - bb)) + (b - bb);
    return vec2<f32>(s, e);
}

fn quick_two_sum(a: f32, b: f32) -> vec2<f32> {
    let s = keep(a + b);
    let e = b - (s - a);
    return vec2<f32>(s, e);
}

fn split(a: f32) -> vec2<f32> {
    let c = keep(4097.0 * a);
    let hi = keep(c - keep(c - a));
    return vec2<f32>(hi, a - hi);
}

fn two_prod(a: f32, b: f32) -> vec2<f32> {
    let p = keep(a * b);
    if (FMA_EXACT) {
        return vec2<f32>(p, fma(a, b, -p));
    }
    let sa = split(a);
    let sb = split(b);
    let e = ((sa.x * sb.x - p) + sa.x * sb.y + sa.y * sb.x) + sa.y * sb.y;
    return vec2<f32>(p, e);
}

fn rf(x: f32) -> R { return vec2<f32>(x, 0.0); }
fn rto(a: R) -> f32 { return a.x + a.y; }

fn radd(a: R, b: R) -> R {
    let s = two_sum(a.x, b.x);
    let t = two_sum(a.y, b.y);
    let u = quick_two_sum(s.x, s.y + t.x);
    return quick_two_sum(u.x, u.y + t.y);
}

fn rneg(a: R) -> R { return -a; }
fn rsub(a: R, b: R) -> R { return radd(a, -b); }

fn rmul(a: R, b: R) -> R {
    let p = two_prod(a.x, b.x);
    return quick_two_sum(p.x, p.y + fma(a.x, b.y, a.y * b.x));
}

fn rmulf(a: R, b: f32) -> R {
    let p = two_prod(a.x, b);
    return quick_two_sum(p.x, fma(a.y, b, p.y));
}

fn rdiv(a: R, b: R) -> R {
    // A quotient from the high parts, then two corrections from the remainder.
    let q1 = a.x / b.x;
    let r1 = rsub(a, rmulf(b, q1));
    let q2 = r1.x / b.x;
    let r2 = rsub(r1, rmulf(b, q2));
    let q3 = r2.x / b.x;
    let q = quick_two_sum(q1, q2);
    return radd(q, vec2<f32>(q3, 0.0));
}

fn rsqrt(a: R) -> R {
    if (a.x <= 0.0) {
        return vec2<f32>(0.0, 0.0);
    }
    // Karp's trick: one Newton correction of the f32 square root, in df64.
    let x = 1.0 / sqrt(a.x);
    let y = a.x * x;
    let yy = two_prod(y, y);
    let d = rsub(a, yy);
    return radd(vec2<f32>(y, 0.0), vec2<f32>(d.x * x * 0.5, 0.0));
}

fn rabs(a: R) -> R {
    if (a.x < 0.0) {
        return -a;
    }
    return a;
}

fn rlt(a: R, b: R) -> bool {
    return a.x < b.x || (a.x == b.x && a.y < b.y);
}

fn rfloor(a: R) -> R {
    let hi = floor(a.x);
    if (hi == a.x) {
        return quick_two_sum(hi, floor(a.y));
    }
    return vec2<f32>(hi, 0.0);
}

fn rtrunc(a: R) -> R {
    if (a.x < 0.0) {
        return -rfloor(-a);
    }
    return rfloor(a);
}

fn rround(a: R) -> R {
    return rfloor(radd(a, vec2<f32>(0.5, 0.0)));
}
