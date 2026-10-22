// IAS15 (Rein & Spiegel 2015) for test particles in the field of massive bodies whose trajectory
// is given by their own IAS15 steps: each step's start (x0, v0, a0) and Gauss–Radau coefficients
// b0..b6, from which the bodies' positions follow anywhere in the step (as in
// `IAS15::interpolate`). One invocation per test particle, each with its own adaptive timestep,
// all running from time 0 to T_END (coef[0]) of the block; a dispatch takes at most
// P.max_steps steps and the host dispatches again until every particle has arrived.
//
// The step is the CPU IAS15's (nbody/integrators/ias15.rs) for a single particle: the
// predictor-corrector error and the timestep come from this particle's own coefficients. Needs
// the constants H, RRI, C, D, XF, VF, the series factors and SAFETY, EPS_FAC, PC_TOL (f32).

struct Params {
    n: u32,
    mask: u32,
    // Rows of the bodies' table.
    nrows: u32,
    // Massive bodies.
    nb: u32,
    // Most steps per particle in this dispatch.
    max_steps: u32,
    c: u32,
    d: u32,
    e: u32,
}

@group(0) @binding(0) var<uniform> P: Params;
// Per particle: x, v, t, dt, dt_last, then b_last[0..7] and e_last[0..7] (vectors).
@group(0) @binding(1) var<storage, read_write> st: array<R>;
// Per row: t0, 1 / dt, then per body x0, v0, a0, b0..b6.
@group(0) @binding(2) var<storage, read> tab: array<R>;
// T_END, then G m of each body.
@group(0) @binding(3) var<storage, read> coef: array<R>;
// Per particle: status (0 running, 1 arrived, 2 failed) and the table row it is in.
@group(0) @binding(4) var<storage, read_write> status: array<u32>;
// Particles still running after this dispatch, and steps taken in it (accepted ones).
@group(0) @binding(5) var<storage, read_write> counts: array<atomic<u32>, 2>;

const ST: u32 = 51u;

fn MASK() -> u32 { return P.mask; }

fn rowbase(k: u32) -> u32 {
    return k * (2u + 30u * P.nb);
}

fn tabv(i: u32) -> V3 {
    return V3(tab[i], tab[i + 1u], tab[i + 2u]);
}

fn stv(i: u32) -> V3 {
    return V3(st[i], st[i + 1u], st[i + 2u]);
}

fn put(i: u32, v: V3) {
    st[i] = v.x;
    st[i + 1u] = v.y;
    st[i + 2u] = v.z;
}

// |v| in f32 without overflow or underflow of the squares.
fn norm32(v: vec3<f32>) -> f32 {
    let m = max(max(abs(v.x), abs(v.y)), abs(v.z));
    if (m == 0.0) {
        return 0.0;
    }
    let w = v / m;
    return m * sqrt(dot(w, w));
}

// Acceleration of a test particle at x at time tau from the bodies; `cur` is the table row of
// the last call, moved to the row holding tau.
fn accel(x: V3, tau: R, cur: ptr<function, u32>) -> V3 {
    loop {
        if (*cur + 1u < P.nrows && !rlt(tau, tab[rowbase(*cur + 1u)])) {
            *cur = *cur + 1u;
        } else {
            break;
        }
    }
    loop {
        if (*cur > 0u && rlt(tau, tab[rowbase(*cur)])) {
            *cur = *cur - 1u;
        } else {
            break;
        }
    }
    let base = rowbase(*cur);
    let s0 = rsub(tau, tab[base]);
    let h = rmul(s0, tab[base + 1u]);
    let s1 = rmulf(rmul(s0, s0), 0.5);
    let s2 = rmul(rmul(s1, h), THIRD);
    let s3 = rmulf(rmul(s2, h), 0.5);
    let s4 = rmul(rmul(s3, h), THREE_FIFTHS);
    let s5 = rmul(rmul(s4, h), TWO_THIRDS);
    let s6 = rmul(rmul(s5, h), FIVE_SEVENTHS);
    let s7 = rmulf(rmul(s6, h), 0.75);
    let s8 = rmul(rmul(s7, h), SEVEN_NINTHS);
    var acc = vzero();
    for (var j = 0u; j < P.nb; j = j + 1u) {
        let o = base + 2u + 30u * j;
        var xb = vscale(s8, tabv(o + 27u));
        xb = vaxpy(xb, s7, tabv(o + 24u));
        xb = vaxpy(xb, s6, tabv(o + 21u));
        xb = vaxpy(xb, s5, tabv(o + 18u));
        xb = vaxpy(xb, s4, tabv(o + 15u));
        xb = vaxpy(xb, s3, tabv(o + 12u));
        xb = vaxpy(xb, s2, tabv(o + 9u));
        xb = vaxpy(xb, s1, tabv(o + 6u));
        xb = vaxpy(xb, s0, tabv(o + 3u));
        xb = vadd(xb, tabv(o));
        let d = vsub(xb, x);
        let r2 = vdot(d, d);
        acc = vaxpy(acc, rdiv(coef[1u + j], rmul(r2, rsqrt(r2))), d);
    }
    return acc;
}

// Coefficients for a step `ratio` times the last one, from that step's (as
// `predict_next_coefficients`).
fn predict(ratio: R, bl: ptr<function, array<V3, 7>>, el: ptr<function, array<V3, 7>>, b: ptr<function, array<V3, 7>>, e: ptr<function, array<V3, 7>>) {
    if (abs(rto(ratio)) > 20.0) {
        for (var k = 0u; k < 7u; k = k + 1u) {
            (*b)[k] = vzero();
            (*e)[k] = vzero();
        }
        return;
    }
    let q1 = ratio;
    let q2 = rmul(q1, q1);
    let q3 = rmul(q1, q2);
    let q4 = rmul(q2, q2);
    let q5 = rmul(q2, q3);
    let q6 = rmul(q3, q3);
    let q7 = rmul(q3, q4);
    let p = *bl;
    // Binomial weights: e_k = q^(k+1) sum_j C(j+1, k+1) b_j.
    let e0 = vscale(q1, vadd(vadd(vadd(vscale(rf(7.0), p[6]), vscale(rf(6.0), p[5])), vadd(vscale(rf(5.0), p[4]), vscale(rf(4.0), p[3]))), vadd(vadd(vscale(rf(3.0), p[2]), vscale(rf(2.0), p[1])), p[0])));
    let e1 = vscale(q2, vadd(vadd(vadd(vscale(rf(21.0), p[6]), vscale(rf(15.0), p[5])), vadd(vscale(rf(10.0), p[4]), vscale(rf(6.0), p[3]))), vadd(vscale(rf(3.0), p[2]), p[1])));
    let e2 = vscale(q3, vadd(vadd(vscale(rf(35.0), p[6]), vscale(rf(20.0), p[5])), vadd(vadd(vscale(rf(10.0), p[4]), vscale(rf(4.0), p[3])), p[2])));
    let e3 = vscale(q4, vadd(vadd(vscale(rf(35.0), p[6]), vscale(rf(15.0), p[5])), vadd(vscale(rf(5.0), p[4]), p[3])));
    let e4 = vscale(q5, vadd(vadd(vscale(rf(21.0), p[6]), vscale(rf(6.0), p[5])), p[4]));
    let e5 = vscale(q6, vadd(vscale(rf(7.0), p[6]), p[5]));
    let e6 = vscale(q7, p[6]);
    let en = array<V3, 7>(e0, e1, e2, e3, e4, e5, e6);
    for (var k = 0u; k < 7u; k = k + 1u) {
        (*b)[k] = vadd(en[k], vsub(p[k], (*el)[k]));
        (*e)[k] = en[k];
    }
}

@compute @workgroup_size(64)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    if (i >= P.n || status[2u * i] != 0u) {
        return;
    }
    let o = ST * i;
    var x0 = stv(o);
    var v0 = stv(o + 3u);
    var t = st[o + 6u];
    var dt = st[o + 7u];
    var dt_last = st[o + 8u];
    var bl: array<V3, 7>;
    var el: array<V3, 7>;
    for (var k = 0u; k < 7u; k = k + 1u) {
        bl[k] = stv(o + 9u + 3u * k);
        el[k] = stv(o + 30u + 3u * k);
    }
    var cur = status[2u * i + 1u];
    let t_end = coef[0];
    var state = 0u;

    var b: array<V3, 7>;
    var e: array<V3, 7>;
    var g: array<V3, 7>;
    var taken = 0u;
    for (var steps = 0u; steps < P.max_steps; steps = steps + 1u) {
        if (!rlt(t, t_end)) {
            state = 1u;
            break;
        }
        // Land on the end of the block, then go back to the natural timestep.
        let natural = dt;
        var clipped = false;
        if (rlt(rsub(t_end, t), dt)) {
            dt = rsub(t_end, t);
            clipped = true;
        }
        let a0 = accel(x0, t, &cur);
        var rejections = 0u;
        var accepted = false;
        loop {
            // Predicted coefficients for this step.
            if (rto(dt_last) != 0.0) {
                predict(rdiv(dt, dt_last), &bl, &el, &b, &e);
            } else {
                for (var k = 0u; k < 7u; k = k + 1u) {
                    b[k] = vzero();
                    e[k] = vzero();
                }
            }
            g[0] = vadd(vadd(vadd(vscale(D[15], b[6]), vscale(D[10], b[5])), vadd(vscale(D[6], b[4]), vscale(D[3], b[3]))), vadd(vadd(vscale(D[1], b[2]), vscale(D[0], b[1])), b[0]));
            g[1] = vadd(vadd(vadd(vscale(D[16], b[6]), vscale(D[11], b[5])), vadd(vscale(D[7], b[4]), vscale(D[4], b[3]))), vadd(vscale(D[2], b[2]), b[1]));
            g[2] = vadd(vadd(vscale(D[17], b[6]), vscale(D[12], b[5])), vadd(vadd(vscale(D[8], b[4]), vscale(D[5], b[3])), b[2]));
            g[3] = vadd(vadd(vscale(D[18], b[6]), vscale(D[13], b[5])), vadd(vscale(D[9], b[4]), b[3]));
            g[4] = vadd(vadd(vscale(D[19], b[6]), vscale(D[14], b[5])), b[4]);
            g[5] = vadd(vscale(D[20], b[6]), b[5]);
            g[6] = b[6];

            var pc_err = 1e30;
            var pc_err_last = 2.0;
            var iterations = 0u;
            loop {
                if (pc_err < PC_TOL || (iterations > 2u && pc_err_last <= pc_err) || iterations >= 12u) {
                    break;
                }
                pc_err_last = pc_err;
                pc_err = 0.0;
                iterations = iterations + 1u;
                for (var sub = 1u; sub < 8u; sub = sub + 1u) {
                    let hh = H[sub];
                    let xf7 = rmulf(rmul(dt, hh), 0.5);
                    let xf8 = rmul(dt, hh);
                    // Position and velocity at the substep from the series.
                    let f = 7u * (sub - 1u);
                    var dx = vscale(XF[f], b[6]);
                    var dv = vscale(VF[f], b[6]);
                    for (var k = 1u; k < 7u; k = k + 1u) {
                        dx = vscale(XF[f + k], vadd(dx, b[6u - k]));
                        dv = vscale(VF[f + k], vadd(dv, b[6u - k]));
                    }
                    dx = vscale(xf8, vadd(vscale(xf7, vadd(dx, a0)), v0));
                    dv = vscale(xf8, vadd(dv, a0));
                    let a = accel(vadd(x0, dx), radd(t, xf8), &cur);
                    // Update g[sub - 1] from the new acceleration, and the b's with it.
                    let m = sub - 1u;
                    let base = sub * (sub - 1u) / 2u;
                    var gg = vscale(RRI[base], vsub(a, a0));
                    for (var j = 0u; j < m; j = j + 1u) {
                        gg = vscale(RRI[base + j + 1u], vsub(gg, g[j]));
                    }
                    let tmp = vsub(gg, g[m]);
                    g[m] = gg;
                    if (m > 0u) {
                        let cb = (m - 1u) * m / 2u;
                        for (var j = 0u; j < m; j = j + 1u) {
                            b[j] = vaxpy(b[j], C[cb + j], tmp);
                        }
                    }
                    b[m] = vadd(b[m], tmp);
                    if (sub == 7u) {
                        let err = norm32(vto(tmp)) / norm32(vto(a));
                        if (err == err && abs(err) < 3.0e38) {
                            pc_err = err;
                        }
                    }
                }
            }

            // The next timestep from this particle's coefficients.
            let dtf = rto(dt);
            var dt_new = dtf / SAFETY;
            let a0f = vto(a0);
            let bf = array<vec3<f32>, 7>(vto(b[0]), vto(b[1]), vto(b[2]), vto(b[3]), vto(b[4]), vto(b[5]), vto(b[6]));
            let n2 = norm32(a0f + bf[0] + bf[1] + bf[2] + bf[3] + bf[4] + bf[5] + bf[6]);
            let n3 = norm32(bf[0] + 2.0 * bf[1] + 3.0 * bf[2] + 4.0 * bf[3] + 5.0 * bf[4] + 6.0 * bf[5] + 7.0 * bf[6]);
            let n4 = norm32(2.0 * bf[1] + 6.0 * bf[2] + 12.0 * bf[3] + 20.0 * bf[4] + 30.0 * bf[5] + 42.0 * bf[6]);
            if (norm32(a0f) > 0.0 && n2 > 0.0) {
                let r3 = n3 / n2;
                let r4 = n4 / n2;
                let ts2 = 2.0 / (r3 * r3 + r4);
                if (ts2 > 0.0 && ts2 < 3.0e38) {
                    dt_new = sqrt(ts2) * dtf * EPS_FAC;
                }
            }
            let ratio = abs(dt_new / dtf);
            if (ratio < SAFETY) {
                // Rejected: again with the shorter step.
                rejections = rejections + 1u;
                if (rejections > 30u || !(dt_new == dt_new) || dt_new == 0.0) {
                    break;
                }
                dt = rf(dt_new);
                clipped = false;
                continue;
            }
            if (ratio > 1.0 / SAFETY) {
                dt_new = dtf / SAFETY;
            }
            // Accepted.
            let dt2 = rmul(dt, dt);
            var px = vadd(vscale(INV72, b[6]), vscale(INV56, b[5]));
            px = vadd(px, vadd(vscale(INV42, b[4]), vscale(INV30, b[3])));
            px = vadd(px, vadd(vscale(INV20, b[2]), vscale(INV12, b[1])));
            px = vadd(px, vadd(vscale(INV6, b[0]), vscale(rf(0.5), a0)));
            var pv = vadd(vscale(INV8, b[6]), vscale(INV7, b[5]));
            pv = vadd(pv, vadd(vscale(INV6, b[4]), vscale(INV5, b[3])));
            pv = vadd(pv, vadd(vscale(INV4, b[2]), vscale(INV3, b[1])));
            pv = vadd(pv, vadd(vscale(rf(0.5), b[0]), a0));
            x0 = vadd(vadd(x0, vscale(dt, v0)), vscale(dt2, px));
            v0 = vadd(v0, vscale(dt, pv));
            if (clipped) {
                t = t_end;
                dt_new = rto(natural);
            } else {
                t = radd(t, dt);
            }
            dt_last = dt;
            dt = rf(dt_new);
            for (var k = 0u; k < 7u; k = k + 1u) {
                bl[k] = b[k];
                el[k] = e[k];
            }
            accepted = true;
            taken = taken + 1u;
            break;
        }
        if (!accepted || !(rto(x0.x) == rto(x0.x))) {
            state = 2u;
            break;
        }
    }
    if (state == 0u && !rlt(t, t_end)) {
        state = 1u;
    }

    put(o, x0);
    put(o + 3u, v0);
    st[o + 6u] = t;
    st[o + 7u] = dt;
    st[o + 8u] = dt_last;
    for (var k = 0u; k < 7u; k = k + 1u) {
        put(o + 9u + 3u * k, bl[k]);
        put(o + 30u + 3u * k, el[k]);
    }
    status[2u * i] = state;
    status[2u * i + 1u] = cur;
    if (state == 0u) {
        atomicAdd(&counts[0], 1u);
    }
    atomicAdd(&counts[1], taken);
}
