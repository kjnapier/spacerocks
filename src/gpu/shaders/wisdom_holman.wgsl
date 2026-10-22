// Wisdom–Holman steps of test particles through a block of recorded steps of the massive bodies
// (see `MassiveRun` in nbody/integrators/wisdom_holman.rs). One invocation per test particle,
// state in democratic heliocentric coordinates (heliocentric q, barycentric u). Each step:
//
// kick(h/2) jump(h/2) kepler(h) jump(h/2) kick(h/2)
//
// where the kicks are the pull of the massive bodies other than the central one at their
// recorded positions, and the jumps the recorded shifts. Needs FAST_KICK (bool): kicks in f32.

struct Params {
    n: u32,
    mask: u32,
    // Steps k0..k1 of the table.
    k0: u32,
    k1: u32,
    // Bodies pulling in the kicks.
    np: u32,
    c: u32,
    d: u32,
    e: u32,
}

@group(0) @binding(0) var<uniform> P: Params;
// q, u per particle.
@group(0) @binding(1) var<storage, read_write> st: array<R>;
// Per step: dq1, dq2, then the pulling bodies' q at the first kick, then at the second.
@group(0) @binding(2) var<storage, read> tab: array<R>;
// G m of the central body, h, h / 2, then G m of each pulling body.
@group(0) @binding(3) var<storage, read> coef: array<R>;
// Nonzero once a particle's Kepler solve failed; it isn't stepped after that.
@group(0) @binding(4) var<storage, read_write> flags: array<u32>;

fn MASK() -> u32 { return P.mask; }

fn tab3(i: u32) -> V3 {
    return V3(tab[i], tab[i + 1u], tab[i + 2u]);
}

// u after a kick of length `half` from the bodies at tab[at..].
fn kick(q: V3, u: V3, at: u32, half: R) -> V3 {
    if (FAST_KICK) {
        var acc = vec3<f32>(0.0);
        for (var p = 0u; p < P.np; p = p + 1u) {
            let d = vto(vsub(q, tab3(at + 3u * p)));
            let r2 = dot(d, d);
            acc = acc - rto(coef[3u + p]) / (r2 * sqrt(r2)) * d;
        }
        return V3(radd(u.x, rmulf(half, acc.x)), radd(u.y, rmulf(half, acc.y)), radd(u.z, rmulf(half, acc.z)));
    }
    var acc = vzero();
    for (var p = 0u; p < P.np; p = p + 1u) {
        let d = vsub(q, tab3(at + 3u * p));
        let r2 = vdot(d, d);
        let s = rdiv(coef[3u + p], rmul(r2, rsqrt(r2)));
        acc = vsub(acc, vscale(s, d));
    }
    return vaxpy(u, half, acc);
}

@compute @workgroup_size(64)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    if (i >= P.n || flags[i] != 0u) {
        return;
    }
    let gm = coef[0];
    let h = coef[1];
    let half = coef[2];
    let stride = 6u + 6u * P.np;
    var q = V3(st[6u * i], st[6u * i + 1u], st[6u * i + 2u]);
    var u = V3(st[6u * i + 3u], st[6u * i + 4u], st[6u * i + 5u]);
    for (var k = P.k0; k < P.k1; k = k + 1u) {
        let base = k * stride;
        u = kick(q, u, base + 6u, half);
        q = vadd(q, tab3(base));
        let kp = kepler_step(q, u, gm, h);
        if (!kp.ok) {
            flags[i] = k + 1u;
            break;
        }
        q = vadd(kp.q, tab3(base + 3u));
        u = kick(q, kp.u, base + 6u + 3u * P.np, half);
    }
    st[6u * i] = q.x;
    st[6u * i + 1u] = q.y;
    st[6u * i + 2u] = q.z;
    st[6u * i + 3u] = u.x;
    st[6u * i + 4u] = u.y;
    st[6u * i + 5u] = u.z;
}
