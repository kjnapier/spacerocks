// One Kepler step per case: in = r0, v0, gm, dt; out = r, v, 1 if it converged (else 0).
struct Params {
    n: u32,
    mask: u32,
    a: u32,
    b: u32,
    c: u32,
    d: u32,
    e: u32,
    f: u32,
}

@group(0) @binding(0) var<uniform> P: Params;
@group(0) @binding(1) var<storage, read> cases: array<R>;
@group(0) @binding(2) var<storage, read_write> out: array<R>;

fn MASK() -> u32 { return P.mask; }

@compute @workgroup_size(64)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    if (i >= P.n) {
        return;
    }
    let b = 8u * i;
    let k = kepler_step(V3(cases[b], cases[b + 1u], cases[b + 2u]), V3(cases[b + 3u], cases[b + 4u], cases[b + 5u]), cases[b + 6u], cases[b + 7u]);
    let o = 7u * i;
    out[o] = k.q.x;
    out[o + 1u] = k.q.y;
    out[o + 2u] = k.q.z;
    out[o + 3u] = k.u.x;
    out[o + 4u] = k.u.y;
    out[o + 5u] = k.u.z;
    out[o + 6u] = rf(select(0.0, 1.0, k.ok));
}
