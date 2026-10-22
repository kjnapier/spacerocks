// Arithmetic self-test: for each pair (a[i], b[i]), out = a + b, a * b, a / b, sqrt(|a|).
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
@group(0) @binding(1) var<storage, read> xa: array<R>;
@group(0) @binding(2) var<storage, read> xb: array<R>;
@group(0) @binding(3) var<storage, read_write> out: array<R>;

fn MASK() -> u32 { return P.mask; }

@compute @workgroup_size(64)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    if (i >= P.n) {
        return;
    }
    let a = xa[i];
    let b = xb[i];
    out[4u * i] = radd(a, b);
    out[4u * i + 1u] = rmul(a, b);
    out[4u * i + 2u] = rdiv(a, b);
    out[4u * i + 3u] = rsqrt(rabs(a));
}
