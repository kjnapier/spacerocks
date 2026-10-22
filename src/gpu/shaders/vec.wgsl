// Three-vectors of R, for every precision.
struct V3 {
    x: R,
    y: R,
    z: R,
}

fn v3(x: R, y: R, z: R) -> V3 { return V3(x, y, z); }
fn vadd(a: V3, b: V3) -> V3 { return V3(radd(a.x, b.x), radd(a.y, b.y), radd(a.z, b.z)); }
fn vsub(a: V3, b: V3) -> V3 { return V3(rsub(a.x, b.x), rsub(a.y, b.y), rsub(a.z, b.z)); }
fn vscale(s: R, a: V3) -> V3 { return V3(rmul(s, a.x), rmul(s, a.y), rmul(s, a.z)); }
fn vdot(a: V3, b: V3) -> R { return radd(radd(rmul(a.x, b.x), rmul(a.y, b.y)), rmul(a.z, b.z)); }
fn vzero() -> V3 { return V3(rf(0.0), rf(0.0), rf(0.0)); }
// a + s b
fn vaxpy(a: V3, s: R, b: V3) -> V3 { return vadd(a, vscale(s, b)); }
// f32 copy of a vector, for estimates that don't need the full precision.
fn vto(a: V3) -> vec3<f32> { return vec3<f32>(rto(a.x), rto(a.y), rto(a.z)); }
