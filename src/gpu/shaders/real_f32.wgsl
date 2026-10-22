// Real numbers as plain f32 (about 7 significant digits).
alias R = f32;

fn rf(x: f32) -> R { return x; }
fn rto(a: R) -> f32 { return a; }
fn radd(a: R, b: R) -> R { return a + b; }
fn rsub(a: R, b: R) -> R { return a - b; }
fn rmul(a: R, b: R) -> R { return a * b; }
fn rmulf(a: R, b: f32) -> R { return a * b; }
fn rdiv(a: R, b: R) -> R { return a / b; }
fn rsqrt(a: R) -> R { return sqrt(a); }
fn rneg(a: R) -> R { return -a; }
fn rabs(a: R) -> R { return abs(a); }
fn rlt(a: R, b: R) -> bool { return a < b; }
fn rtrunc(a: R) -> R { return trunc(a); }
fn rround(a: R) -> R { return round(a); }
