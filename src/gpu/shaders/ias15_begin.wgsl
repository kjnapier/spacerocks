// Start a block: every test particle at time 0 of it, running, in the table's first row.
@compute @workgroup_size(64)
fn begin(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    if (i >= P.n) {
        return;
    }
    st[ST * i + 6u] = rf(0.0);
    status[2u * i] = 0u;
    status[2u * i + 1u] = 0u;
}
