// The universal-variable Kepler step of `transforms::universal_kepler_step`, in R. The starting
// guess is computed in f32 (only Halley's iteration needs the full precision); on unbound orbits
// it starts from the hyperbolic Kepler equation instead of a series, since this solver has no
// bracketing fallback. A step that doesn't converge comes back with ok = false. Needs the
// constants KEPLER_K_TOL, KEPLER_DX_TOL and KEPLER_STALL (f32, scaled to the precision's
// round-off) and TAU (R).

struct Stumpff {
    c0: R,
    c1: R,
    c2: R,
    c3: R,
}

// Stumpff functions c0..c3 at z, by series at a reduced argument and the quadruple-argument
// identities.
fn stumpff(z0: R) -> Stumpff {
    var z = z0;
    var n = 0u;
    while (abs(rto(z)) > 0.1) {
        z = rmulf(z, 0.25);
        n = n + 1u;
    }
    let c2 = rsub(INV_F2, rmul(z, rsub(INV_F4, rmul(z, rsub(INV_F6, rmul(z, rsub(INV_F8, rmul(z, rsub(INV_F10, rmul(z, rsub(INV_F12, rmul(z, INV_F14))))))))))));
    let c3 = rsub(INV_F3, rmul(z, rsub(INV_F5, rmul(z, rsub(INV_F7, rmul(z, rsub(INV_F9, rmul(z, rsub(INV_F11, rmul(z, rsub(INV_F13, rmul(z, INV_F15))))))))))));
    var s = Stumpff(rsub(rf(1.0), rmul(z, c2)), rsub(rf(1.0), rmul(z, c3)), c2, c3);
    for (var k = 0u; k < n; k = k + 1u) {
        s.c3 = rmulf(radd(s.c2, rmul(s.c0, s.c3)), 0.25);
        s.c2 = rmulf(rmul(s.c1, s.c1), 0.5);
        s.c1 = rmul(s.c1, s.c0);
        s.c0 = rsub(rmulf(rmul(s.c0, s.c0), 2.0), rf(1.0));
    }
    return s;
}

struct Kepler {
    q: V3,
    u: V3,
    ok: bool,
}

fn kepler_step(r0: V3, v0: V3, gm: R, dt: R) -> Kepler {
    let r0n = rsqrt(vdot(r0, r0));
    let v2 = vdot(v0, v0);
    let eta0 = vdot(r0, v0);
    let beta = rsub(rdiv(rmulf(gm, 2.0), r0n), v2);
    let zeta0 = rsub(gm, rmul(beta, r0n));

    // Whole periods of a bound orbit change nothing.
    var dt_red = dt;
    let tau_mu = rmul(TAU, gm);
    let bf = rto(beta);
    let tf = rto(dt);
    let tmf = rto(tau_mu);
    if (bf > 0.0 && tf * tf * bf * bf * bf > tmf * tmf) {
        let period = rdiv(tau_mu, rmul(beta, rsqrt(beta)));
        dt_red = rsub(dt, rmul(period, rtrunc(rdiv(dt, period))));
    }

    // Starting guess, in f32.
    let t = rto(dt_red);
    let r0f = rto(r0n);
    let ir = 1.0 / r0f;
    let etaf = rto(eta0);
    let zetaf = rto(zeta0);
    let gmf = rto(gm);
    var xg: f32;
    if (t * t * rto(v2) < 0.09 * r0f * r0f) {
        xg = t * ir * (1.0 - t * ir * ir * (0.5 * etaf - t * ir * (0.5 * etaf * etaf * ir - zetaf / 6.0)));
    } else if (bf > 0.0) {
        let sb = sqrt(bf);
        let ecos = 1.0 - r0f * bf / gmf;
        let esin = etaf * sb / gmf;
        let e0 = atan2(esin, ecos);
        let n_dt = bf * sb / gmf * t;
        var m = e0 - esin + n_dt;
        let tau = 6.2831853071795864;
        m = m - tau * round(m / tau);
        let e1 = m + 0.85 * sqrt(ecos * ecos + esin * esin) * sign(sin(m));
        var de = e1 - e0;
        de = de + tau * round((n_dt - de) / tau);
        xg = de / sb;
    } else {
        // Unbound: the usual starter for the hyperbolic Kepler equation e sinh F - F = N.
        let sb = sqrt(-bf);
        let ecosh = 1.0 - r0f * bf / gmf;
        let esinh = etaf * sb / gmf;
        let e = sqrt(max(ecosh * ecosh - esinh * esinh, 1e-12));
        let f0 = asinh(esinh / e);
        let n_end = esinh - f0 + sb * sb * sb / gmf * t;
        let f1 = sign(n_end) * log(2.0 * abs(n_end) / e + 1.8);
        xg = (f1 - f0) / sb;
    }

    var x = rf(xg);
    var converged = false;
    var last_dx = 3.0e38;
    var g1 = rf(0.0);
    var g2 = rf(0.0);
    var g3 = rf(0.0);
    for (var it = 0; it < 30; it = it + 1) {
        let c = stumpff(rmul(beta, rmul(x, x)));
        let g0 = c.c0;
        let x2 = rmul(x, x);
        g1 = rmul(x, c.c1);
        g2 = rmul(x2, c.c2);
        g3 = rmul(rmul(x2, x), c.c3);
        let f = rsub(radd(radd(rmul(r0n, x), rmul(eta0, g2)), rmul(zeta0, g3)), dt_red);
        let fp = radd(radd(r0n, rmul(eta0, g1)), rmul(zeta0, g2));
        let fpp = radd(rmul(eta0, g0), rmul(zeta0, g1));
        let ifp = rdiv(rf(1.0), fp);
        let dx = rdiv(f, rsub(fp, rmul(rmulf(f, 0.5), rmul(fpp, ifp))));
        x = rsub(x, dx);
        let adx = abs(rto(dx));
        let ax = abs(rto(x));
        // Halley's error after this correction is about K dx^3; once that is below round-off,
        // carry the g-functions to the new x by their Taylor series.
        if (adx <= 1e-4 * ax) {
            let fppp = rto(zeta0) * rto(g0) - bf * etaf * rto(g1);
            let ifpf = rto(ifp);
            let fppf = rto(fpp);
            let k = abs(ifpf * (0.25 * fppf * fppf * ifpf - fppp / 6.0));
            if (10.0 * k * adx * adx * adx <= KEPLER_K_TOL * ax) {
                let d = rneg(dx);
                let d2 = rmulf(rmul(d, d), 0.5);
                let d3 = rdiv(rmul(rmul(d, d), d), rf(6.0));
                let n1 = rsub(radd(g1, rmul(d, g0)), rmul(beta, radd(rmul(d2, g1), rmul(d3, g0))));
                let n2 = rsub(radd(radd(g2, rmul(d, g1)), rmul(d2, g0)), rmul(d3, rmul(beta, g1)));
                let n3 = radd(radd(radd(g3, rmul(d, g2)), rmul(d2, g1)), rmul(d3, g0));
                g1 = n1;
                g2 = n2;
                g3 = n3;
                converged = true;
                break;
            }
        }
        if (adx <= KEPLER_DX_TOL * ax || (adx >= last_dx && adx <= KEPLER_STALL * ax)) {
            let c = stumpff(rmul(beta, rmul(x, x)));
            let x2 = rmul(x, x);
            g1 = rmul(x, c.c1);
            g2 = rmul(x2, c.c2);
            g3 = rmul(rmul(x2, x), c.c3);
            converged = true;
            break;
        }
        last_dx = adx;
    }

    let irn = rdiv(rf(1.0), radd(radd(r0n, rmul(eta0, g1)), rmul(zeta0, g2)));
    let gm_r0 = rdiv(gm, r0n);
    let f = rsub(rf(1.0), rmul(gm_r0, g2));
    let g = rsub(dt_red, rmul(gm, g3));
    let fdot = rneg(rmul(rmul(gm_r0, irn), g1));
    let gdot = rsub(rf(1.0), rmul(rmul(gm, irn), g2));
    let q = vadd(vscale(f, r0), vscale(g, v0));
    let u = vadd(vscale(fdot, r0), vscale(gdot, v0));
    let ok = converged && abs(rto(q.x)) < 3.0e38 && abs(rto(u.x)) < 3.0e38 && rto(q.x) == rto(q.x);
    return Kepler(q, u, ok);
}
