"""Per-arc (piecewise-constant) non-gravitational parameters: layup vs spacerocks.

    LAYUP_TESTS=/path/to/layup/tests SPACEROCKS_KERNELS=/path/to/kernels python compare_per_arc.py

Arcs as in layup's tests/layup/test_nongrav_per_arc.py: 48 geocentric detections over +-3 years
around the fit epoch, generated from one state with one set of A1, A2, A3 before the epoch and
another after. Both codes run layup's `_orbitfit` sequence from the same seed: a gravity-only fit,
then the joint per-arc fit starting from it.
"""
import os, sys
import numpy as np

sys.path.insert(0, os.path.join(os.environ["LAYUP_TESTS"], "layup"))
import test_nongrav_per_arc as T
from layup.routines import Observation, get_ephem, run_from_vector_with_initial_guess

from spacerocks import SpaceRock, orbfit
from spacerocks.spice import SpiceKernel
from spacerocks.time import Time

K = os.environ["SPACEROCKS_KERNELS"]
kernel = SpiceKernel()
for f in ["latest_leapseconds.tls", "de440.bsp", "sb441-n16.bsp", "earth_1962_240827_2124_combined.bpc"]:
    kernel.load(os.path.join(K, f))
ephem = get_ephem(T.CACHE)
SIG = 0.1 / 206265
BITS = {"A1": 1, "A2": 2, "A3": 4}


def build(a_before, a_after, noise, rng):
    """layup's _build_piecewise_arc, for any (A1, A2, A3) on each side, plus optional noise."""
    import assist, rebound
    ep = assist.Ephem(os.path.join(T.CACHE, T._EPHEM[0]), os.path.join(T.CACHE, T._EPHEM[1]))
    jr = ep.jd_ref

    def ast_at(t, a123):
        sim = rebound.Simulation()
        sim.t = T._EPOCH - jr
        sim.add(x=T._STATE[0], y=T._STATE[1], z=T._STATE[2], vx=T._STATE[3], vy=T._STATE[4], vz=T._STATE[5])
        ax = assist.Extras(sim, ep)
        ax.forces = ["SUN", "PLANETS", "ASTEROIDS", "NON_GRAVITATIONAL", "GR_SIMPLE"]
        for k, v in T._GR.items():
            setattr(ax, k, v)
        ax.particle_params = np.array(a123, dtype=float)
        sim.integrate(t - jr)
        p = sim.particles[0]
        return np.array([p.x, p.y, p.z])

    obs = []
    for dt in np.linspace(-3 * 365.25, 3 * 365.25, 48):
        if abs(dt) < 1.0:
            continue
        a123 = a_before if dt < 0 else a_after
        t = T._EPOCH + dt
        e = ep.get_particle("Earth", t - jr)
        r_obs, v_obs = np.array([e.x, e.y, e.z]), np.array([e.vx, e.vy, e.vz])
        lt = 0.0
        for _ in range(3):
            rho = ast_at(t - lt, a123) - r_obs
            lt = np.linalg.norm(rho) / T._C
        rho /= np.linalg.norm(rho)
        ra, dec = np.arctan2(rho[1], rho[0]), np.arcsin(rho[2])
        if noise:
            dra, ddec = rng.normal(0, SIG, 2)
            ra, dec = ra + dra / np.cos(dec), dec + ddec
        o = Observation.from_astrometry_with_id("c", ra, dec, t, list(r_obs), list(v_obs))
        o.ra_unc = o.dec_unc = SIG
        obs.append(o)
    return obs


rng = np.random.default_rng(351)
cases = [("A2", (0, -5e-14, 0), (0, -5e-14, 0)), ("A2", (0, 1e-12, 0), (0, 5e-12, 0)),
         ("A1", (2e-11, 0, 0), (6e-11, 0, 0)), ("A3", (0, 0, 2e-11), (0, 0, -2e-11)),
         ("A1A2A3", (3e-11, 2e-12, 1e-11), (1e-11, 6e-12, -1e-11))]
print(f"{'fit':7} {'noise':5} {'code':10} {'flag':>4} {'chi2':>10}  arc A {'':26} arc B")
for mask_name, before, after in cases:
    for noise in (False, True):
        obs = build(before, after, noise, rng)
        mask = sum(BITS[p] for p in BITS if p in mask_name)
        seed = T._seed()
        g = run_from_vector_with_initial_guess(ephem, seed, obs, 100, 0)
        L = run_from_vector_with_initial_guess(ephem, g, obs, 100, mask, [], True)

        ra = np.array([np.arctan2(o.rho_hat[1], o.rho_hat[0]) for o in obs])
        dec = np.array([np.arcsin(o.rho_hat[2]) for o in obs])
        jd = np.array([o.epoch for o in obs])
        pos = np.array([list(o.observer_position) + list(o.observer_velocity) for o in obs])
        rock = SpaceRock.from_xyz("seed", *seed.state, Time(seed.epoch, "tdb", "jd"), "J2000", "SSB")
        S = orbfit.fit(ra, dec, jd, pos, kernel, sigma_ra=SIG, sigma_dec=SIG, timescale="tdb", initial=rock,
                       nongrav=mask_name, per_arc=True)

        la, lb = np.array([L.a1, L.a2, L.a3]), np.array([L.a1_arc2, L.a2_arc2, L.a3_arc2])
        lsa, lsb = np.array([L.a1_unc, L.a2_unc, L.a3_unc]), np.array([L.a1_arc2_unc, L.a2_arc2_unc, L.a3_arc2_unc])
        on = [k for k in range(3) if mask & (1 << k)]
        da = max(abs(S.nongrav[k] - la[k]) / lsa[k] for k in on)
        db = max(abs(S.nongrav_arc2[k] - lb[k]) / lsb[k] for k in on)
        cov = np.array(L.cov).reshape(6, 6) if len(L.cov) == 36 else None
        ds = S.state - np.array(L.state)
        dstate = np.sqrt(ds @ np.linalg.solve(np.array(S.state_covariance), ds))
        fmt = lambda v: " ".join(f"{v[k]:10.3e}" for k in on)
        print(f"{mask_name:7} {str(noise):5} {'layup':10} {L.flag:4d} {L.csq:10.4g}  {fmt(la):32} {fmt(lb)}")
        print(f"{'':7} {'':5} {'spacerocks':10} {S.flag:4d} {S.chi2:10.4g}  {fmt(S.nongrav):32} {fmt(S.nongrav_arc2)}")
        print(f"{'':13} differences: arc A {da:.1e} sigma, arc B {db:.1e} sigma, state {dstate:.1e} sigma, "
              f"sigma ratio {np.nanmax(np.abs(S.nongrav_sigma[on] / lsa[on] - 1)):.1e}")
