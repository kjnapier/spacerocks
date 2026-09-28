"""Non-gravitational fits: layup vs spacerocks on layup's synthetic Apophis-like arcs
(tests/layup/test_nongrav_a2.py), fitting A1, A2, A3 and all three.

    LAYUP_TESTS=/path/to/layup/tests SPACEROCKS_KERNELS=/path/to/kernels python compare_nongrav.py

Both codes run layup's `_orbitfit` sequence from a perturbed seed: a gravity-only fit, then the
joint state + non-grav fit starting from it.
"""
import os, sys
import numpy as np

sys.path.insert(0, os.path.join(os.environ["LAYUP_TESTS"], "layup"))
import test_nongrav_a2 as T
from layup.routines import get_ephem, run_from_vector_with_initial_guess

from spacerocks import SpaceRock, orbfit
from spacerocks.spice import SpiceKernel
from spacerocks.time import Time

K = os.environ["SPACEROCKS_KERNELS"]
kernel = SpiceKernel()
for f in ["latest_leapseconds.tls", "de440.bsp", "sb441-n16.bsp", "earth_1962_240827_2124_combined.bpc"]:
    kernel.load(os.path.join(K, f))
ephem = get_ephem(T.CACHE)

cases = [("A2", (0.0, T._TRUE_A2, 0.0)), ("A1", (T._TRUE_A1, 0.0, 0.0)), ("A3", (0.0, 0.0, T._TRUE_A3)),
         ("A1A2A3", (T._TRUE_A1, T._TRUE_A2, T._TRUE_A3))]
print(f"{'fit':8} {'code':10} {'flag':>4} {'niter':>5} {'chi2':>11} {'A1':>12} {'A2':>12} {'A3':>12}  sigma(A)")
for name, a123 in cases:
    obs = T._build_arc(a123=a123)
    mask = sum(T._BIT[p] for p in ("A1", "A2", "A3") if p in name)
    seed = T._perturbed_seed()
    g = run_from_vector_with_initial_guess(ephem, seed, obs, 100, 0)
    L = run_from_vector_with_initial_guess(ephem, g, obs, 100, mask)

    ra = np.array([np.arctan2(o.rho_hat[1], o.rho_hat[0]) for o in obs])
    dec = np.array([np.arcsin(o.rho_hat[2]) for o in obs])
    ep = np.array([o.epoch for o in obs])
    pos = np.array([o.observer_position for o in obs])
    s = np.array(seed.state)
    rock = SpaceRock.from_xyz("seed", *s, Time(seed.epoch, "tdb", "jd"), "J2000", "SSB")
    S = orbfit.fit(ra, dec, ep, pos, kernel, sigma_ra=0.1 / 206265, sigma_dec=0.1 / 206265, timescale="tdb", initial=rock, nongrav=name)

    lsig = [L.a1_unc, L.a2_unc, L.a3_unc]
    print(f"{name:8} {'layup':10} {L.flag:4d} {L.niter:5d} {L.csq:11.4e} {L.a1:12.5e} {L.a2:12.5e} {L.a3:12.5e}  " + " ".join(f"{v:.3e}" for v in lsig))
    print(f"{'':8} {'spacerocks':10} {S.flag:4d} {S.niter:5d} {S.chi2:11.4e} " + " ".join(f"{v:12.5e}" for v in S.nongrav) + "  " + " ".join(f"{v:.3e}" for v in S.nongrav_sigma))
    cov = np.array(L.cov).reshape(6, 6)
    dx = S.state - np.array(L.state)
    print(f"{'':8} state difference {np.sqrt(dx @ np.linalg.solve(cov, dx)):.2e} sigma (Mahalanobis), {np.abs(dx[:3]).max():.2e} AU")
