"""Automatic non-gravitational model selection (fit_nongrav="auto"): layup vs spacerocks.

    LAYUP_TESTS=/path/to/layup/tests SPACEROCKS_KERNELS=/path/to/kernels python compare_nongrav_auto.py

Layup's synthetic Apophis-like arcs (tests/layup/test_nongrav_a2.py: 40 geocentric detections
over four years) with A1, A2, A3 scaled from its test values, plus 0.1" Gaussian noise so that
the reduced-chi-square gate is exercised. Both codes start from layup's perturbed seed: a
gravity-only fit, then the model ladder (layup's `_select_nongrav_auto`, spacerocks'
`nongrav="auto"`). Each case runs with layup's default thresholds and with the gate off
(accept_reduced_chi2 = 0).
"""
import os, sys
import numpy as np

sys.path.insert(0, os.path.join(os.environ["LAYUP_TESTS"], "layup"))
import test_nongrav_a2 as T
from layup import orbitfit as O
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
NAMES = ("A1", "A2", "A3")


def label(mask_bits):
    return "".join(n for n, b in zip(NAMES, (1, 2, 4)) if mask_bits & b) or "gravity"


rng = np.random.default_rng(357)
cases = []
for base in ((0, 0, 0), (0, 1, 0), (1, 0, 0), (0, 0, 1), (1, 1, 0), (1, 1, 1)):
    for scale in ((0.0,) if base == (0, 0, 0) else (30.0, 100.0, 300.0, 1000.0)):
        cases.append((base, scale))

agree = total = 0
print(f"{'truth':>28} {'gate':>5} | {'layup':>8} {'chi2':>10} | {'spacerocks':>10} {'chi2':>10} | max |dA|/sigma")
for base, scale in cases:
    a123 = tuple(scale * b * t for b, t in zip(base, (T._TRUE_A1, T._TRUE_A2, T._TRUE_A3)))
    clean = T._build_arc(a123=a123)
    noisy = []
    for o in clean:
        ra = np.arctan2(o.rho_hat[1], o.rho_hat[0])
        dec = np.arcsin(o.rho_hat[2])
        dra, ddec = rng.normal(0, SIG, 2)
        n = Observation.from_astrometry_with_id("x", ra + dra / np.cos(dec), dec + ddec, o.epoch,
                                                list(o.observer_position), list(o.observer_velocity))
        n.ra_unc = n.dec_unc = SIG
        noisy.append(n)
    ra = np.array([np.arctan2(o.rho_hat[1], o.rho_hat[0]) for o in noisy])
    dec = np.array([np.arcsin(o.rho_hat[2]) for o in noisy])
    ep = np.array([o.epoch for o in noisy])
    pos = np.array([o.observer_position for o in noisy])
    seed = T._perturbed_seed()
    rock = SpaceRock.from_xyz("seed", *seed.state, Time(seed.epoch, "tdb", "jd"), "J2000", "SSB")
    g = run_from_vector_with_initial_guess(ephem, seed, noisy, 100, 0)
    for gate in (1.5, 0.0):
        th = O.NongravAutoThresholds(accept_reduced_chi2=gate)
        L = O._select_nongrav_auto(ephem, g, noisy, th) if g.flag == 0 else g
        S = orbfit.fit(ra, dec, ep, pos, kernel, sigma_ra=SIG, sigma_dec=SIG, timescale="tdb", initial=rock,
                       nongrav="auto", nongrav_thresholds=(gate, 9.0, 3.0))
        lmask = L.nongrav_mask
        smask = sum(b for b, s in zip((1, 2, 4), S.nongrav_sigma) if np.isfinite(s))
        la = np.array([L.a1, L.a2, L.a3]) if lmask else np.zeros(3)
        ls = np.array([L.a1_unc, L.a2_unc, L.a3_unc])
        fitted = [k for k in range(3) if lmask & (1, 2, 4)[k]]
        da = max((abs(S.nongrav[k] - la[k]) / ls[k] for k in fitted), default=0.0)
        same = lmask == smask
        agree += same
        total += 1
        truth = " ".join(f"{v:.1e}" for v in a123)
        print(f"{truth:>28} {gate:5.1f} | {label(lmask):>8} {L.csq:10.3f} | {label(smask):>10} {S.chi2:10.3f} | {da:.2e}"
              + ("" if same else "   <-- differs"))
print(f"\nmodel choice agrees in {agree} of {total} cases")
