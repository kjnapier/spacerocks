"""Herget IOD: layup's herget_iod vs spacerocks' orbfit.herget_iod on the primary arc (seq[0]) of
every object in layup runs (npz files from run_layup.py), as layup's do_fit calls it.

    LAYUP_CACHE=/path SPACEROCKS_KERNELS=/path python compare_herget.py run1.npz [run2.npz ...]

layup integrates with ASSIST's default force model (GR_EIH, Earth and Sun harmonics) and a plain
REBOUND IAS15; spacerocks uses its fitter's force model, so the seeds agree to the tolerance of the
iteration (0.003 AU in range), not to round-off. Reported: whether each side finds a seed, and the
relative position and velocity differences where both do.
"""
import os, sys, time
import numpy as np
import assist
import layup.orbitfit as LO
from layup.routines import Observation
from layup.utilities.herget_iod import herget_with_assist
from layup.utilities.universal_kepler import KeplerConvergenceError
from spacerocks import orbfit
from spacerocks.spice import SpiceKernel

K = os.environ["SPACEROCKS_KERNELS"]
C = os.environ["LAYUP_CACHE"]
kernel = SpiceKernel()
for f in ["latest_leapseconds.tls", "de440.bsp", "sb441-n16.bsp"]:
    kernel.load(os.path.join(K, f))
ephem = assist.Ephem(os.path.join(C, "linux_p1550p2650.440"), os.path.join(C, "sb441-n16.bsp"))


raised = []


def layup_herget(oid, O):
    """layup's herget_iod, with the ephemeris passed in (it otherwise builds one from its default
    cache). layup lets a KeplerConvergenceError from the two-body solver escape, which ends its
    whole orbitfit run; here, as in spacerocks, it only fails that starting range."""
    for rho in (2, 5, 40):
        try:
            # fresh copies: an escaping error leaves layup's observations light-time shifted
            s = herget_with_assist(O(), [list(range(len(O())))], ephem, initial_rho=rho)
        except KeplerConvergenceError:
            raised.append((oid, rho))
            continue
        if s:
            return s[0]
    return None


rows, t_l, t_s = [], 0.0, 0.0
for path in sys.argv[1:]:
    d = np.load(path, allow_pickle=True)
    obs, et, pv = d["obs"], d["et"], d["pv"]
    ids = obs["provID"].astype(str)
    for oid in dict.fromkeys(ids):
        m = np.where(ids == oid)[0]
        m = m[np.argsort(obs["obsTime"][m].astype(str), kind="mergesort")]
        jd = 2451545.0 + et[m] / 86400.0
        seq = LO._build_sequence(jd, 90.0)
        p = m[seq[0]]
        if len(p) < 3:
            continue
        pjd = 2451545.0 + et[p] / 86400.0
        O = lambda: [Observation.from_astrometry_with_id(oid, obs["ra"][i] * np.pi / 180.0, obs["dec"][i] * np.pi / 180.0, t,
                                                 [pv["x"][i], pv["y"][i], pv["z"][i]], [0, 0, 0]) for i, t in zip(p, pjd)]
        t0 = time.perf_counter()
        L = layup_herget(oid, O)
        t_l += time.perf_counter() - t0
        pos = np.array([[pv["x"][i], pv["y"][i], pv["z"][i]] for i in p])
        t0 = time.perf_counter()
        S = orbfit.herget_iod(obs["ra"][p] * np.pi / 180.0, obs["dec"][p] * np.pi / 180.0, pjd, pos, kernel, timescale="tdb")
        t_s += time.perf_counter() - t0
        if L is None or S is None:
            rows.append((oid, L is not None, S is not None, np.nan, np.nan, np.nan, np.nan))
            continue
        ls = np.array(L.state)
        ss = np.concatenate([S.position, S.velocity])
        rows.append((oid, True, True, np.linalg.norm(ss[:3] - ls[:3]) / np.linalg.norm(ls[:3]),
                     np.linalg.norm(ss[3:] - ls[3:]) / np.linalg.norm(ls[3:]), S.epoch.jd() - L.epoch, L.epoch))

both = [r for r in rows if r[1] and r[2]]
print(f"{len(rows)} primary arcs: layup finds {sum(r[1] for r in rows)}, spacerocks {sum(r[2] for r in rows)}, both {len(both)}")
for r in rows:
    if r[1] != r[2]:
        print(f"  {r[0]:>12}: only {'layup' if r[1] else 'spacerocks'}")
dp = np.array([r[3] for r in both]); dv = np.array([r[4] for r in both])
print(f"relative differences where both succeed: position median {np.median(dp):.1e}, max {dp.max():.1e}; "
      f"velocity median {np.median(dv):.1e}, max {dv.max():.1e}; epochs differ by at most {max(abs(r[5]) for r in both):.1e} d")
for r in sorted(both, key=lambda r: -abs(r[5]))[:3]:
    print(f"  epoch {r[0]:>12}  {r[5]:.3e} at {r[6]}")
for r in sorted(both, key=lambda r: -r[3])[:5]:
    print(f"  {r[0]:>12}  {r[3]:.2e}  {r[4]:.2e}")
print(f"layup raised KeplerConvergenceError {len(raised)} times (arc, starting range): {raised[:8]}")
print(f"time: layup {t_l:.1f} s, spacerocks {t_s:.1f} s")
