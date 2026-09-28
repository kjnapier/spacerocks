"""Bernstein-Khushalani IOD: layup's run_bk_iod vs spacerocks' orbfit.bk_iod on the primary arc
(seq[0]) of every object in layup runs (npz files from run_layup.py), at the arc's middle
detection, as layup's do_fit calls it.

    SPACEROCKS_KERNELS=/path python compare_bk.py run1.npz [run2.npz ...]
"""
import os, sys
import numpy as np
import layup.orbitfit as LO
from layup.routines import Observation, run_bk_iod
from layup.constants import MU_SUN
from spacerocks import orbfit
from spacerocks.spice import SpiceKernel

kernel = SpiceKernel()
kernel.load(os.path.join(os.environ["SPACEROCKS_KERNELS"], "de440.bsp"))

worst = []
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
        O = [Observation.from_astrometry_with_id(oid, obs["ra"][i] * np.pi / 180.0, obs["dec"][i] * np.pi / 180.0, t,
                                                 [pv["x"][i], pv["y"][i], pv["z"][i]], [0, 0, 0]) for i, t in zip(p, pjd)]
        e = float(O[len(O) // 2].epoch)
        L = run_bk_iod(O, e, MU_SUN)
        pos = np.array([[pv["x"][i], pv["y"][i], pv["z"][i]] for i in p])
        S = orbfit.bk_iod(obs["ra"][p] * np.pi / 180.0, obs["dec"][p] * np.pi / 180.0, pjd, pos, kernel, timescale="tdb", epoch0=e)
        if L.flag != 0 or S is None:
            print(oid, "layup flag", L.flag, "spacerocks", S)
            continue
        ls = np.array(L.state)
        ss = np.concatenate([S.position, S.velocity])
        rel = np.abs(ss - ls) / np.abs(ls).max()
        worst.append((np.linalg.norm(ss[:3] - ls[:3]) / np.linalg.norm(ls[:3]), np.linalg.norm(ss[3:] - ls[3:]) / np.linalg.norm(ls[3:]), oid))
worst.sort(reverse=True)
print(f"{len(worst)} arcs; largest relative differences (position, velocity):")
for w in worst[:5]:
    print(f"  {w[2]:>10}  {w[0]:.2e}  {w[1]:.2e}")
