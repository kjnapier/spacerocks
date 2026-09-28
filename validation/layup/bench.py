"""Benchmark: layup's orbit-fit pipeline vs spacerocks', per object, from Python.

    LAYUP_CACHE=/path SPACEROCKS_KERNELS=/path RAYON_NUM_THREADS=1 python bench.py out.json run1.npz [run2.npz ...]

Both codes get the same inputs: layup's TDB epochs and barycentric observer positions (from
run_layup.py's npz), with RA/Dec converted as layup does (x * pi / 180). Setup is excluded on
both sides: layup's `get_ephem` is memoized (its do_fit otherwise reloads the ASSIST ephemeris on
every call), and spacerocks' kernel is loaded once. Timed per object, best of `REPEATS`:

  pipeline  layup `do_fit(iod="auto")`            vs  `orbfit.fit(...)`          (IOD -> fit)
  lm        `run_from_vector_with_initial_guess`  vs  `orbfit.fit(initial=...)`  from layup's
            converged state perturbed by 1e-6 AU / 1e-8 AU/d, full detection set

Set RAYON_NUM_THREADS=1 for a per-core comparison; layup is single-threaded.
"""
import json, os, sys, time
import numpy as np

import layup.orbitfit as LO
from layup.routines import FitResult, Observation, get_ephem, run_from_vector_with_initial_guess

from spacerocks import SpaceRock, orbfit
from spacerocks.spice import SpiceKernel
from spacerocks.time import Time

REPEATS = 3
CACHE = os.environ["LAYUP_CACHE"]

# --- setup, timed separately ------------------------------------------------------------------
t = time.perf_counter()
EPHEM = get_ephem(CACHE)
t_layup_setup = time.perf_counter() - t
LO.get_ephem = lambda cache_dir: EPHEM  # do_fit calls it twice per object

t = time.perf_counter()
kernel = SpiceKernel()
for f in ["latest_leapseconds.tls", "de440.bsp", "sb441-n16.bsp"]:
    kernel.load(os.path.join(os.environ["SPACEROCKS_KERNELS"], f))
t_sr_setup = time.perf_counter() - t


def best(fn):
    out, ts = None, []
    for _ in range(REPEATS):
        t = time.perf_counter()
        out = fn()
        ts.append(time.perf_counter() - t)
    return out, min(ts)


rows = []
for path in sys.argv[2:]:
    d = np.load(path, allow_pickle=True)
    obs, et, pv, fits = d["obs"], d["et"], d["pv"], d["fits"]
    ids = obs["provID"].astype(str)
    ref = {str(r["provID"]): r for r in fits}
    for oid in dict.fromkeys(ids):
        m = np.where(ids == oid)[0]
        m = m[np.argsort(obs["obsTime"][m].astype(str), kind="mergesort")]
        if len(m) < 3:
            continue
        jd = 2451545.0 + et[m] / 86400.0
        ra, dec = obs["ra"][m] * np.pi / 180.0, obs["dec"][m] * np.pi / 180.0
        pos = np.column_stack([pv["x"][m], pv["y"][m], pv["z"][m]])
        O = [Observation.from_astrometry_with_id(oid, ra[i], dec[i], jd[i], list(pos[i]),
                                                 [pv["vx"][m][i], pv["vy"][m][i], pv["vz"][m][i]]) for i in range(len(m))]

        def layup_pipeline():
            return LO.do_fit(O, LO._build_sequence(jd, 90.0), CACHE, iod="auto")

        L, tl = best(layup_pipeline)
        S, ts = best(lambda: orbfit.fit(ra, dec, jd, pos, kernel, timescale="tdb"))
        row = dict(set=os.path.basename(path), id=oid, n=len(m), arc=float(jd[-1] - jd[0]),
                   layup_s=tl, sr_s=ts, layup_flag=int(L.flag), sr_flag=int(S.flag))

        r = ref[oid]
        if r["flag"] == 0:
            s0 = np.array([r[c] for c in ["x", "y", "z", "xdot", "ydot", "zdot"]])
            s0[:3] += 1e-6
            s0[3:] += 1e-8
            e0 = float(r["epochMJD_TDB"]) + 2400000.5
            g = FitResult()
            g.state, g.epoch, g.flag = list(s0), e0, 0
            rock = SpaceRock.from_xyz("x", *s0, Time(e0, "tdb", "jd"), "J2000", "SSB")
            Ll, tll = best(lambda: run_from_vector_with_initial_guess(EPHEM, g, O, 100))
            Sl, tsl = best(lambda: orbfit.fit(ra, dec, jd, pos, kernel, timescale="tdb", initial=rock))
            row.update(layup_lm_s=tll, sr_lm_s=tsl, layup_lm_niter=int(Ll.niter), sr_lm_niter=int(Sl.niter),
                       lm_flags=[int(Ll.flag), int(Sl.flag)])
        rows.append(row)
        print(f"{row['set'][:18]:18} {oid:>12} n={len(m):5d} pipeline {tl:8.4f} / {ts:8.4f} s"
              + (f"   lm {row['layup_lm_s']:8.4f} / {row['sr_lm_s']:8.4f} s" if "layup_lm_s" in row else ""), flush=True)

json.dump(dict(rows=rows, layup_setup_s=t_layup_setup, sr_setup_s=t_sr_setup, threads=os.environ.get("RAYON_NUM_THREADS")),
          open(sys.argv[1], "w"), indent=1)
print(f"setup: layup get_ephem {t_layup_setup:.3f} s, spacerocks kernel load {t_sr_setup:.3f} s")
