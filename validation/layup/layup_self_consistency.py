"""How reproducible is layup itself? Run layup's do_fit on every object twice, with RA/Dec
converted to radians as `x * pi / 180` (what layup's _orbitfit does) and as `np.radians(x)`
(which can differ in the last bit), and report objects whose outcome changes.

    LAYUP_CACHE=/path python layup_self_consistency.py run.npz out.npz
"""
import os, sys
import numpy as np
import layup.orbitfit as LO
from layup.routines import Observation

d = np.load(sys.argv[1], allow_pickle=True)
obs, et, pv = d["obs"], d["et"], d["pv"]
ids = obs["provID"].astype(str)
res = {}
for oid in dict.fromkeys(ids):
    m = np.where(ids == oid)[0]
    m = m[np.argsort(obs["obsTime"][m].astype(str), kind="mergesort")]
    jd = 2451545.0 + et[m] / 86400.0
    if len(m) < 3:
        continue
    seq = LO._build_sequence(jd, 90.0)
    out = []
    for conv in (lambda x: x * np.pi / 180.0, np.radians):
        O = [Observation.from_astrometry_with_id(oid, conv(obs["ra"][i]), conv(obs["dec"][i]), t,
                                                 [pv["x"][i], pv["y"][i], pv["z"][i]], [pv["vx"][i], pv["vy"][i], pv["vz"][i]]) for i, t in zip(m, jd)]
        x = LO.do_fit(O, seq, os.environ["LAYUP_CACHE"], iod="auto")
        out.append((x.flag, x.csq, x.epoch, np.array(x.state), np.array(x.cov).reshape(6, 6)))
    (f1, c1, e1, s1, cov1), (f2, c2, e2, s2, cov2) = out
    note = ""
    if f1 == f2 == 0:
        dx = s2 - s1
        note = f"{np.sqrt(dx @ np.linalg.solve(cov1, dx)):.2e} sigma" if e1 == e2 else f"epoch differs by {e2 - e1:.2e} d"
    print(f"{oid:>12} flags {f1}/{f2}  csq {c1:.6g}/{c2:.6g}  {note}", flush=True)
    res[oid] = (f1, f2, c1, c2, e1, e2)
np.savez(sys.argv[2], res=res)
