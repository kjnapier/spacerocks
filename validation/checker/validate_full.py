"""Timing on the full MPCORB catalog (~1.4 million orbits), far from its epoch.

Loads MPCORB (parsed once, then from the binary cache), makes a snapshot 16 months after the
catalog's epoch (N-body, every orbit), saves and reloads it, and checks the September 2026
detections of (3666) Holman and (433) Eros against it; then checks the same detections against
the catalog without the snapshot (which integrates it on the fly).
"""
import os
import sys
import tempfile
import time

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from common import detections, kernel
from spacerocks import checker

k = kernel()
path = os.environ.get("MPCORB", os.path.expanduser("~/.spacerocks/mpc/mpcorb_extended.json.gz"))
t = time.time()
cat = checker.Catalog.mpcorb(k, path=path)
print(f"load {cat}: {time.time() - t:.1f} s")
det = pd.concat([detections(x, "2026-09-01", "2026-09-30") for x in ["3666", "433"]]).reset_index(drop=True)
print(len(det), "detections", det.desig.value_counts().to_dict(), "epochs", det.jd_utc.min(), det.jd_utc.max())
args = (det.ra.to_numpy(), det.dec.to_numpy(), det.jd_utc.to_numpy(), det.stn.tolist(), k)
kw = dict(sigma_ra=det.sigma.to_numpy(), sigma_dec=det.sigma.to_numpy())


def report(r, label, dt):
    res = pd.DataFrame({c: r[c] for c in ["detection", "name", "consistent", "distance"]})
    res["truth"] = det.desig.to_numpy()[res.detection]
    first = res.groupby("detection").head(1)
    ok = (first.name == first.truth) & first.consistent
    oth = res[(res.name != res.truth) & res.consistent].assign(sigma_major=lambda x: np.asarray(r["sigma_major"])[x.index])
    print(f"{label}: {dt:.1f} s; right object first and consistent for {ok.sum()}/{len(det)}; "
          f"other consistent pairs {len(oth)} (U {sorted(set(cat.u[np.asarray(r['object'])[oth.index]].tolist()))}, "
          f"median sigma_major {np.median(oth.sigma_major) if len(oth) else float('nan'):.0f}\")")


t = time.time()
r = checker.check(cat, *args, **kw)
report(r, "no snapshot (integrates the catalog)", time.time() - t)

snap = cat.select(np.arange(len(cat)))
t = time.time()
snap.snapshot(float(np.median(det.jd_utc)), k, timescale="utc")
print(f"snapshot of {len(snap)} orbits: {time.time() - t:.1f} s")
f = os.path.join(tempfile.gettempdir(), "mpcorb_snapshot.srcat")
t = time.time(); snap.save(f); ts = time.time() - t
t = time.time(); snap = checker.Catalog.load(f); print(f"save {ts:.1f} s, load {time.time() - t:.1f} s, {os.path.getsize(f) / 1e6:.0f} MB")
t = time.time()
r = checker.check(snap, *args, **kw)
report(r, "with the snapshot", time.time() - t)
