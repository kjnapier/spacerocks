"""Real MPC detections of five objects checked against MPCORB.

Every ground-based detection since 2000 of (3666) Holman, (433) Eros, (90377) Sedna,
(99942) Apophis and (101955) Bennu (a sample of at most 400 per object) is checked against
the five orbits plus 30,000 other MPCORB orbits drawn at random (the full catalog would need a
full-catalog integration per 20-day cluster of epochs; `validate_full.py` times that).
Reports how often the right object is found, is consistent, and ranks first, the distribution
of its Mahalanobis distance, and the other objects found consistent.

Needs: the MPC ADES files cached by the orbit-fitting notebook (MPC_CACHE), the MPCORB JSON
(MPCORB, default ~/.spacerocks/mpc/mpcorb_extended.json.gz) and kernels (SPACEROCKS_KERNELS).
"""
import os
import sys
import time

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from common import ARCSEC, detections, kernel
from spacerocks import checker

OBJECTS = {"3666": "(3666) Holman", "433": "(433) Eros", "90377": "(90377) Sedna", "99942": "(99942) Apophis", "101955": "(101955) Bennu"}

k = kernel()
path = os.environ.get("MPCORB", os.path.expanduser("~/.spacerocks/mpc/mpcorb_extended.json.gz"))
full = checker.Catalog.mpcorb(k, path=path)
rng = np.random.default_rng(1)
idx = [full.index(d) for d in OBJECTS]
others = rng.choice(len(full), 30000, replace=False)
cat = full.select(np.unique(np.concatenate([idx, others])))
print(cat)

det = []
for d in OBJECTS:
    x = detections(d, "2000-01-01")
    if len(x) > 400:
        x = x.sample(400, random_state=1)
    det.append(x)
det = pd.concat(det).sort_values("jd_utc").reset_index(drop=True)
print(len(det), "detections", det.desig.value_counts().to_dict())

t = time.time()
r = checker.check(cat, det.ra.to_numpy(), det.dec.to_numpy(), det.jd_utc.to_numpy(), det.stn.tolist(), k,
                  sigma_ra=det.sigma.to_numpy(), sigma_dec=det.sigma.to_numpy())
print(f"check: {time.time() - t:.1f} s")
cols = ["detection", "name", "distance", "consistent", "log_likelihood", "separation", "sigma_major", "dra", "ddec"]
res = pd.DataFrame({c: r[c] for c in cols})
res["truth"] = det.desig.to_numpy()[res.detection]
res["rank"] = res.groupby("detection").cumcount()
true = res[res.name == res.truth].set_index("detection")

rows = []
for d, label in OBJECTS.items():
    js = det.index[det.desig == d]
    t_ = true.reindex(js)
    found = t_.distance.notna()
    rows.append({
        "object": label, "n": len(js),
        "found": found.sum(), "consistent": (t_.consistent == True).sum(),
        "ranked first": (t_["rank"] == 0).sum(),
        "median distance": t_.distance.median(), "95% distance": t_.distance.quantile(0.95),
        "median sep (\")": t_.separation.median(), "median sigma_det (\")": det.sigma[js].median() / ARCSEC,
    })
summary = pd.DataFrame(rows)
print(summary.to_string(index=False, float_format=lambda x: f"{x:.2f}"))
oth = res[(res.name != res.truth) & res.consistent]
print(f"other objects consistent with a detection: {len(oth)} pairs "
      f"({oth.detection.nunique()} detections); median sigma_major {oth.sigma_major.median():.0f}\"")
print("worst true-object pairs:")
print(true.sort_values("distance").tail(8)[["truth", "distance", "separation", "sigma_major", "dra", "ddec"]].join(det[["jd_utc", "stn", "sigma"]]).assign(sigma=lambda x: x.sigma / ARCSEC))
summary.to_csv(os.path.join(os.path.dirname(__file__), "validate_mpcorb.csv"), index=False)
