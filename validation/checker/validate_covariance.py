"""The checker with orbits that carry a covariance: fits made here, and the MPC's own orbit.

- (3666) Holman fitted to its 2015-2023 detections, then checked against its detections from
  2024 on (prediction) and from 1990-2014 (retrodiction), among 30,000 MPCORB orbits.
- Short arcs: Holman fitted to 30 and to 7 days of 2023 data; the covariance is large and
  elongated, and the later detections should still be consistent with it.
- The MPC's orbit of Holman (mpc_orb JSON, with its covariance) against the same detections.
- Negative control: every detection moved 10" is not consistent with the long-arc fit.
"""
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from common import ARCSEC, CACHE, detections, kernel
from spacerocks import checker, orbfit

k = kernel()
path = os.environ.get("MPCORB", os.path.expanduser("~/.spacerocks/mpc/mpcorb_extended.json.gz"))
full = checker.Catalog.mpcorb(k, path=path)
rng = np.random.default_rng(2)
others = full.select(rng.choice(len(full), 30000, replace=False))
others = others.select(np.array([n != "3666" for n in others.names]))

d = detections("3666", "1990-01-01")


def fit(sel, name):
    x = d[sel]
    f = orbfit.fit(x.ra.to_numpy(), x.dec.to_numpy(), x.jd_utc.to_numpy(), x.stn.tolist(), k,
                   sigma_ra=x.sigma.to_numpy(), sigma_dec=x.sigma.to_numpy(), robust=True, name=name)
    s = np.sqrt(np.diag(f.state_covariance)[:3]) * 1.496e8
    print(f"{name}: {len(x)} detections, flag {f.flag}, chi2/ndof {f.chi2 / f.ndof:.2f}, position sigma {s.round(1)} km")
    return f


yr = lambda s: (d.jd_utc - 2451545.0) / 365.25 + 2000.0
long_fit = fit((yr(d) >= 2015) & (yr(d) < 2024), "fit 2015-2023")
t0 = d.jd_utc[(yr(d) >= 2023) & (yr(d) < 2024)].min()
f30 = fit((d.jd_utc >= t0) & (d.jd_utc < t0 + 30), "fit 30 d")
f7 = fit((d.jd_utc >= t0) & (d.jd_utc < t0 + 7), "fit 7 d")
mpc_orb = json.loads((CACHE / "orb_3666.json").read_text())


def run(cat, x, label, **kw):
    r = checker.check(cat, x.ra.to_numpy(), x.dec.to_numpy(), x.jd_utc.to_numpy(), x.stn.tolist(), k,
                      sigma_ra=x.sigma.to_numpy(), sigma_dec=x.sigma.to_numpy(), **kw)
    res = pd.DataFrame({c: r[c] for c in ["detection", "name", "distance", "consistent", "sigma_major", "separation"]})
    res["rank"] = res.groupby("detection").cumcount()
    t = res[res.name == "3666"].set_index("detection").reindex(range(len(x)))
    oth = res[(res.name != "3666") & res.consistent]
    row = {"case": label, "n": len(x), "found": int(t.distance.notna().sum()), "consistent": int((t.consistent == True).sum()),
           "first": int((t["rank"] == 0).sum()), "median d": t.distance.median(), "95% d": t.distance.quantile(0.95),
           "median sigma_major (\")": t.sigma_major.median(), "others consistent": len(oth)}
    return row


rows = []
later = d[yr(d) >= 2024].reset_index(drop=True)
earlier = d[(yr(d) >= 1990) & (yr(d) < 2015)].sample(300, random_state=1).reset_index(drop=True)
after23 = d[d.jd_utc >= t0 + 30].sample(300, random_state=1).reset_index(drop=True)
# The short-arc fits' uncertainties grow past the default max_uncertainty (600") within a
# year; they are also run with it lifted.
for f, x, label, kw in [(long_fit, later, "fit 2015-23 -> 2024-26", {}), (long_fit, earlier, "fit 2015-23 -> 1990-2014", {}),
                        (f30, after23, "fit 30 d -> later", {}), (f30, after23, "fit 30 d -> later, max_unc 10 deg", {"max_uncertainty": 36000}),
                        (f7, after23, "fit 7 d -> later", {})]:
    cat = checker.Catalog.from_fits({"3666": f})
    cat.extend(others)
    rows.append(run(cat, x, label, **kw))
cat = checker.Catalog.from_mpc_orb(mpc_orb, k)
print(cat, cat.names)
cat.extend(others)
rows.append(run(cat, later, "MPC orbit -> 2024-26"))
rows.append(run(cat, earlier, "MPC orbit -> 1990-2014"))

# Negative control: move every detection 10" (in a random direction).
cat = checker.Catalog.from_fits({"3666": long_fit})
moved = later.copy()
phi = rng.uniform(0, 2 * np.pi, len(moved))
moved["dec"] = moved.dec + 10 * ARCSEC * np.sin(phi)
moved["ra"] = moved.ra + 10 * ARCSEC * np.cos(phi) / np.cos(moved.dec)
rows.append(run(cat, moved, "moved 10\" (should fail)"))

out = pd.DataFrame(rows)
print(out.to_string(index=False, float_format=lambda v: f"{v:.2f}"))
out.to_csv(os.path.join(os.path.dirname(__file__), "validate_covariance.csv"), index=False)
