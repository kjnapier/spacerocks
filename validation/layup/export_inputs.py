"""Write the detections layup fitted (sorted per object, TDB epochs and layup's observer
positions) and layup's results as CSV files for the Rust comparison driver.

    python export_inputs.py layup_run.npz outdir
"""
import sys, os
import numpy as np

d = np.load(sys.argv[1], allow_pickle=True)
out = sys.argv[2]
obs, et, pv, fits = d["obs"], d["et"], d["pv"], d["fits"]
jd = 2451545.0 + et / 86400.0
ids = obs["provID"].astype(str)
with open(os.path.join(out, "detections.csv"), "w") as f:
    f.write("id,epoch_tdb,ra,dec,ox,oy,oz,stn,obstime\n")
    for oid in dict.fromkeys(ids):
        m = np.where(ids == oid)[0]
        # layup sorts each object's rows by obsTime (stable)
        m = m[np.argsort(obs["obsTime"][m].astype(str), kind="mergesort")]
        for i in m:
            f.write(f"{oid},{float(jd[i])!r},{float(obs['ra'][i] * np.pi / 180.0)!r},{float(obs['dec'][i] * np.pi / 180.0)!r},"
                    f"{float(pv['x'][i])!r},{float(pv['y'][i])!r},{float(pv['z'][i])!r},{obs['stn'][i]},{obs['obsTime'][i]}\n")
cols = ["x", "y", "z", "xdot", "ydot", "zdot"]
with open(os.path.join(out, "layup_fits.csv"), "w") as f:
    f.write("id,flag,csq,ndof,niter,epoch_tdb," + ",".join(cols) + "," + ",".join(f"cov_{i}_{j}" for i in range(6) for j in range(6)) + "\n")
    for r in fits:
        vals = [r["flag"], r["csq"], r["ndof"], r["niter"], r["epochMJD_TDB"] + 2400000.5] + [r[c] for c in cols] + [r[f"cov_{i}_{j}"] for i in range(6) for j in range(6)]
        f.write(str(r["provID"]) + "," + ",".join(repr(float(v)) for v in vals) + "\n")
