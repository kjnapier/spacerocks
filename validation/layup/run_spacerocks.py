"""Fit the same CSV of astrometry with spacerocks, end to end (time parsing and observatory
positions included), and compare with a layup run (run_layup.py).

    SPACEROCKS_KERNELS=/path/to/kernels python run_spacerocks.py input.csv layup_run.npz

Set CONV_FRAC to use the scaled convergence test, and WEIGHT_DATA=veres for Vereš et al. (2017)
weights, DUMP=out.json for the per-object differences, IOD=gauss or IOD=herget for the other initial-orbit methods, ENGINE=bk_native for the
Bernstein-Khushalani engine (compare with a layup run
made with the same).
"""
import os, sys, time
import numpy as np

from spacerocks.spice import SpiceKernel
from spacerocks.time import Time
from spacerocks import orbfit

K = os.environ["SPACEROCKS_KERNELS"]
kernel = SpiceKernel()
for f in ["latest_leapseconds.tls", "de440.bsp", "sb441-n16.bsp", "earth_1962_240827_2124_combined.bpc", "earth_latest_high_prec.bpc"]:
    kernel.load(os.path.join(K, f))

ref = np.load(sys.argv[2], allow_pickle=True)
obs, fits = ref["obs"], ref["fits"]
ids = obs["provID"].astype(str)
epochs = [Time.from_isot(t.replace("Z", "")) for t in obs["obsTime"].astype(str)]

# time conversion and observer positions against layup's (SPICE + Sorcha)
jd_layup = 2451545.0 + ref["et"] / 86400.0
jd_sr = np.array([t.tdb().jd() for t in epochs])
print("TDB epoch difference: max %.3g us" % (np.abs(jd_sr - jd_layup).max() * 86400e6))

t0 = time.time()
sigma = None
if os.environ.get("WEIGHT_DATA", "") == "veres":
    sigma = orbfit.veres_sigma(list(obs["stn"].astype(str)), epochs)
res = orbfit.fit_many(ids, obs["ra"] * np.pi / 180.0, obs["dec"] * np.pi / 180.0, epochs, list(obs["stn"].astype(str)), kernel,
                      sigma_ra=sigma, sigma_dec=sigma, conv_frac=float(os.environ.get("CONV_FRAC", "0")),
                      iod=os.environ.get("IOD", "auto"), engine=os.environ.get("ENGINE", "cartesian"))
print("fit_many: %d objects in %.2f s" % (len(res["id"]), time.time() - t0))

cols = ["x", "y", "z", "xdot", "ydot", "zdot"]
lay = {str(r["provID"]): r for r in fits}
nflag, maha, dchi, per = 0, [], [], []
for i, oid in enumerate(res["id"]):
    l = lay[oid]
    if int(l["flag"]) != int(res["flag"][i]):
        nflag += 1
        print("flag mismatch", oid, l["flag"], res["flag"][i])
        continue
    if l["flag"] != 0:
        continue
    cov = np.array([[l[f"cov_{a}_{b}"] for b in range(6)] for a in range(6)])
    dx = res["state"][i] - np.array([l[c] for c in cols])
    maha.append(np.sqrt(dx @ np.linalg.solve(cov, dx)))
    dchi.append(abs(res["chi2"][i] - l["csq"]) / l["csq"])
    dpos = np.linalg.norm(dx[:3]) * 149597870.7
    per.append(dict(id=oid, sigma=float(maha[-1]), dchi2=float(dchi[-1]), dpos_km=float(dpos),
                    n=int((ids == oid).sum()), chi2_ndof=float(l["csq"] / l["ndof"])))
print("flag mismatches: %d; converged in both: %d" % (nflag, len(maha)))
print("state difference (Mahalanobis, layup covariance): median %.2e max %.2e" % (np.median(maha), np.max(maha)))
print("chi2 relative difference: median %.2e max %.2e" % (np.median(dchi), np.max(dchi)))
if os.environ.get("DUMP"):
    import json
    json.dump(dict(per_object=per, flag_mismatches=nflag, time_fit_many=None), open(os.environ["DUMP"], "w"), indent=1)
