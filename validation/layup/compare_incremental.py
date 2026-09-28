"""Incremental fitting: `orbfit.fit_many(..., prior=...)` vs layup's `incremental_orbitfit`.

    LAYUP_CACHE=/path SPACEROCKS_KERNELS=/path python compare_incremental.py input.csv [n_objects]

From the first n objects (default 40) of an MPC CSV, a "previous" data set and a "current" one:
a quarter of the objects unchanged, a quarter with detections added (the last 10%, at least 3,
missing from the previous set), a quarter with one detection removed since, and a quarter new.
Each code fits the previous set, then brings the fits up to date with the current one. Reported:
the routes each code took, and the agreement of the final orbits.
"""
import os, sys
from collections import Counter
import numpy as np
from layup.utilities.file_io.CSVReader import CSVDataReader
from layup.orbitfit import orbitfit, incremental_orbitfit
from spacerocks import orbfit
from spacerocks.spice import SpiceKernel
from spacerocks.time import Time


def main():
    cache = os.environ["LAYUP_CACHE"]
    K = os.environ["SPACEROCKS_KERNELS"]
    kernel = SpiceKernel()
    for f in ["latest_leapseconds.tls", "de440.bsp", "sb441-n16.bsp", "earth_1962_240827_2124_combined.bpc", "earth_latest_high_prec.bpc"]:
        kernel.load(os.path.join(K, f))

    data = CSVDataReader(sys.argv[1], "csv", primary_id_column_name="provID").read_rows()
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 40
    ids = list(dict.fromkeys(data["provID"]))[:n]
    data = data[np.isin(data["provID"], ids)]
    q = n // 4
    unchanged, appended, removed, new = ids[:q], ids[q:2 * q], ids[2 * q:3 * q], ids[3 * q:]

    keep_prev = np.zeros(len(data), bool)
    keep_cur = np.ones(len(data), bool)
    for oid in ids:
        m = np.where(data["provID"] == oid)[0]
        m = m[np.argsort(data["obsTime"][m].astype(str), kind="mergesort")]
        if oid in unchanged:
            keep_prev[m] = True
        elif oid in appended:
            keep_prev[m[: len(m) - max(3, len(m) // 10)]] = True
        elif oid in removed:
            keep_prev[m] = True
            keep_cur[m[len(m) // 2]] = False
    prev, cur = data[keep_prev], data[keep_cur]

    # layup
    prior_catalog = orbitfit(prev, cache)
    L, routing = incremental_orbitfit(cur, cache, prior_catalog, prior_obs=prev)
    print("layup routes:     ", dict(routing))


    # spacerocks
    def arrays(d):
        ep = [Time.from_isot(t.replace("Z", "")) for t in d["obsTime"].astype(str)]
        return (list(d["provID"].astype(str)), np.radians(d["ra"]), np.radians(d["dec"]), ep, list(d["stn"].astype(str)))


    P = orbfit.fit_many(*arrays(prev), kernel)
    S = orbfit.fit_many(*arrays(cur), kernel, prior=P)
    print("spacerocks routes:", dict(Counter(S["route"])))

    cols = ["x", "y", "z", "xdot", "ydot", "zdot"]
    lay = {str(r["provID"]): r for r in L}
    mismatch, dist = [], []
    for i, oid in enumerate(S["id"]):
        r = lay[oid]
        if int(r["flag"]) != int(S["flag"][i]):
            mismatch.append((oid, int(r["flag"]), int(S["flag"][i])))
            continue
        if int(r["flag"]) != 0:
            continue
        ls = np.array([r[c] for c in cols])
        lc = np.array([[r[f"cov_{a}_{b}"] for b in range(6)] for a in range(6)])
        dx = S["state"][i] - ls
        dist.append((np.sqrt(dx @ np.linalg.solve(lc, dx)), oid, S["route"][i]))
    print(f"flag mismatches: {len(mismatch)} {mismatch}")
    dist.sort(reverse=True)
    print(f"final states (Mahalanobis, layup covariance): median {np.median([x[0] for x in dist]):.2e}, max {dist[0][0]:.2e} ({dist[0][1]}, {dist[0][2]})")


if __name__ == "__main__":  # layup's workers re-import this script
    main()
