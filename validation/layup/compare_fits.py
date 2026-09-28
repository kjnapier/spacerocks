"""Compare spacerocks fits with layup's: python compare_fits.py layup_fits.csv sr_fits.csv"""
import sys
import numpy as np

def load(p):
    rows = [l.rstrip("\n").split(",") for l in open(p)]
    head = rows[0]
    return {r[0]: dict(zip(head[1:], map(float, r[1:]))) for r in rows[1:]}

L, S = load(sys.argv[1]), load(sys.argv[2])
cols = ["x", "y", "z", "xdot", "ydot", "zdot"]
worst = []
nflag = 0
for oid, l in L.items():
    s = S[oid]
    if int(l["flag"]) != int(s["flag"]):
        nflag += 1
        print(f"{oid:>10} FLAG layup {int(l['flag'])} spacerocks {int(s['flag'])}  csq {l['csq']:.6g} / {s['csq']:.6g}")
        continue
    if int(l["flag"]) != 0:
        continue
    dx = np.array([s[c] - l[c] for c in cols])
    cov = np.array([[l[f"cov_{i}_{j}"] for j in range(6)] for i in range(6)])
    sig = np.sqrt(np.diag(cov))
    maha = float(np.sqrt(dx @ np.linalg.solve(cov, dx)))
    worst.append((maha, oid, np.abs(dx[:3]).max(), abs(s["csq"] - l["csq"]) / l["csq"], s["epoch_tdb"] - l["epoch_tdb"], int(l["niter"]), int(s["niter"]), np.abs(dx / sig).max()))
worst.sort(reverse=True)
print(f"{len(L)} objects, {nflag} flag mismatches, {len(worst)} both converged")
print("   id     mahalanobis   max|dpos| AU   dcsq/csq   depoch     niter L/S  max|dx|/sigma")
for w in worst[:12]:
    print(f"{w[1]:>10} {w[0]:12.3e} {w[2]:12.3e} {w[3]:10.2e} {w[4]:9.2e}   {w[5]:3d}/{w[6]:3d}   {w[7]:.2e}")
m = np.array([w[0] for w in worst])
print("mahalanobis distance: median %.2e  max %.2e" % (np.median(m), m.max()))
print("niter equal: %d / %d" % (sum(w[5] == w[6] for w in worst), len(worst)))
if "seconds" in next(iter(S.values())):
    print("spacerocks seconds per object: median %.3f  total %.1f" % (np.median([s["seconds"] for s in S.values()]), sum(s["seconds"] for s in S.values())))
