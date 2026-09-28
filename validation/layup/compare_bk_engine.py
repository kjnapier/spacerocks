"""The Bernstein-Khushalani fitting engine alone: `fit(..., initial=seed, engine="bk_native")` vs
layup's `run_bk_native_fit`, from the same seed and inputs, on layup runs (npz from
run_layup.py). The seed is layup's Bernstein-Khushalani IOD (`run_bk_iod`) over all of an
object's detections at the middle one, the pipeline's fallback seed.

    LAYUP_CACHE=/path SPACEROCKS_KERNELS=/path python compare_bk_engine.py run.npz [run2.npz ...]

(An end-to-end run of layup with ENGINE=bk_native on the synthetic short arcs stalls: IAS15
grinds on one main-belt arc for over half an hour, so the engine is compared directly.)
layup reports a fit that does not converge with flag 2 (its chi-square stays infinite); that is
read as 1 here. spacerocks' flag 9 (converged, but a hyperbolic excess speed above 200 km/s, a
check of the pipeline that layup's engine doesn't make) is counted as converged and listed.
"""
import os, sys
import numpy as np
from layup.routines import FitResult, Observation, get_ephem, run_bk_iod, run_bk_native_fit
from layup.constants import MU_SUN
from spacerocks import SpaceRock, orbfit
from spacerocks.spice import SpiceKernel
from spacerocks.time import Time

K = os.environ["SPACEROCKS_KERNELS"]
kernel = SpiceKernel()
for f in ["latest_leapseconds.tls", "de440.bsp", "sb441-n16.bsp"]:
    kernel.load(os.path.join(K, f))
ephem = get_ephem(os.environ["LAYUP_CACHE"])
SIG = 1.0 / 206265.0

rows = []
for path in sys.argv[1:]:
    d = np.load(path, allow_pickle=True)
    obs, et, pv = d["obs"], d["et"], d["pv"]
    ids = obs["provID"].astype(str)
    for oid in dict.fromkeys(ids):
        m = np.where(ids == oid)[0]
        m = m[np.argsort(et[m], kind="mergesort")]
        if len(m) < 4:
            continue
        jd = 2451545.0 + et[m] / 86400.0
        ra, dec = np.radians(obs["ra"][m]), np.radians(obs["dec"][m])
        pos = np.column_stack([pv["x"][m], pv["y"][m], pv["z"][m]])
        O = []
        for i in range(len(m)):
            o = Observation.from_astrometry_with_id(oid, ra[i], dec[i], jd[i], list(pos[i]), [0, 0, 0])
            o.ra_unc = o.dec_unc = SIG
            O.append(o)
        e0 = float(jd[len(jd) // 2])
        seed = run_bk_iod(O, e0, MU_SUN)
        if seed.flag != 0:
            continue
        L = run_bk_native_fit(ephem, seed, O, MU_SUN)
        rock = SpaceRock.from_xyz(oid, *seed.state, Time(seed.epoch, "tdb", "jd"), "J2000", "SSB")
        S = orbfit.fit(ra, dec, jd, pos, kernel, sigma_ra=SIG, sigma_dec=SIG, timescale="tdb", initial=rock, engine="bk_native")
        sflag = 0 if S.flag == 9 else S.flag
        lflag = 1 if (L.flag == 2 and not np.isfinite(L.csq)) else L.flag
        cov = np.array(L.cov).reshape(6, 6)
        dx = S.state - np.array(L.state)
        dm = np.sqrt(dx @ np.linalg.solve(cov, dx)) if lflag in (0, 2) else np.nan
        rows.append((path, oid, lflag, sflag, dm, abs(S.chi2 / L.csq - 1) if lflag in (0, 2) else np.nan, L.niter, S.niter, S.flag == 9))

for path in sys.argv[1:]:
    r = [x for x in rows if x[0] == path]
    both = [x for x in r if x[2] in (0, 2) and x[2] == x[3]]
    mism = [x for x in r if x[2] != x[3]]
    dm = np.array([x[4] for x in both])
    print(f"{os.path.basename(path)}: {len(r)} objects; flag mismatches {len(mism)}; converged in both {sum(x[2] == 0 for x in both)}; "
          f"same iterations {sum(x[6] == x[7] for x in both)}/{len(both)}")
    if len(both):
        print(f"   states (Mahalanobis, layup covariance): median {np.median(dm):.2e}, max {dm.max():.2e}; "
              f"chi2 relative max {np.nanmax([x[5] for x in both]):.2e}")
    print(f"   flagged implausible (9) by spacerocks: {[x[1] for x in r if x[8]]}")
    for x in mism[:8]:
        print(f"   mismatch {x[1]}: layup {x[2]} ({x[6]} iterations), spacerocks {x[3]} ({x[7]})")
