"""Sequential updates: spacerocks' `orbfit.sequential_update` vs layup's `run_sequential_update`,
on the 99 objects of a layup run (npz from run_layup.py).

    LAYUP_CACHE=/path SPACEROCKS_KERNELS=/path python compare_sequential.py run.npz

Per object, the detections are split in time: the last 10% (at least 3) are "new". Each code
fits the old ones from layup's full-arc orbit (the prior), updates that prior with the new ones
alone, and refits everything from the prior. Reported: the update's agreement between the codes
(Mahalanobis distance in layup's posterior covariance, chi-square, and the size of the move that
the max_update_sigma gate tests), and how close each code's update comes to its own full refit.
"""
import os, sys
import numpy as np
import layup.orbitfit as LO
from layup.routines import Observation, get_ephem, run_from_vector_with_initial_guess, run_sequential_update
from layup.utilities.data_processing_utilities import parse_fit_result
from spacerocks import SpaceRock, orbfit
from spacerocks.spice import SpiceKernel
from spacerocks.time import Time

K = os.environ["SPACEROCKS_KERNELS"]
kernel = SpiceKernel()
for f in ["latest_leapseconds.tls", "de440.bsp", "sb441-n16.bsp"]:
    kernel.load(os.path.join(K, f))
ephem = get_ephem(os.environ["LAYUP_CACHE"])
SIG = 1.0 / 206265.0


def mahal(a, b, cov):
    d = np.asarray(a) - np.asarray(b)
    return float(np.sqrt(d @ np.linalg.solve(np.asarray(cov).reshape(6, 6), d)))


d = np.load(sys.argv[1], allow_pickle=True)
obs, et, pv, fits = d["obs"], d["et"], d["pv"], d["fits"]
ids = obs["provID"].astype(str)
rows = []
for row in fits:
    oid = str(row["provID"])
    if int(row["flag"]) != 0:
        continue
    m = np.where(ids == oid)[0]
    m = m[np.argsort(et[m], kind="mergesort")]
    n_new = max(3, len(m) // 10)
    if len(m) - n_new < 6:
        continue
    jd = 2451545.0 + et[m] / 86400.0
    ra, dec = np.radians(obs["ra"][m]), np.radians(obs["dec"][m])
    pos = np.column_stack([pv["x"][m], pv["y"][m], pv["z"][m]])
    O = []
    for i in range(len(m)):
        o = Observation.from_astrometry_with_id(oid, ra[i], dec[i], jd[i], list(pos[i]), [0, 0, 0])
        o.ra_unc = o.dec_unc = SIG
        O.append(o)
    old, new = slice(0, len(m) - n_new), slice(len(m) - n_new, len(m))

    guess = parse_fit_result(row)
    Lp = run_from_vector_with_initial_guess(ephem, guess, O[old])
    if Lp.flag != 0:
        continue
    Ls = run_sequential_update(ephem, Lp, O[new])
    Lf = run_from_vector_with_initial_guess(ephem, Lp, O)

    rock = SpaceRock.from_xyz(oid, *guess.state, Time(guess.epoch, "tdb", "jd"), "J2000", "SSB")
    args = dict(sigma_ra=SIG, sigma_dec=SIG, timescale="tdb")
    Sp = orbfit.fit(ra[old], dec[old], jd[old], pos[old], kernel, initial=rock, **args)
    if Sp.flag != 0:
        rows.append((oid, Lp.flag, Sp.flag, None))
        continue
    Ss = orbfit.sequential_update(Sp, ra[new], dec[new], jd[new], pos[new], kernel, **args)
    Sf = orbfit.fit(ra, dec, jd, pos, kernel, initial=Sp.rock, **args)
    Sr = orbfit.fit(ra, dec, jd, pos, kernel, prior=Sp, **args)
    rows.append((oid, Ls.flag, Ss.flag, dict(
        d_update=mahal(Ss.state, Ls.state, Ls.cov) if Ls.flag == 0 else np.nan,
        chi2=(Ss.chi2, Ls.csq), ndof=(Ss.ndof, Ls.ndof),
        gate=(mahal(Ss.state, Sp.state, Sp.state_covariance), mahal(Ls.state, Lp.state, Lp.cov)),
        vs_full=(mahal(Ss.state, Sf.state, Sf.state_covariance), mahal(Ls.state, Lf.state, Lf.cov) if Lf.flag == 0 else np.nan),
        route=Sr.route)))

ok = [r for r in rows if r[3] is not None]
print(f"{len(rows)} objects with a converged prior in layup; spacerocks prior converged for {len(ok)}")
print(f"flag mismatches in the update: {sum(r[1] != r[2] for r in ok)}; converged in both: {sum(r[1] == 0 and r[2] == 0 for r in ok)}")
both = [r for r in ok if r[1] == 0 and r[2] == 0]
du = np.array([r[3]["d_update"] for r in both])
chi = np.array([abs(r[3]["chi2"][0] / r[3]["chi2"][1] - 1) for r in both])
g = np.array([r[3]["gate"] for r in both])
vf = np.array([r[3]["vs_full"] for r in both])
print(f"updated state, spacerocks vs layup: median {np.median(du):.2e} max {du.max():.2e} sigma; chi2 relative max {chi.max():.2e}; "
      f"ndof equal: {all(r[3]['ndof'][0] == r[3]['ndof'][1] for r in both)}")
print(f"move tested by the gate (prior sigmas): spacerocks vs layup differ by at most {np.abs(g[:, 0] - g[:, 1]).max():.2e}; "
      f"{np.sum(g[:, 1] > 4)} of {len(g)} exceed 4 in layup, {np.sum(g[:, 0] > 4)} in spacerocks")
print(f"update vs its own full refit (sigma): spacerocks median {np.median(vf[:, 0]):.2e} max {vf[:, 0].max():.2e}; "
      f"layup median {np.nanmedian(vf[:, 1]):.2e} max {np.nanmax(vf[:, 1]):.2e}")
from collections import Counter
print("fit(..., prior=) routes:", dict(Counter(r[3]["route"] for r in ok)))
for r in rows:
    if r[3] is None or r[1] != r[2]:
        print("  ", r[0], "layup", r[1], "spacerocks", r[2])
