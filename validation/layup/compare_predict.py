"""Predictions with uncertainty: `orbfit.predict` vs layup's `predict_sequence`, from the same
orbit and covariance, on layup runs (npz from run_layup.py).

    LAYUP_CACHE=/path SPACEROCKS_KERNELS=/path python compare_predict.py run.npz

Per object: spacerocks fits the detections (from layup's orbit), and both codes predict that
orbit, with its covariance, at every detection (its observer) and at 12 geocentric epochs up to
three years past the last detection. Reported: the angle between the predicted directions and
the relative difference of the 2x2 sky covariances.
"""
import os, sys
import numpy as np
from layup.routines import FitResult, Observation, get_ephem, numpy_to_eigen, predict_sequence
from spacerocks import SpaceRock, orbfit
from spacerocks.spice import SpiceKernel
from spacerocks.time import Time

K = os.environ["SPACEROCKS_KERNELS"]
kernel = SpiceKernel()
for f in ["latest_leapseconds.tls", "de440.bsp", "sb441-n16.bsp", "earth_1962_240827_2124_combined.bpc", "earth_latest_high_prec.bpc"]:
    kernel.load(os.path.join(K, f))
ephem = get_ephem(os.environ["LAYUP_CACHE"])
SIG = 1.0 / 206265.0

d = np.load(sys.argv[1], allow_pickle=True)
obs, et, pv, fits = d["obs"], d["et"], d["pv"], d["fits"]
ids = obs["provID"].astype(str)
sep, dcov, big = [], [], []
for row in fits:
    oid = str(row["provID"])
    if int(row["flag"]) != 0:
        continue
    m = np.where(ids == oid)[0]
    jd = 2451545.0 + et[m] / 86400.0
    ra, dec = np.radians(obs["ra"][m]), np.radians(obs["dec"][m])
    pos = np.column_stack([pv["x"][m], pv["y"][m], pv["z"][m]])
    rock = SpaceRock.from_xyz(oid, row["x"], row["y"], row["z"], row["xdot"], row["ydot"], row["zdot"],
                              Time(row["epochMJD_TDB"] + 2400000.5, "tdb", "jd"), "J2000", "SSB")
    S = orbfit.fit(ra, dec, jd, pos, kernel, sigma_ra=SIG, sigma_dec=SIG, timescale="tdb", initial=rock)
    if S.flag != 0:
        continue
    future = jd.max() + np.linspace(10, 3 * 365.25, 12)
    geo = orbfit.observers(["500"] * len(future), list(future), kernel, timescale="tdb")[:, :3]
    t_all = np.concatenate([jd, future])
    p_all = np.vstack([pos, geo])

    P = orbfit.predict(S, t_all, p_all, kernel, timescale="tdb")
    fr = FitResult()
    fr.state, fr.epoch = list(S.state), S.epoch.jd()
    O = []
    for t, p in zip(t_all, p_all):
        o = Observation()
        o.observer_position, o.observer_velocity, o.epoch = list(p), [0.0, 0.0, 0.0], float(t)
        O.append(o)
    L = predict_sequence(ephem, fr, O, numpy_to_eigen(list(np.asarray(S.state_covariance).ravel()), 6, 6))
    for i, l in enumerate(L):
        u = np.array([np.cos(P["dec"][i]) * np.cos(P["ra"][i]), np.cos(P["dec"][i]) * np.sin(P["ra"][i]), np.sin(P["dec"][i])])
        lr = np.array(l.rho)
        sep.append(np.degrees(np.arctan2(np.linalg.norm(np.cross(u, lr)), u @ lr)) * 3.6e6)
        lc = np.array(l.obs_cov).reshape(2, 2)
        dcov.append(np.linalg.norm(P["covariance"][i] - lc) / np.linalg.norm(lc))
        big.append(P["sigma_major"][i])
sep, dcov, big = map(np.array, (sep, dcov, big))
print(f"{len(sep)} predictions: direction differs by median {np.median(sep):.2e} mas, max {sep.max():.2e} mas; "
      f"covariance relative difference median {np.median(dcov):.1e}, max {dcov.max():.1e}; "
      f"1-sigma major axes from {big.min():.2g}\" to {big.max():.2g}\"")
