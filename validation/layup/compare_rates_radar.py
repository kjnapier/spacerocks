"""Streak (rate) and radar fits: layup vs spacerocks.

    LAYUP_TESTS=/path/to/layup/tests LAYUP_CACHE=/path SPACEROCKS_KERNELS=/path python compare_rates_radar.py

Cases:
  1. layup's synthetic streak arc (tests/data/streak_synthetic.json), and with one rate
     corrupted by 50 sigma.
  2. layup's synthetic radar arc (tests/data/radar_synthetic.json): delay + Doppler, delay
     only, Doppler only.
  3. Real JPL radar of (99942) Apophis, 2013 (tests/data/apophis_2013_radar.json): with layup's
     own station states and transmitter model, and end to end with spacerocks' observatories and
     its exact transmitter state (the two differ by design; see README).
  4. layup's streak-format ADES file (tests/data/streak_observations_complete.csv) through both
     full pipelines (IOD included).

For each, the model is compared by chi-square at fixed states (a 1-iteration fit returns the
chi-square at its starting state), then the fits from the same perturbed start.
"""
import csv, json, os
import numpy as np

import layup.orbitfit as LO
from layup.orbitfit import orbitfit, _get_result_dtypes, _radar_observation, _append_observer_acceleration
from layup.routines import FitResult, Observation, get_ephem, run_from_vector_with_initial_guess
from layup.utilities.data_processing_utilities import LayupObservatory, get_cov_columns
from layup.utilities.file_io.CSVReader import CSVDataReader
import numpy.lib.recfunctions as rfn
import spiceypy as spice

from spacerocks import SpaceRock, orbfit
from spacerocks.spice import SpiceKernel
from spacerocks.time import Time

TESTS = os.environ["LAYUP_TESTS"]
CACHE = os.environ["LAYUP_CACHE"]
EPHEM = get_ephem(CACHE)
kernel = SpiceKernel()
for f in ["latest_leapseconds.tls", "de440.bsp", "sb441-n16.bsp", "earth_1962_240827_2124_combined.bpc", "earth_latest_high_prec.bpc"]:
    kernel.load(os.path.join(os.environ["SPACEROCKS_KERNELS"], f))
C = 2.99792458e8 * 86400.0 / 149597870700.0
F = 2.38e9  # an arbitrary transmit frequency for the synthetic Doppler (Hz)


def seed(state, epoch):
    g = FitResult()
    g.state, g.epoch, g.flag = list(state), epoch, 0
    return g


def rock(state, epoch):
    return SpaceRock.from_xyz("x", *state, Time(epoch, "tdb", "jd"), "J2000", "SSB")


def report(name, obs_layup, sr_kwargs, states, start, epoch, iters=100):
    """chi-square at each of `states` in both codes, then both fits from `start`."""
    print(f"\n== {name}")
    for label, s in states:
        L = run_from_vector_with_initial_guess(EPHEM, seed(s, epoch), obs_layup, 1)
        S = orbfit.fit(**sr_kwargs, initial=rock(s, epoch), max_iter=1)
        print(f"   chi2 at {label:10}: layup {L.csq:.10e}  spacerocks {S.chi2:.10e}  rel diff {abs(S.chi2 - L.csq) / max(L.csq, 1e-300):.1e}")
    L = run_from_vector_with_initial_guess(EPHEM, seed(start, epoch), obs_layup, iters)
    S = orbfit.fit(**sr_kwargs, initial=rock(start, epoch), max_iter=iters)
    print(f"   fit: layup flag {L.flag} niter {L.niter} ndof {L.ndof} chi2 {L.csq:.6e} | spacerocks flag {S.flag} niter {S.niter} ndof {S.ndof} chi2 {S.chi2:.6e}")
    if L.flag == 0 and S.flag == 0:
        cov = np.array(L.cov).reshape(6, 6)
        dx = S.state - np.array(L.state)
        print(f"   state difference {np.sqrt(dx @ np.linalg.solve(cov, dx)):.2e} sigma, {np.abs(dx[:3]).max():.2e} AU")
    return L, S


def main():
    # --- 1. synthetic streaks ----------------------------------------------------------------------
    d = json.load(open(os.path.join(TESTS, "data/streak_synthetic.json")))
    truth, epoch = np.array(d["true_state"]), d["epoch"]
    O = d["observations"]
    for corrupt in (False, True):
        obs_l, rr = [], []
        for k, o in enumerate(O):
            r = o["ra_rate"] + (50 * d["rate_unc_radday"] if corrupt and k == 0 else 0.0)
            rr.append(r)
            ob = Observation.from_streak_with_id("s", o["ra"], o["dec"], r, o["dec_rate"], o["epoch"], o["observer_position"],
                                                 o["observer_velocity"], d["rate_unc_radday"], d["rate_unc_radday"])
            ob.ra_unc = ob.dec_unc = d["ra_unc_rad"]
            obs_l.append(ob)
        kw = dict(ra=[o["ra"] for o in O], dec=[o["dec"] for o in O], epoch=[o["epoch"] for o in O],
                  observer=np.array([o["observer_position"] + o["observer_velocity"] for o in O]), kernel=kernel, timescale="tdb",
                  sigma_ra=d["ra_unc_rad"], sigma_dec=d["ra_unc_rad"], ra_rate=rr, dec_rate=[o["dec_rate"] for o in O],
                  sigma_ra_rate=d["rate_unc_radday"], sigma_dec_rate=d["rate_unc_radday"])
        start = truth.copy()
        start[3:] *= 1.001
        report("synthetic streaks" + (" (one rate corrupted by 50 sigma)" if corrupt else ""), obs_l, kw,
               [("truth", truth), ("start", start)], start, epoch, 50)

    # --- 2. synthetic radar --------------------------------------------------------------------------
    d = json.load(open(os.path.join(TESTS, "data/radar_synthetic.json")))
    truth, epoch = np.array(d["true_state"]), d["epoch"]
    O = d["observations"]
    for has_delay, has_doppler in ((True, True), (True, False), (False, True)):
        obs_l = [Observation.from_radar_with_id("r", o["delay"], o["doppler"], has_delay, has_doppler, o["epoch"], o["observer_position"],
                                                o["observer_velocity"], d["delay_unc_days"], d["doppler_unc_audy"]) for o in O]
        nan = np.full(len(O), np.nan)
        kw = dict(ra=nan, dec=nan, epoch=[o["epoch"] for o in O], observer=np.array([o["observer_position"] + o["observer_velocity"] for o in O]),
                  kernel=kernel, timescale="tdb",
                  delay=np.array([o["delay"] * 86400.0 for o in O]) if has_delay else None,
                  doppler=np.array([-o["doppler"] * F / C for o in O]) if has_doppler else None,
                  sigma_delay=d["delay_unc_days"] * 86400.0, sigma_doppler=d["doppler_unc_audy"] * F / C, frequency=F)
        start = truth.copy()
        start[:3] += 1e-6
        start[3:] += 1e-8
        report(f"synthetic radar (delay={has_delay}, doppler={has_doppler})", obs_l, kw, [("truth", truth), ("start", start)], start, epoch)

    # --- 3. Apophis 2013, real radar -----------------------------------------------------------------
    fx = json.load(open(os.path.join(TESTS, "data/apophis_2013_radar.json")))
    ref_epoch = fx["reference"]["epoch_jd_tdb"]
    ref = np.array(fx["reference"]["state_au_au_per_day"])
    rows = fx["observations"]
    data = np.empty(len(rows), dtype=[("provID", "U8"), ("obsTime", "U32"), ("stn", "U4"), ("delay", "f8"), ("rmsDelay", "f8"),
                                      ("doppler", "f8"), ("rmsDoppler", "f8"), ("freqTx", "f8")])
    for i, o in enumerate(rows):
        data[i] = ("apophis", o["obsTime"], o["stn"], np.nan, np.nan, np.nan, np.nan, o["freqTx_hz"])
        if o["units"] == "us":
            data[i]["delay"], data[i]["rmsDelay"] = o["value"], o["sigma"]
        else:
            data[i]["doppler"], data[i]["rmsDoppler"] = o["value"], o["sigma"]
    start = ref.copy()
    start[:3] += 1e-4
    start[3:] += 1e-6

    # layup's own inputs, exactly as orbitfit() prepares them
    obsy = LayupObservatory(cache_dir=CACHE)
    et = np.array([spice.str2et(t) for t in data["obsTime"]])
    prep = rfn.append_fields(data, "et", et, usemask=False, asrecarray=True)
    prep = rfn.merge_arrays([prep, obsy.obscodes_to_barycentric(prep)], flatten=True, asrecarray=True, usemask=False)
    prep = _append_observer_acceleration(prep, obsy)
    jd = 2451545.0 + et / 86400.0
    obs_l = [_radar_observation("apophis", r, j, prep.dtype.names) for r, j in zip(prep, jd)]
    obs_radar = obs_l
    common = dict(ra=np.full(len(rows), np.nan), dec=np.full(len(rows), np.nan), kernel=kernel,
                  delay=data["delay"] * 1e-6, sigma_delay=data["rmsDelay"] * 1e-6, doppler=data["doppler"], sigma_doppler=data["rmsDoppler"],
                  frequency=data["freqTx"])
    layup_states = np.column_stack([prep[c] for c in ("x", "y", "z", "vx", "vy", "vz", "ax", "ay", "az")])
    # layup's transmitter: the receiving station extrapolated back by the round trip tau with
    # x - v tau - a tau^2 / 2 (sic; the Taylor expansion has +a tau^2 / 2) and v - a tau. Passed to
    # spacerocks as explicit transmitter states, with tau from the model at the JPL orbit, so the
    # rest of the radar model can be compared like for like.
    r0 = orbfit.residuals(rock(ref, ref_epoch), common["ra"], common["dec"], jd, layup_states, kernel, timescale="tdb", delay=np.zeros(len(rows)))
    tau = (-r0[:, 4] / 86400.0)[:, None]
    p0, v0, a0 = layup_states[:, :3], layup_states[:, 3:6], layup_states[:, 6:]
    layup_tx = np.hstack([p0 - v0 * tau - 0.5 * a0 * tau ** 2, v0 - a0 * tau])
    L, S = report("Apophis 2013 radar, layup's station states and transmitter model", obs_l,
                  dict(common, epoch=jd, timescale="tdb", observer=layup_states, transmitter=layup_tx),
                  [("JPL", ref), ("start", start)], start, ref_epoch)
    times = [Time.from_isot(t) for t in data["obsTime"]]
    L2, S2 = report("Apophis 2013 radar, spacerocks end to end (exact transmitter state)", obs_l,
                    dict(common, epoch=times, observer=list(data["stn"])), [("JPL", ref), ("start", start)], start, ref_epoch)
    for lab, f in (("layup", L), ("spacerocks", S2)):
        s = np.array(f.state)
        print(f"   {lab:10} vs JPL Horizons: position {np.linalg.norm(s[:3] - ref[:3]) * 1.495978707e8:6.2f} km, velocity {np.linalg.norm(s[3:] - ref[3:]) * 1.495978707e11 / 86400:.2e} m/s")
    res = S2.residuals
    print(f"   spacerocks residuals: delay rms {np.sqrt(np.nanmean(res[:, 4] ** 2)) * 1e6:.2f} us, Doppler rms {np.sqrt(np.nanmean(res[:, 5] ** 2)):.3f} Hz")

    # --- 4. streak ADES file through both full pipelines --------------------------------------------
    path = os.path.join(TESTS, "data/streak_observations_complete.csv")
    lay = orbitfit(CSVDataReader(path, "csv", primary_id_column_name="orbitID").read_rows(), CACHE, primary_id_column_name="orbitID")
    r = list(csv.DictReader(open(path)))
    k = np.pi / 180 / 3600 * 24  # arcsec/hour -> rad/day
    f = lambda c: np.array([float(x[c]) for x in r])
    S = orbfit.fit_many([x["orbitID"] for x in r], f("ra") * np.pi / 180, f("dec") * np.pi / 180, [Time.from_isot(x["obsTime"]) for x in r],
                        [x["stn"] for x in r], kernel, ra_rate=f("raRate") * k, dec_rate=f("decRate") * k,
                        sigma_ra_rate=f("rmsRArate") * k, sigma_dec_rate=f("rmsDecrate") * k)
    print("\n== streak ADES file (2 objects), full pipelines (IOD included)")
    for i, oid in enumerate(S["id"]):
        l = lay[lay["orbitID"].astype(str) == oid][0]
        print(f"   {oid}: layup flag {l['flag']} niter {l['niter']} ndof {l['ndof']} chi2 {l['csq']:.6e} epoch {l['epochMJD_TDB'] + 2400000.5:.8f}")
        print(f"   {oid}: spacerocks flag {S['flag'][i]} niter {S['niter'][i]} ndof {S['ndof'][i]} chi2 {S['chi2'][i]:.6e} epoch {S['epoch'][i]:.8f}")
        if l["flag"] == 0 and S["flag"][i] == 0:
            cov = np.array([[l[f"cov_{a}_{b}"] for b in range(6)] for a in range(6)])
            dx = S["state"][i] - np.array([l[c] for c in ("x", "y", "z", "xdot", "ydot", "zdot")])
            print(f"   {oid}: state difference {np.sqrt(dx @ np.linalg.solve(cov, dx)):.2e} sigma")

    # --- 5. full pipelines on the synthetic streak arc, and on Apophis optical + radar ----------
    def pipelines(name, obs_l, jd, kw):
        order = np.argsort(jd, kind="mergesort")
        obs_l = [obs_l[i] for i in order]
        kw = {k: (np.asarray(v)[order] if isinstance(v, (list, np.ndarray)) and np.ndim(v) >= 1 and len(v) == len(jd) else v) for k, v in kw.items()}
        L = LO.do_fit(obs_l, LO._build_sequence(np.asarray(jd)[order], 90.0), CACHE, iod="auto")
        S = orbfit.fit(**kw)
        print(f"\n== {name}: full pipelines (IOD included)")
        print(f"   layup flag {L.flag} niter {L.niter} ndof {L.ndof} chi2 {L.csq:.6e} epoch {L.epoch:.8f}")
        print(f"   spacerocks flag {S.flag} niter {S.niter} ndof {S.ndof} chi2 {S.chi2:.6e} epoch {S.epoch.jd():.8f}")
        if L.flag == 0 and S.flag == 0:
            cov = np.array(L.cov).reshape(6, 6)
            dx = S.state - np.array(L.state)
            print(f"   state difference {np.sqrt(dx @ np.linalg.solve(cov, dx)):.2e} sigma, {np.abs(dx[:3]).max():.2e} AU")
        return L, S

    d = json.load(open(os.path.join(TESTS, "data/streak_synthetic.json")))
    O = d["observations"]
    obs_l = []
    for o in O:
        ob = Observation.from_streak_with_id("s", o["ra"], o["dec"], o["ra_rate"], o["dec_rate"], o["epoch"], o["observer_position"],
                                             o["observer_velocity"], d["rate_unc_radday"], d["rate_unc_radday"])
        ob.ra_unc = ob.dec_unc = d["ra_unc_rad"]
        obs_l.append(ob)
    jd_s = [o["epoch"] for o in O]
    pipelines("synthetic streaks", obs_l, jd_s,
              dict(ra=[o["ra"] for o in O], dec=[o["dec"] for o in O], epoch=jd_s, timescale="tdb",
                   observer=np.array([o["observer_position"] + o["observer_velocity"] for o in O]), kernel=kernel,
                   sigma_ra=d["ra_unc_rad"], sigma_dec=d["ra_unc_rad"], ra_rate=[o["ra_rate"] for o in O], dec_rate=[o["dec_rate"] for o in O],
                   sigma_ra_rate=d["rate_unc_radday"], sigma_dec_rate=d["rate_unc_radday"]))

    # Apophis: synthetic optical astrometry (from the JPL state, spacerocks' ephemeris, 0.3"
    # noise, station G96) over Dec 2012 - Apr 2013, plus the real 2013 radar.
    from spacerocks import RockCollection
    from spacerocks.observing import Observatory
    rng = np.random.default_rng(7)
    t_opt = np.sort(np.concatenate([2456270.5 + n + np.array([0.0, 0.02, 0.04]) for n in (0, 1, 8, 20, 35, 50, 70, 90, 110, 130)]))
    rc = RockCollection()
    rc.add(rock(ref, ref_epoch))
    eph = rc.ephemeris([Time(t, "tdb", "jd") for t in t_opt], Observatory.from_obscode("G96"), kernel)
    sig = 0.3 / 206265
    ra_o = eph["ra"][0] + rng.normal(0, sig, len(t_opt)) / np.cos(eph["dec"][0])
    dec_o = eph["dec"][0] + rng.normal(0, sig, len(t_opt))
    g96 = Observatory.from_obscode("G96")
    st_o = np.array([np.concatenate([o.position, o.velocity, [0, 0, 0]]) for o in (g96.at(Time(t, "tdb", "jd"), kernel) for t in t_opt)])
    obs_o = [Observation.from_astrometry_with_id("apophis", ra_o[i], dec_o[i], t_opt[i], list(st_o[i, :3]), list(st_o[i, 3:6])) for i in range(len(t_opt))]
    for ob in obs_o:
        ob.ra_unc = ob.dec_unc = sig
    n_o, n_r = len(t_opt), len(rows)
    nan_o, nan_r = np.full(n_o, np.nan), np.full(n_r, np.nan)
    jd_all = np.concatenate([t_opt, jd])
    pipelines("Apophis: synthetic optical + real radar", obs_o + obs_radar, jd_all,
              dict(ra=np.concatenate([ra_o, nan_r]), dec=np.concatenate([dec_o, nan_r]), epoch=jd_all, timescale="tdb", kernel=kernel,
                   observer=np.vstack([st_o, layup_states]), sigma_ra=np.concatenate([np.full(n_o, sig), nan_r]),
                   sigma_dec=np.concatenate([np.full(n_o, sig), nan_r]),
                   delay=np.concatenate([nan_o, data["delay"] * 1e-6]), sigma_delay=np.concatenate([nan_o, data["rmsDelay"] * 1e-6]),
                   doppler=np.concatenate([nan_o, data["doppler"]]), sigma_doppler=np.concatenate([nan_o, data["rmsDoppler"]]),
                   frequency=np.concatenate([nan_o, data["freqTx"]]),
                   transmitter=np.vstack([np.full((n_o, 6), np.nan), layup_tx])))

if __name__ == "__main__":
    main()
