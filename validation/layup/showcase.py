"""Five real objects fitted by spacerocks and layup from the same MPC astrometry, for comparison
with JPL (Horizons states, if given) and the MPC's own orbit.

    LAYUP_CACHE=/path SPACEROCKS_KERNELS=/path python showcase.py mpc_cache_dir out.json [jpl.json]

mpc_cache_dir holds the notebook's cache (obs_<id>.json from the MPC's get-obs API, orb_3666.json,
radar_99942.json); jpl.json maps each id to JPL's barycentric J2000 state at JD 2460600.5 TDB.

Both codes get the same detections (2010 on, as the notebook's table; Apophis from its discovery
in 2004, since the Yarkovsky drift needs the whole arc), observer states and
weights (the MPC-reported rmsRA/rmsDec, floored at 0.1"). Each orbit is then refitted with its
epoch set to JD 2460600.5 TDB (2024-10-20), so every state and covariance refers to one instant.
- spacerocks: the robust pipeline, and layup's pipeline (`robust=False`);
- layup: `do_fit(iod="auto")`, then `run_from_vector_with_initial_guess` at the common epoch;
- Apophis: A1 and A2 fitted jointly (JPL's model), optical only for the cross-code comparison, and
  once more with the radar astrometry (spacerocks only; layup's Python API takes radar through
  ADES files).
"""
import json, os, sys, time
import numpy as np
import pandas as pd

import layup.orbitfit as LO
from layup.routines import FitResult, Observation, get_ephem, run_from_vector_with_initial_guess
from spacerocks import SpaceRock, orbfit
from spacerocks.observing import Observatory
from spacerocks.spice import SpiceKernel
from spacerocks.time import Time

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE, OUT = sys.argv[1], sys.argv[2]
JPL = json.load(open(sys.argv[3])) if len(sys.argv) > 3 else {}
EPOCH = 2460600.5
AU_KM = 149597870.7
ARCSEC = np.pi / 180 / 3600
K = os.environ["SPACEROCKS_KERNELS"]
kernel = SpiceKernel()
for f in ["latest_leapseconds.tls", "de440.bsp", "sb441-n16.bsp", "earth_1962_240827_2124_combined.bpc",
          "earth_latest_high_prec.bpc", "pck00010.tpc"]:
    kernel.load(os.path.join(K, f))
ephem = get_ephem(os.environ["LAYUP_CACHE"])
NAMES = {"3666": "Holman", "433": "Eros", "101955": "Bennu", "99942": "Apophis", "90377": "Sedna"}


_ground = {}
def is_ground(code):
    if code not in _ground:
        try:
            _ground[code] = Observatory.from_obscode(code).lat is not None
        except Exception:
            _ground[code] = False
    return _ground[code]


def table(desig, start="2010-01-01"):
    """Optical MPC astrometry since 2010 with observer states, as the notebook's `astrometry`."""
    df = pd.DataFrame(json.load(open(os.path.join(CACHE, f"obs_{desig}.json"))))
    df.columns = df.columns.str.lower()
    num = lambda c: pd.to_numeric(df[c], errors="coerce") if c in df else pd.Series(np.nan, index=df.index)
    t = pd.to_datetime(df["obstime"], format="mixed", utc=True)
    keep = num("ra").notna() & num("dec").notna()
    if start:
        keep &= t >= pd.Timestamp(start, tz="UTC")
    if "deprecated" in df:                            # detections the MPC has withdrawn
        keep &= df["deprecated"].isna()
    # ground stations only (the same observer model in both codes)
    keep &= df["stn"].astype(str).map(is_ground)
    df, t = df[keep], t[keep]
    jd_utc = (t - pd.Timestamp("1970-01-01", tz="UTC")).dt.total_seconds().values / 86400 + 2440587.5
    epochs = [Time(j, "utc", "jd") for j in jd_utc]
    stn = list(df["stn"].astype(str))
    obs = orbfit.observers(stn, epochs, kernel)
    ok = np.isfinite(obs).all(axis=1)
    sig = lambda c: np.maximum(num(c).fillna(1.0).values, 0.1) * ARCSEC
    out = pd.DataFrame({"jd_tdb": [e.tdb().jd() for e in epochs], "ra": np.radians(num("ra").values),
                        "dec": np.radians(num("dec").values), "sigma_ra": sig("rmsra"), "sigma_dec": sig("rmsdec"),
                        "stn": stn})
    out[["ox", "oy", "oz", "ovx", "ovy", "ovz"]] = obs
    return out[ok].sort_values("jd_tdb").reset_index(drop=True)


def sr_fit(d, **kw):
    return orbfit.fit(d["ra"].values, d["dec"].values, d["jd_tdb"].values, d[["ox", "oy", "oz", "ovx", "ovy", "ovz"]].values,
                      kernel, sigma_ra=d["sigma_ra"].values, sigma_dec=d["sigma_dec"].values, timescale="tdb", **kw)


def at_epoch(rock, nongrav=None):
    r = SpaceRock.from_xyz("x", *rock.position, *rock.velocity, rock.epoch, "J2000", "SSB")
    if nongrav is not None:
        r.set_nongrav(*nongrav)
    r.propagate(Time(EPOCH, "tdb", "jd"), kernel)
    return r


def record(state, cov, nongrav=None, nongrav_sigma=None, **extra):
    return dict(state=list(map(float, state)), cov=None if cov is None else np.asarray(cov).tolist(),
                nongrav=None if nongrav is None else list(map(float, nongrav)),
                nongrav_sigma=None if nongrav_sigma is None else [None if not np.isfinite(v) else float(v) for v in nongrav_sigma],
                **extra)


results = {}
for desig, name in NAMES.items():
    d = table(desig, None if desig == "99942" else "2010-01-01")
    ng = "A1A2" if desig == "99942" else None
    mask = 3 if ng else 0
    res = dict(name=name, n=len(d), span=[float(d["jd_tdb"].min()), float(d["jd_tdb"].max())])
    print(f"{name}: {len(d)} detections", flush=True)

    # spacerocks, robust and plain, refitted at the common epoch
    for label, kw in (("spacerocks_robust", dict(robust=True)), ("spacerocks", {})):
        t0 = time.perf_counter()
        f = sr_fit(d, nongrav=ng, **kw)
        dt = time.perf_counter() - t0
        if f.flag not in (0, 2):
            res[label] = dict(flag=int(f.flag), time=dt)
            print(f"  {label}: flag {f.flag}", flush=True)
            continue
        r = at_epoch(f.rock, f.nongrav if ng else None)
        g = sr_fit(d, nongrav=ng, initial=r, **kw)
        res[label] = record(g.state, g.state_covariance, g.nongrav if ng else None, g.nongrav_sigma if ng else None,
                            flag=int(f.flag), time=dt, used=int(g.used.sum()), chi2=float(g.chi2), ndof=int(g.ndof),
                            rms=float(np.sqrt(np.nanmean((g.residuals[g.used, :2] / ARCSEC) ** 2))))
        if label == "spacerocks_robust":
            robust_used, robust_rock = np.array(g.used), g.rock
        print(f"  {label}: flag {f.flag}, {dt:.2f} s, used {g.used.sum()}", flush=True)

    # layup, same inputs
    O = []
    for _, row in d.iterrows():
        o = Observation.from_astrometry_with_id(desig, row["ra"], row["dec"], row["jd_tdb"],
                                                [row["ox"], row["oy"], row["oz"]], [row["ovx"], row["ovy"], row["ovz"]])
        o.ra_unc, o.dec_unc = row["sigma_ra"], row["sigma_dec"]
        O.append(o)
    t0 = time.perf_counter()
    L = LO.do_fit(O, LO._build_sequence(d["jd_tdb"].values, 90.0), os.environ["LAYUP_CACHE"], iod="auto")
    if L.flag == 0 and mask:
        L = run_from_vector_with_initial_guess(ephem, L, O, 100, mask)
    dt = time.perf_counter() - t0
    if L.flag in (0, 2):
        rock = SpaceRock.from_xyz("x", *L.state, Time(L.epoch, "tdb", "jd"), "J2000", "SSB")
        ngL = [L.a1, L.a2, L.a3] if mask else None
        r = at_epoch(rock, ngL)
        g = FitResult()
        g.state, g.epoch, g.flag = list(r.position) + list(r.velocity), EPOCH, 0
        if mask:
            g.a1, g.a2, g.a3 = ngL
        G = run_from_vector_with_initial_guess(ephem, g, O, 100, mask)
        res["layup"] = record(G.state, np.array(G.cov).reshape(-1, int(np.sqrt(len(G.cov))))[:6, :6],
                              [G.a1, G.a2, G.a3] if mask else None, [G.a1_unc, G.a2_unc, G.a3_unc] if mask else None,
                              flag=int(L.flag), refit_flag=int(G.flag), time=dt, chi2=float(G.csq), ndof=int(G.ndof))
    else:
        res["layup"] = dict(flag=int(L.flag), time=dt)
    print(f"  layup: flag {L.flag}, {dt:.2f} s", flush=True)
    if L.flag not in (0, 2) and "state" in res.get("spacerocks_robust", {}):
        # layup's pipeline has no outlier rejection. Give both codes the detections the robust fit
        # kept, starting from its orbit at the common epoch, and fit once more (layup's
        # initial-guess path; spacerocks' fit from `initial`).
        du = d[robust_used].reset_index(drop=True)
        Ou = [o for o, u in zip(O, robust_used) if u]
        g = FitResult()
        g.state, g.epoch, g.flag = list(robust_rock.position) + list(robust_rock.velocity), EPOCH, 0
        if mask:
            g.a1, g.a2, g.a3 = list(res["spacerocks_robust"]["nongrav"])
        G = run_from_vector_with_initial_guess(ephem, g, Ou, 100, mask)
        S = sr_fit(du, nongrav=ng, initial=robust_rock)
        if G.flag in (0, 2):
            res["layup_same_detections"] = record(G.state, np.array(G.cov).reshape(-1, int(np.sqrt(len(G.cov))))[:6, :6],
                                                  [G.a1, G.a2, G.a3] if mask else None, [G.a1_unc, G.a2_unc, G.a3_unc] if mask else None,
                                                  flag=int(G.flag), n=len(Ou))
        res["spacerocks_same_detections"] = record(S.state, S.state_covariance, S.nongrav if ng else None,
                                                   S.nongrav_sigma if ng else None, flag=int(S.flag), n=len(du))
        print(f"  same detections: layup flag {G.flag}, spacerocks flag {S.flag}", flush=True)

    if desig == "99942":
        # spacerocks with the radar astrometry too
        rad = json.load(open(os.path.join(CACHE, "radar_99942.json")))
        rdf = pd.DataFrame(rad["data"], columns=rad["fields"])
        rdf = rdf[rdf["bp"] == "C"]
        MAP = {-1: "251", -2: "254", -9: "256", -13: "252", -14: "253", -25: "257", -38: "255", -73: "259"}
        rc, tx = rdf["rcvr"].astype(int).map(MAP), rdf["xmit"].astype(int).map(MAP)
        rdf, rc, tx = rdf[rc.notna() & tx.notna()], rc[rc.notna() & tx.notna()], tx[rc.notna() & tx.notna()]
        ep = [Time.from_isot(s.replace(" ", "T")) for s in rdf["epoch"]]
        is_delay = (rdf["units"] == "us").values
        val, sg = rdf["value"].astype(float).values, rdf["sigma"].astype(float).values
        n0, n1 = len(d), len(rdf)
        nan = lambda n: np.full(n, np.nan)
        ro = orbfit.observers(list(rc), ep, kernel)
        all_ra = np.concatenate([d["ra"].values, nan(n1)])
        all_dec = np.concatenate([d["dec"].values, nan(n1)])
        all_t = np.concatenate([d["jd_tdb"].values, [e.tdb().jd() for e in ep]])
        all_o = np.vstack([d[["ox", "oy", "oz", "ovx", "ovy", "ovz"]].values, ro])
        kw = dict(sigma_ra=np.concatenate([d["sigma_ra"].values, nan(n1)]), sigma_dec=np.concatenate([d["sigma_dec"].values, nan(n1)]),
                  delay=np.concatenate([nan(n0), np.where(is_delay, val * 1e-6, np.nan)]),
                  sigma_delay=np.concatenate([nan(n0), np.where(is_delay, sg * 1e-6, np.nan)]),
                  doppler=np.concatenate([nan(n0), np.where(is_delay, np.nan, val)]),
                  sigma_doppler=np.concatenate([nan(n0), np.where(is_delay, np.nan, sg)]),
                  frequency=np.concatenate([nan(n0), rdf["freq"].astype(float).values * 1e6]),
                  transmitter=["500"] * n0 + list(tx), timescale="tdb", nongrav="A1A2", robust=True)
        f = orbfit.fit(all_ra, all_dec, all_t, all_o, kernel, **kw)
        r = at_epoch(f.rock, f.nongrav)
        g = orbfit.fit(all_ra, all_dec, all_t, all_o, kernel, initial=r, **kw)
        res["spacerocks_radar"] = record(g.state, g.state_covariance, g.nongrav, g.nongrav_sigma, flag=int(g.flag),
                                         n_radar=int(n1), used=int(g.used.sum()))
        # 2029 close approach from each orbit
        def approach(state, ngp):
            rock = SpaceRock.from_xyz("x", *state, Time(EPOCH, "tdb", "jd"), "J2000", "SSB")
            rock.set_nongrav(*ngp)
            best = (np.inf, None)
            t = Time.from_isot("2029-04-13T20:00:00").tdb()
            for _ in range(240):                                  # every minute for 4 hours
                rock.propagate(t, kernel)
                e = SpaceRock.from_spice("earth", t, "J2000", "SSB", kernel)
                dist = np.linalg.norm(np.subtract(rock.position, e.position)) * AU_KM
                if dist < best[0]:
                    best = (dist, t.utc().iso())
                t = t + 1 / 1440
            return best
        for label in ("spacerocks_robust", "spacerocks_radar", "layup", "layup_same_detections"):
            if "state" in res.get(label, {}):
                res[label]["approach_2029"] = approach(res[label]["state"], res[label]["nongrav"])
        if "99942" in JPL:
            res["jpl_approach_2029"] = approach(JPL["99942"]["state"], JPL["99942"].get("nongrav", [0, 0, 0]))

    if desig == "3666" and os.path.exists(os.path.join(CACHE, "orb_3666.json")):
        orb = json.load(open(os.path.join(CACHE, "orb_3666.json")))
        orb = orb[0] if isinstance(orb, list) else orb
        x = orb["CAR"]["coefficient_values"][:6]
        ep = orb["epoch_data"]
        scale = {"TDT": "tt", "TT": "tt", "TDB": "tdb", "UTC": "utc"}[ep.get("timesystem", "TT").upper()]
        rock = SpaceRock.from_xyz("mpc", *x, Time(ep["epoch"], scale, "mjd"), "ECLIPJ2000", "SUN")
        rock.to_ssb(kernel)
        rock.change_reference_plane("J2000")
        r = at_epoch(rock)
        res["mpc"] = record(list(r.position) + list(r.velocity), None)
    if desig in JPL:
        res["jpl"] = JPL[desig]
    results[desig] = res

json.dump(dict(epoch=EPOCH, objects=results), open(OUT, "w"), indent=1)
print("wrote", OUT)
