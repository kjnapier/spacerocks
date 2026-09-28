"""Shared helpers for the checker validation: kernels, MPC ADES detections."""
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

import spacerocks
from spacerocks import SpiceKernel
from spacerocks.observing import Observatory
from spacerocks import orbfit

ARCSEC = np.pi / (180 * 3600)
CACHE = Path(os.environ.get("MPC_CACHE", "notebooks/examples/mpc_cache"))


def kernel():
    d = Path(os.environ.get("SPACEROCKS_KERNELS", Path.home() / ".spacerocks" / "spice"))
    k = SpiceKernel()
    for f in ["latest_leapseconds.tls", "de440s.bsp", "gm_de440.tpc", "sb441-n16.bsp"]:
        k.load(str(d / f))
    for f in d.glob("earth_*.bpc"):
        k.load(str(f))
    return k


_ground = {}
def is_ground(code):
    if code not in _ground:
        try:
            _ground[code] = Observatory.from_obscode(code).lat is not None
        except Exception:
            _ground[code] = False
    return _ground[code]


def detections(desig, start=None, end=None):
    """Ground-based optical MPC detections of an object: a DataFrame with jd_utc, ra, dec
    (radians), sigma (radians, Veres et al. 2017 unless the observer reported rmsRA/rmsDec) and
    stn."""
    df = pd.DataFrame(json.loads((CACHE / f"obs_{desig}.json").read_text()))
    df.columns = df.columns.str.lower()
    df = df[(df.get("obstype", "optical") == "optical") & df.ra.notna() & df.dec.notna()].copy()
    df = df[df.stn.map(is_ground)]
    t = pd.to_datetime(df.obstime, format="ISO8601", utc=True)
    df["jd_utc"] = 2451545.0 + (t - pd.Timestamp("2000-01-01T12:00:00", tz="UTC")) / pd.Timedelta(days=1)
    if start is not None:
        df = df[t >= pd.Timestamp(start, tz="UTC")]
        t = t[t >= pd.Timestamp(start, tz="UTC")]
    if end is not None:
        df = df[t < pd.Timestamp(end, tz="UTC")]
    df["ra"] = np.radians(pd.to_numeric(df.ra))
    df["dec"] = np.radians(pd.to_numeric(df.dec))
    veres = orbfit.veres_sigma(df.stn.tolist(), df.jd_utc.to_numpy(), catalog=df.astcat.fillna("").tolist(), program=df.prog.fillna("").tolist())
    rms = np.hypot(pd.to_numeric(df.rmsra, errors="coerce"), pd.to_numeric(df.rmsdec, errors="coerce")) / np.sqrt(2) * ARCSEC
    df["sigma"] = np.where(np.isfinite(rms) & (rms > 0), np.maximum(rms, 0.1 * ARCSEC), veres)
    df["desig"] = desig
    return df[["desig", "jd_utc", "ra", "dec", "sigma", "stn"]].reset_index(drop=True)
