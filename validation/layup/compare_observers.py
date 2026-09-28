"""Observer states from ADES positions: spacerocks' `orbfit.observers` vs layup's
`LayupObservatory.obscodes_to_barycentric`, and ground stations before 1962 (IAU_EARTH).

Cases: satellites with ICRF_KM positions (with and without velocities), ICRF_AU positions and
velocities, roving observers (247) with WGS84 positions, ground stations, and ground stations
between 1901 and 1962. Both codes get the same TDB instants (their pre-1972 UTC conventions
differ; see the README).

Before 1962 spacerocks corrects IAU_EARTH for the Earth's actual rotation (UT1, via Delta T),
which layup does not, so there the stations are also compared with layup's positions turned to
IAU_EARTH evaluated `iau_earth_rotation_delay` earlier (needs skyfield for Delta T).

    LAYUP_CACHE=/path/to/layup/cache SPACEROCKS_KERNELS=/path/to/kernels python compare_observers.py

The layup cache's meta-kernel must include pck00010 (layup's own kernel list does).
"""
import os
import numpy as np
import numpy.lib.recfunctions as rfn
import spiceypy as spice

from layup.utilities.data_processing_utilities import LayupObservatory
from spacerocks import orbfit
from spacerocks.spice import SpiceKernel
from spacerocks.time import Time

K = os.environ["SPACEROCKS_KERNELS"]
kernel = SpiceKernel()
for f in ["latest_leapseconds.tls", "de440.bsp", "earth_1962_240827_2124_combined.bpc", "earth_latest_high_prec.bpc", "pck00010.tpc"]:
    kernel.load(os.path.join(K, f))
lo = LayupObservatory(cache_dir=os.environ["LAYUP_CACHE"])
rng = np.random.default_rng(147)
AU_KM = 149597870.7
n = 60

rows = []
def add(stn, t, sys="", pos=(np.nan,) * 3, vel=(np.nan,) * 3):
    rows.append(("x", t, stn, sys, 399, *pos, *vel))

times = [f"{y}-{m:02d}-{d:02d}T{h:02d}:{mi:02d}:{s:02d}" for y, m, d, h, mi, s in
         zip(rng.integers(1990, 2026, n), rng.integers(1, 13, n), rng.integers(1, 28, n), rng.integers(0, 24, n),
             rng.integers(0, 60, n), rng.integers(0, 60, n))]
for t in times:
    u = rng.normal(size=3); u /= np.linalg.norm(u)
    r_km = u * rng.uniform(6800, 42000)                                  # LEO to GEO
    v_kms = np.cross(u, rng.normal(size=3)); v_kms *= rng.uniform(3, 7.5) / np.linalg.norm(v_kms)
    add("C51", t, "ICRF_KM", r_km)                                        # position only
    add("C53", t, "ICRF_KM", r_km * 1.01, v_kms)                          # with velocity (layup caches per station and epoch)
    add("C57", t, "ICRF_AU", r_km / AU_KM * 3, v_kms * 86400 / AU_KM)     # AU, AU/day
    add("247", t, "WGS84", (rng.uniform(-180, 180), rng.uniform(-70, 70), rng.uniform(0, 3000)))
    add("G96", t)
old = [f"{y}-{m:02d}-{d:02d}T{h:02d}:13:07" for y, m, d, h in
       zip(rng.integers(1901, 1962, n), rng.integers(1, 13, n), rng.integers(1, 28, n), rng.integers(0, 24, n))]
for t in old:
    add(["024", "675", "511", "807"][rng.integers(4)], t)

dt = [("provID", "U4"), ("obsTime", "U32"), ("stn", "U4"), ("sys", "U8"), ("ctr", "i8"),
      ("pos1", "f8"), ("pos2", "f8"), ("pos3", "f8"), ("vel1", "f8"), ("vel2", "f8"), ("vel3", "f8")]
data = np.array(rows, dtype=dt)
et = np.array([spice.str2et(t) for t in data["obsTime"]])
data = rfn.append_fields(data, "et", et, usemask=False, asrecarray=True)
L = np.array([[r[i] for i in range(6)] for r in lo.obscodes_to_barycentric(data)])

ep = [Time(2451545.0 + e / 86400.0, "tdb", "jd") for e in et]
pos = np.column_stack([data["pos1"], data["pos2"], data["pos3"]])
vel = np.column_stack([data["vel1"], data["vel2"], data["vel3"]])
sys = [s or None for s in data["sys"]]
S = orbfit.observers(list(data["stn"]), ep, kernel, sys=sys, pos=pos, vel=vel)

dp = np.linalg.norm(S[:, :3] - L[:, :3], axis=1) * AU_KM * 1e3
dv = np.linalg.norm(S[:, 3:] - L[:, 3:], axis=1) * AU_KM / 86400 * 1e3
labels = np.array([f"{s} {y or 'station'}" for s, y in zip(data["stn"], data["sys"])], dtype=object)
labels[et < spice.str2et("1962-Jan-20")] = "stations before 1962"

# Before 1962: layup's station, rotated from IAU_EARTH at et to IAU_EARTH at et - delay.
try:
    from skyfield import timelib as tl
    s15 = tl.Splines(tl.load_bundled_npy("delta_t.npz")["Table-S15.2020.txt"])
    old_rows = np.where(et < spice.str2et("1962-Jan-20"))[0]
    dq = []
    for i in old_rows:
        year = 2000.0 + et[i] / (365.25 * 86400.0)
        delay = s15(year) - 75.025 - 0.5620 * (year - 2000.0)
        e = np.array(spice.spkez(399, et[i], "J2000", "NONE", 0)[0][:3]) / AU_KM
        fixed = spice.pxform("J2000", "IAU_EARTH", et[i]) @ (L[i, :3] - e)
        pred = e + spice.pxform("IAU_EARTH", "J2000", et[i] - delay) @ fixed
        dq.append(np.linalg.norm(S[i, :3] - pred) * AU_KM * 1e3)
    print(f"before 1962, against layup with the UT1 correction: position max {max(dq):.3f} m")
except ImportError:
    pass

for lab in dict.fromkeys(labels):
    m = labels == lab
    print(f"{lab:22} {m.sum():4d} rows: position max {dp[m].max():8.3f} m, velocity max {dv[m].max():.2e} m/s")
