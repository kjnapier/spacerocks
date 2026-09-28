"""UTC -> TT in spacerocks.time against ERFA (the C version of SOFA), and TT - UT before 1960
against the Delta T table it follows.

    python compare_erfa.py

1. 1960-2030, one epoch every ~3.7 days plus every table boundary +-1 ms: TT - UTC from
   `Time(jd, "utc", "jd").tt()` against `erfa.utctai` + 32.184 s (TAI - UTC from 1960 to 1972
   runs at an offset rate, then leap seconds). Days that end in a step of TAI - UTC (leap
   seconds, and the 1961-1968 and 1972 steps) are skipped: ERFA reads the UTC Julian date as a
   quasi-JD that spreads the step over that day, while spacerocks applies it at midnight, when
   it happened. The day before 1960 is skipped too (ERFA has no UTC then).
2. Round trips UTC -> TT -> UTC and UTC -> TDB -> UTC, 1600-2030.
3. Before 1960 (no UTC; the epoch is UT), TT - UT against Stephenson, Morrison & Hohenkerk's
   Table S15.2020 as Skyfield bundles it (`pip install skyfield`), if available.
"""
import numpy as np
import erfa
from spacerocks.time import Time

rng = np.random.default_rng(1960)
grid = np.arange(2436934.5, 2462502.5, 3.7)
jd = grid + rng.uniform(0, 1, len(grid))
bounds = [2436934.5, 2437300.5, 2437512.5, 2437665.5, 2438334.5, 2438395.5, 2438486.5, 2438639.5, 2438761.5,
          2438820.5, 2438942.5, 2439004.5, 2439126.5, 2439887.5, 2441317.5]
jd = np.concatenate([jd, np.array(bounds) + 1e-3 / 86400, np.array(bounds[1:]) - 1e-3 / 86400])

# leap-second days (UTC dates whose TAI - UTC changes at the next midnight) are skipped
d1 = erfa.dat(*erfa.jd2cal(jd, 0.0)[:3], 0.0)
d2 = erfa.dat(*erfa.jd2cal(jd + 1, 0.0)[:3], 0.0)
keep = (jd < 2441317.5) | (d1 == d2)
for b in bounds[1:]:
    keep &= ~((jd >= b - 1) & (jd < b))
jd = jd[keep]

# ERFA takes UTC as a quasi-JD (two parts); give it the same instant
ref = np.empty(len(jd))
for i, j in enumerate(jd):
    t1, t2 = erfa.utctai(j, 0.0)
    ref[i] = ((t1 - j) + t2) * 86400 + 32.184
ours = np.array([(Time(j, "utc", "jd").tt().jd() - j) * 86400 for j in jd])
err = np.abs(ours - ref)
print(f"1. TT - UTC, 1960-2030, {len(jd)} epochs: max difference from ERFA {err.max() * 1e6:.1f} us "
      f"(1960-1972: {err[jd < 2441317.5].max() * 1e6:.1f} us); one ulp of a JD is {np.spacing(2.44e6) * 86400e6:.1f} us")

# 2. round trips
jd2 = rng.uniform(2305447.5, 2462502.5, 20000)
rt1 = np.array([Time(j, "utc", "jd").tt().utc().jd() - j for j in jd2]) * 86400e6
rt2 = np.array([Time(j, "utc", "jd").tdb().utc().jd() - j for j in jd2]) * 86400e6
print(f"2. round trips 1600-2030, {len(jd2)} epochs: UTC->TT->UTC max {np.abs(rt1).max():.1f} us, "
      f"UTC->TDB->UTC max {np.abs(rt2).max():.1f} us")

# 3. before 1960
try:
    from skyfield.api import load
    import skyfield.timelib as tl
    s15 = tl.load_bundled_npy("delta_t.npz")["Table-S15.2020.txt"]
    sp = tl.Splines(s15)
    jd3 = rng.uniform(1458085.5, 2436934.5, 20000)                  # -720 to 1960
    year = 2000.0 + (jd3 - 2451545.0) / 365.25
    ours = np.array([(Time(j, "utc", "jd").tt().jd() - j) * 86400 for j in jd3])
    print(f"3. TT - UT before 1960, {len(jd3)} epochs: max difference from Table S15.2020 "
          f"{np.abs(ours - sp(year)).max() * 1e6:.1f} us")
except ImportError:
    print("3. skipped (no skyfield)")

for y, j in [(1850, 2396758.5), (1900, 2415020.5), (1950, 2433282.5), (1960, 2436934.5), (1965, 2438761.5),
             (1972, 2441317.5), (2000, 2451544.5)]:
    print(f"   {y}: TT - UTC = {(Time(j, 'utc', 'jd').tt().jd() - j) * 86400:9.4f} s")
