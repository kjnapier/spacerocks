"""Synthetic short arcs (TNOs, Centaurs, main-belt objects; 1-4 nights over up to 14 days) in
layup's CSV format, to exercise the Bernstein-Khushalani fallback. Astrometry from spacerocks'
IAS15 ephemerides (ASSIST force model, light time), Gaussian noise of 0.1", station X05.

    SPACEROCKS_KERNELS=/path python make_short_arcs.py out.csv [seed]
"""
import os, sys
import numpy as np
from spacerocks import SpaceRock, RockCollection
from spacerocks.spice import SpiceKernel
from spacerocks.time import Time
from spacerocks.observing import Observatory

K = os.environ["SPACEROCKS_KERNELS"]
kernel = SpiceKernel()
for f in ["latest_leapseconds.tls", "de440.bsp", "sb441-n16.bsp", "earth_1962_240827_2124_combined.bpc"]:
    kernel.load(os.path.join(K, f))
rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 1)

patterns = {"2n": [0, 1], "3n": [0, 1, 3], "4n": [0, 2, 5, 9], "1n": [0], "wk": [0, 7, 14]}
pops = {"tno": (38, 46, 0.1, 25), "cen": (8, 20, 0.3, 25), "mba": (1.9, 3.0, 0.2, 20)}
site = Observatory.from_obscode("X05")
start = Time.from_isot("2025-09-10T04:00:00.000").utc().jd()
rows = []
for pop, (qlo, qhi, emax, imax) in pops.items():
    for pat, nights in patterns.items():
        for k in range(6):
            q, e = rng.uniform(qlo, qhi), rng.uniform(0, emax)
            rock = SpaceRock.from_kepler(f"{pop}-{pat}-{k}", q, e, np.radians(rng.uniform(0, imax)), rng.uniform(0, 2 * np.pi),
                                         rng.uniform(0, 2 * np.pi), rng.uniform(0, 2 * np.pi), Time(start, "utc", "jd"), "ECLIPJ2000", "SUN")
            per_night = 4 if pat == "1n" else 3
            step = 0.05 if pat == "1n" else 0.03
            jds = [start + n + i * step for n in nights for i in range(per_night)]
            # round to whole milliseconds so the ISO strings are exact
            isos = [Time(jd, "utc", "jd").iso() for jd in jds]
            times = [Time.from_isot(s.rstrip("Z")) for s in isos]
            rc = RockCollection()
            rc.add(rock)
            eph = rc.ephemeris(times, site, kernel)
            ra = eph["ra"][0] + rng.normal(0, 0.1 / 206265, len(jds)) / np.cos(eph["dec"][0])
            dec = eph["dec"][0] + rng.normal(0, 0.1 / 206265, len(jds))
            for s, a, d in zip(isos, ra, dec):
                rows.append(f"{pop}-{pat}-{k},{float(np.degrees(a) % 360)!r},{float(np.degrees(d))!r},{s},X05")
with open(sys.argv[1], "w") as f:
    f.write("provID,ra,dec,obsTime,stn\n" + "\n".join(rows) + "\n")
print(len(rows), "rows")
