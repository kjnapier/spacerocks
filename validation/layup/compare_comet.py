"""Original orbits of long-period comets: `orbfit.comet_orbits` against the CODE catalogue and
layup's `comet`, on layup's test sample (tests/data/code_LPCs.csv: osculating elements of 369
comets, and code_LPCs_originals.csv: the catalogue's original 1/a with its uncertainty).

    LAYUP_TESTS=/path/to/layup/tests LAYUP_CACHE=/path SPACEROCKS_KERNELS=/path python compare_comet.py [n]

(layup's comet converts elements without passing its cache directory; the script passes it.)

The osculating elements are turned into barycentric states by layup's `convert` (the same
states for both codes). Reported: the difference from the catalogue in units of its quoted
uncertainty (as layup's test_apply_comet), for spacerocks and layup, and between the two codes.
The two comets layup's test excludes for their non-gravitational acceleration are listed apart.
"""
import os, sys, time
import numpy as np
from layup.utilities.file_io.CSVReader import CSVDataReader
from layup.convert import convert
import functools
import layup.comet as layup_comet
from layup.comet import _apply_comet
from layup.utilities.layup_configs import LayupConfigs
from spacerocks import SpaceRock, orbfit
from spacerocks.spice import SpiceKernel
from spacerocks.time import Time

T = os.path.join(os.environ["LAYUP_TESTS"], "data")
K = os.environ["SPACEROCKS_KERNELS"]
kernel = SpiceKernel()
for f in ["latest_leapseconds.tls", "de440.bsp", "sb441-n16.bsp"]:
    kernel.load(os.path.join(K, f))
n = int(sys.argv[1]) if len(sys.argv) > 1 else None
KNOWN_NONGRAV = {"C/2020 S4", "C/2007 Q3"}

data = CSVDataReader(os.path.join(T, "code_LPCs.csv"), "csv", primary_id_column_name="ObjID").read_rows()
if n:
    data = data[:n]
ref = CSVDataReader(os.path.join(T, "code_LPCs_originals.csv"), "csv", primary_id_column_name="ObjID").read_rows()
ref = {str(r["ObjID"]): r for r in ref}
layup_comet.convert = functools.partial(convert, cache_dir=os.environ["LAYUP_CACHE"])
bc = convert(data, convert_to="BCART_EQ", cache_dir=os.environ["LAYUP_CACHE"], primary_id_column_name="ObjID")

t0 = time.time()
ours = {}
for r in bc:
    rock = SpaceRock.from_xyz(str(r["ObjID"]), r["x"], r["y"], r["z"], r["xdot"], r["ydot"], r["zdot"],
                              Time(r["epochMJD_TDB"] + 2400000.5, "tdb", "jd"), "J2000", "SSB")
    ours[str(r["ObjID"])] = orbfit.comet_orbits(rock, kernel)["original"]
t_ours = time.time() - t0


class Args:
    primary_id_column_name = "ObjID"; n = 1; chunk = 10000; ar_data_file_path = os.environ["LAYUP_CACHE"]; force = True; code_format = True


t0 = time.time()
lay = _apply_comet(data, Args(), LayupConfigs().auxiliary, primary_id_column_name="ObjID")
t_lay = time.time() - t0
lay = {str(r["ObjID"]): r for r in lay}

rows = []
for oid, o in ours.items():
    if oid not in ref or oid not in lay or o is None:
        continue
    c, s = float(ref[oid]["inv_ao"]), float(ref[oid]["dinv_ao"])
    if s <= 0:
        continue
    mine, theirs = o["inv_a"] * 1e6, float(lay[oid]["inv_ao_CODE"])
    rows.append((oid, abs(mine - c) / s, abs(theirs - c) / s, abs(mine - theirs), o["reached"], o["distance"]))
z_s = np.array([r[1] for r in rows if r[0] not in KNOWN_NONGRAV])
z_l = np.array([r[2] for r in rows if r[0] not in KNOWN_NONGRAV])
print(f"{len(rows)} comets with a catalogue uncertainty ({sum(r[4] for r in rows)} integrated all the way to 250 AU by spacerocks; "
      f"the rest stopped where the ephemeris ends, at {min(r[5] for r in rows):.0f} AU or more)")
print(f"|spacerocks - CODE| / sigma: median {np.median(z_s):.4f}, 90th percentile {np.percentile(z_s, 90):.3f}, max {z_s.max():.2f}")
print(f"|layup - CODE| / sigma:      median {np.median(z_l):.4f}, 90th percentile {np.percentile(z_l, 90):.3f}, max {z_l.max():.2f}")
d = np.array([r[3] for r in rows])
print(f"|spacerocks - layup| (1e-6/AU): median {np.median(d):.4f}, max {d.max():.3f}")
for r in rows:
    if r[0] in KNOWN_NONGRAV:
        print(f"  {r[0]} (non-gravitational, excluded above): {r[1]:.1f} sigma spacerocks, {r[2]:.1f} sigma layup")
print(f"time: spacerocks {t_ours:.1f} s, layup {t_lay:.1f} s")
