"""Compare the Rust SPICE implementation with CSPICE (spiceypy) on synthetic and real kernels."""
import os, subprocess, sys, numpy as np, spiceypy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
EXE = os.environ.get("SPICE_CHECK_EXE", os.path.join(HERE, "..", "..", "target", "release", "examples", "spice_check"))
rng = np.random.default_rng(7)


def run(cmds):
    p = subprocess.run([EXE], input="\n".join(cmds) + "\n", capture_output=True, text=True)
    if p.returncode != 0:
        print(p.stderr); raise SystemExit(1)
    out = []
    for line in p.stdout.splitlines():
        if line.startswith("ok"):
            out.append(np.array([float(x) for x in line.split()[1:]]))
        else:
            out.append(line)
    return out


def rel_err(a, b):
    a = np.asarray(a); b = np.asarray(b)
    return np.linalg.norm(a - b) / max(np.linalg.norm(b), 1e-300)


failures = 0
report = []


def compare_states(label, loads, pairs, epochs, tol=1e-12, frame=None):
    """pairs: list of (target, observer)."""
    global failures
    sp.kclear()
    for f in loads:
        sp.furnsh(f)
    cmds = [f"load {f}" for f in loads]
    refs = []
    for t, o in pairs:
        for et in epochs:
            try:
                if frame:
                    st, _ = sp.spkezr(str(t), et, frame, "NONE", str(o))
                else:
                    st, _ = sp.spkgeo(t, et, "J2000", o)
                refs.append(np.array(st))
            except Exception as e:
                refs.append(None)
                sp.reset()
            cmds.append(f"state {t} {o} {frame} {float(et)!r}" if frame else f"spkgeo {t} {o} {float(et)!r}")
    res = run(cmds)[len(loads):]
    worst_p = worst_v = 0.0
    nerr = 0; nboth_err = 0
    for r, ref in zip(res, refs):
        if ref is None:
            if isinstance(r, str):
                nboth_err += 1
                continue
            nerr += 1
            continue
        if isinstance(r, str):
            nerr += 1
            if nerr <= 3:
                print("   rust error:", r[:200])
            continue
        worst_p = max(worst_p, rel_err(r[:3], ref[:3]))
        worst_v = max(worst_v, rel_err(r[3:], ref[3:]))
    ok = worst_p < tol and worst_v < tol and nerr == 0
    if not ok:
        failures += 1
    line = f"{'PASS' if ok else 'FAIL'}  {label:<42} n={len(refs):5d}  max rel err pos={worst_p:.2e} vel={worst_v:.2e}  mismatched errors={nerr} (both error: {nboth_err})"
    print(line); report.append(line)


def compare_xform(label, loads, pairs, epochs, tol=1e-12):
    global failures
    sp.kclear()
    for f in loads:
        sp.furnsh(f)
    cmds = [f"load {f}" for f in loads]
    refs = []
    for a, b in pairs:
        for et in epochs:
            try:
                m = np.array(sp.sxform(a, b, et))
                refs.append(np.concatenate([m[:3, :3].flatten(), m[3:, :3].flatten()]))
            except Exception as e:
                refs.append(None); sp.reset()
            cmds.append(f"sxform {a} {b} {float(et)!r}")
    res = run(cmds)[len(loads):]
    worst_r = worst_d = 0.0; nerr = 0
    for r, ref in zip(res, refs):
        if ref is None or isinstance(r, str):
            if not (ref is None and isinstance(r, str)):
                nerr += 1
                if nerr <= 3:
                    print("   mismatch:", r if isinstance(r, str) else "cspice error")
            continue
        worst_r = max(worst_r, np.abs(r[:9] - ref[:9]).max())
        worst_d = max(worst_d, rel_err(r[9:], ref[9:]))
    ok = worst_r < tol and worst_d < 1e-9 and nerr == 0
    if not ok:
        failures += 1
    line = f"{'PASS' if ok else 'FAIL'}  {label:<42} n={len(refs):5d}  max abs err R={worst_r:.2e}  rel err dR={worst_d:.2e}  mismatched errors={nerr}"
    print(line); report.append(line)


K = os.environ.get("SPICE_SYNTH_DIR", "/tmp/spacerocks_synthetic_kernels")
T0, T1 = -1.0e8, 1.0e8
ep = np.concatenate([rng.uniform(T0, T1, 300), [T0, T1, 0.0], T0 + (T1 - T0) / 50 * np.arange(51)])
ep = np.clip(ep, T0, T1)

synth = [f"{K}/synthetic.bsp"]
for tgt, lab in [(1002, "type 2"), (1003, "type 3"), (1005, "type 5"), (1008, "type 8 (odd window)"),
                 (1018, "type 8 (even window)"), (1012, "type 12"), (1122, "type 12 deg5"), (1009, "type 9 (odd)"),
                 (1019, "type 9 (even)"), (1013, "type 13"), (1023, "type 13 deg3"), (1014, "type 14"),
                 (1180, "type 18 subtype 0"), (1181, "type 18 subtype 1"), (1020, "type 20"),
                 (1001, "type 1"), (1021, "type 21"), (1190, "type 19 (select last)"), (1191, "type 19 (select first)")]:
    compare_states(f"SPK {lab}", synth, [(tgt, 10)], ep)
compare_states("SPK type 15", synth, [(1015, 399)], ep, tol=1e-11)
compare_states("SPK type 17", synth, [(1017, 399)], ep, tol=1e-11)
compare_states("SPK frames + segment priority (1022)", synth, [(1022, 10)], ep)
compare_states("SPK chaining (1002 wrt 1003, 1015 wrt 1002)", synth, [(1002, 1003), (1015, 1002), (10, 1013)], ep[:120])

# Boundary epochs for the type 19 intervals and type 9/13 interior epochs
bnd = np.array([-1e7, 4e7, -1e7 + 1e-6, 4e7 - 1e-6])
compare_states("SPK type 19 interval boundaries", synth, [(1190, 10), (1191, 10)], bnd)

# Frames and orientation
meta = [f"{K}/meta.tm"]
compare_xform("PCK type 2 (ITRF93 <- ECLIPJ2000)", meta, [("ITRF93", "J2000"), ("J2000", "ITRF93")], ep[:200])
compare_xform("PCK type 20 (kernel-defined frame)", meta, [("MYPCK20", "J2000")], ep[:200])
compare_xform("PCK type 3 (generic segment)", meta, [("MYPCK3", "ECLIPJ2000")], ep[:200])
compare_xform("TK frames (matrix/angles/quaternion)", meta,
              [("MYTK_MAT", "J2000"), ("MYTK_ANG", "J2000"), ("MYTK_Q", "ITRF93"), ("MYTK_ANG", "GALACTIC")], ep[:100])
compare_xform("text PCK IAU model (IAU_MARS)", meta, [("IAU_MARS", "J2000"), ("IAU_MARS", "ECLIPJ2000")],
              rng.uniform(-3e9, 3e9, 200), tol=5e-12)
compare_xform("built-in inertial frames", meta,
              [("B1950", "J2000"), ("FK4", "ECLIPJ2000"), ("GALACTIC", "J2000"), ("MARSIAU", "DE-143"),
               ("ECLIPB1950", "DE-118")], [0.0])
compare_states("states in non-inertial frame (ITRF93)", meta, [(1002, 10), (1013, 1002)], ep[:100], frame="ITRF93")
compare_states("states in TK frame (MYTK_ANG)", meta, [(1002, 10)], ep[:100], frame="MYTK_ANG")

# Meta-kernel name definitions
r = run([f"load {K}/meta.tm", "bodyid MY ROCK", "bodyid other_rock", "bodyid EARTH"])
print(("PASS" if [x[0] for x in r[1:]] == [1002, 1003, 399] else "FAIL"), " kernel-defined body names", [list(x) if not isinstance(x, str) else x for x in r[1:]])

# ---------------------------------------------------------------------------------------------
# Real kernels
# ---------------------------------------------------------------------------------------------
R = os.environ.get("SPACEROCKS_KERNELS")
if not R:
    print("SPACEROCKS_KERNELS not set; skipping real-kernel comparisons"); print("FAILURES:", failures); raise SystemExit(failures)
de = f"{R}/de440s.bsp"
sp.kclear(); sp.furnsh(de)
cov = sp.spkcov(de, 399); lo, hi = cov[0], cov[-1]
ep_de = np.concatenate([rng.uniform(lo, hi, 400), [lo, hi]])
bodies = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 199, 299, 301, 399]
pairs = [(t, 0) for t in bodies if t != 0] + [(399, 10), (301, 399), (10, 301), (5, 399), (199, 301), (0, 399)]
compare_states("de440s (type 2) many pairs", [de], pairs, ep_de[:150])

jw = f"{R}/jwst_pred.bsp"
sp.kclear(); sp.furnsh(jw)
cj = sp.spkcov(jw, -170)
jlo, jhi = cj[0], cj[-1]
# include segment boundaries
hs = sp.dafopr(jw); sp.dafbfs(hs); seg_b = []
while sp.daffna():
    dc, ic = sp.dafus(sp.dafgs(), 2, 6); seg_b += [dc[0], dc[1]]
sp.dafcls(hs)
seg_b = np.array(seg_b)
ep_j = np.concatenate([rng.uniform(jlo, jhi, 500), seg_b[rng.choice(len(seg_b), 200)]])
compare_states("jwst_pred (445 type 13 segments)", [jw], [(-170, 399)], ep_j)
compare_states("jwst_pred + de440s (JWST wrt SSB, Sun)", [de, jw], [(-170, 0), (-170, 10)], ep_j[:200])

pm = f"{R}/polymele_20220626a.bsp"
sp.kclear(); sp.furnsh(pm)
cp = sp.spkcov(pm, 2015094)
ep_p = rng.uniform(cp[0], cp[-1], 300)
compare_states("polymele (type 3)", [pm], [(2015094, 0)], ep_p)

for bpc in ["earth_1962_240827_2124_combined.bpc", "earth_latest_high_prec.bpc"]:
    f = f"{R}/{bpc}"
    sp.kclear(); sp.furnsh(f)
    c = sp.pckcov(f, 3000)
    ep_b = np.concatenate([rng.uniform(c[0], c[-1], 300), [c[0], c[-1]]])
    compare_xform(f"{bpc[:28]} ITRF93", [f], [("ITRF93", "J2000"), ("ITRF93", "ECLIPJ2000")], ep_b)

print()
print("FAILURES:", failures)
open(os.path.join(K, "validation_report.txt"), "w").write("\n".join(report) + f"\nFAILURES: {failures}\n")
