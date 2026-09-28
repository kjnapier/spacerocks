"""Generate synthetic kernels of every supported segment type with CSPICE's own writers."""
import ctypes, os, numpy as np, spiceypy as sp
from spiceypy.utils.libspicehelper import libspice as L

OUT = os.environ.get("SPICE_SYNTH_DIR", "/tmp/spacerocks_synthetic_kernels")
os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(12345)

c_int = ctypes.c_int  # f2c integer (32-bit)
c_dbl = ctypes.c_double
def I(v): return ctypes.byref(c_int(v))
def D(v): return ctypes.byref(c_dbl(v))
def DA(arr):
    a = np.ascontiguousarray(arr, dtype=np.float64)
    return a.ctypes.data_as(ctypes.POINTER(c_dbl)), a
def IA(arr):
    a = np.ascontiguousarray(arr, dtype=np.int32)
    return a.ctypes.data_as(ctypes.POINTER(c_int)), a
def S(s): return ctypes.c_char_p(s.encode())

def fresh(path):
    if os.path.exists(path): os.remove(path)
    return path

def check_err():
    if sp.failed():
        msg = sp.getmsg("LONG", 400); sp.reset(); raise RuntimeError(msg)

sp.erract("SET", 10, "RETURN")

# --------------------------------------------------------------------------------------------
# SPK file with many segment types. Targets 1001.. relative to center 10 (Sun) unless noted.
# --------------------------------------------------------------------------------------------
spk = fresh(f"{OUT}/synthetic.bsp")
h = sp.spkopn(spk, "synthetic", 4000)
T0, T1 = -1.0e8, 1.0e8
GM = 1.32712440041e11

def kepler_states(epochs, a=2.5e8, e=0.3, inc=0.4):
    out = []
    for t in epochs:
        el = [a*(1-e), e, inc, 0.7, 1.1, 0.3, 0.0, GM]
        out.append(sp.conics(el, t))
    return np.array(out)

# Type 2 and 3: random Chebyshev coefficients
for typ, body in [(2, 1002), (3, 1003)]:
    n = 50; intlen = (T1 - T0) / n; deg = 11
    ncomp = 3 if typ == 2 else 6
    coeffs = rng.normal(size=(n, ncomp, deg + 1)) * np.array([1e8 / (k + 1) ** 2 for k in range(deg + 1)])
    if typ == 2:
        sp.spkw02(h, body, 10, "J2000", T0, T1, f"type{typ}", intlen, n, deg, coeffs.flatten(), T0)
    else:
        sp.spkw03(h, body, 10, "J2000", T0, T1, f"type{typ}", intlen, n, deg, coeffs.flatten(), T0)
# Type 2 in a non-J2000 inertial frame (tests frame rotation), and an overlapping
# higher-priority segment for the same body covering part of the interval.
coeffs = rng.normal(size=(10, 3, 8)) * 1e7
sp.spkw02(h, 1022, 10, "ECLIPJ2000", T0, T1, "type2 eclip", (T1 - T0) / 10, 10, 7, coeffs.flatten(), T0)
coeffs = rng.normal(size=(10, 3, 8)) * 1e7
sp.spkw02(h, 1022, 10, "GALACTIC", -2e7, 3e7, "type2 galactic override", 5e6, 10, 7, coeffs.flatten(), -2e7)

# Type 5: two-body between discrete states
ep5 = np.sort(rng.uniform(T0, T1, 237)); ep5[0] = T0; ep5[-1] = T1
sp.spkw05(h, 1005, 10, "J2000", T0, T1, "type5", GM, len(ep5), kepler_states(ep5, e=0.7), ep5)

# Types 8 and 12 (equal spacing), odd and even window sizes
for typ, body, deg in [(8, 1008, 7), (8, 1018, 6), (12, 1012, 7), (12, 1022 + 100, 5)]:
    n = 311; step = (T1 - T0) / (n - 1)
    ep = T0 + step * np.arange(n)
    st = kepler_states(ep) + rng.normal(size=(n, 6)) * [1e3, 1e3, 1e3, 1e-3, 1e-3, 1e-3]
    if typ == 8:
        sp.spkw08(h, body, 10, "J2000", T0, T1, f"type8 deg{deg}", deg, n, st, T0, step)
    else:
        sp.spkw12(h, body, 10, "J2000", T0, T1, f"type12 deg{deg}", deg, n, st, T0, step)

# Types 9 and 13 (unequal spacing), including >100 epochs so directories are exercised
for typ, body, deg in [(9, 1009, 7), (9, 1019, 8), (13, 1013, 7), (13, 1023, 3)]:
    n = 263
    ep = np.sort(rng.uniform(T0, T1, n)); ep[0] = T0; ep[-1] = T1
    st = kepler_states(ep) + rng.normal(size=(n, 6)) * [1e3, 1e3, 1e3, 1e-3, 1e-3, 1e-3]
    if typ == 9:
        sp.spkw09(h, body, 10, "J2000", T0, T1, f"type9 deg{deg}", deg, n, st, ep)
    else:
        sp.spkw13(h, body, 10, "J2000", T0, T1, f"type13 deg{deg}", deg, n, st, ep)

# Type 14: Chebyshev with unequal steps (generic segment)
sp.spk14b(h, "type14", 1014, 10, "J2000", T0, T1, 9)
bounds = np.sort(rng.uniform(T0, T1, 120)); bounds[0] = T0; bounds[-1] = T1
recs, starts = [], []
for i in range(len(bounds) - 1):
    mid = 0.5 * (bounds[i] + bounds[i + 1]); rad = 0.5 * (bounds[i + 1] - bounds[i])
    rec = [mid, rad] + list(rng.normal(size=6 * 10) * 1e6)
    recs.extend(rec); starts.append(bounds[i])
sp.spk14a(h, len(starts), recs, starts)
sp.spk14e(h)

# Type 15: precessing conic
sp.spkw15(h, 1015, 399, "J2000", T0, T1, "type15", 1000.0, [0.0, -0.99, 0.1], [1.0, 0.0, 0.0], 9000.0, 0.02, 0.0,
          [0.0, 0.0, 1.0], 398600.4418, 1.0826e-3, 6378.137)

# Type 17: equinoctial elements
sp.spkw17(h, 1017, 399, "J2000", T0, T1, "type17", 0.0,
          [42164.0, 0.01, 0.02, 0.5, 0.03, 0.04, 1e-7, 7.29e-5, -2e-8], 0.1, 1.2)

# Type 18: subtypes 0 (Hermite, 12-element packets) and 1 (Lagrange, 6-element packets)
for sub, body, deg18 in [(0, 1180, 7), (1, 1181, 9)]:
    n = 211
    ep = np.sort(rng.uniform(T0, T1, n)); ep[0] = T0; ep[-1] = T1
    st = kepler_states(ep)
    if sub == 0:
        pk = np.hstack([st, st[:, 3:] + rng.normal(size=(n, 3)) * 1e-3, rng.normal(size=(n, 3)) * 1e-9])
    else:
        pk = st
    sp.spkw18(h, sub, body, 10, "J2000", T0, T1, f"type18 sub{sub}", deg18, pk, ep)

# Type 20: Chebyshev velocity only
n = 40; intlen_days = (T1 - T0) / 86400.0 / n; deg = 9
cdata = []
for r in range(n):
    for c in range(3):
        cdata.extend(list(rng.normal(size=deg + 1) * 1e-3))  # velocity coeffs (in DSCALE/TSCALE units)
        cdata.append(rng.normal() * 1.0)  # position at midpoint (DSCALE units)
initjd = 2451545.0 - 1158.0; initfr = T0 / 86400.0 + 1158.0
sp.spkw20(h, 1020, 10, "J2000", T0, T1, "type20", intlen_days, n, deg, cdata, 1.495978707e8, 86400.0, initjd, initfr)

# Type 1 and 21 (modified difference arrays) via the Fortran interface
def mda_records(n, maxdim, epochs):
    size = 4 * maxdim + 11
    recs = np.zeros((n, size))
    for i in range(n):
        r = recs[i]
        r[0] = epochs[i] - 0.5 * (epochs[1] - epochs[0] if n > 1 else 1e5)  # TL
        r[1:1 + maxdim] = rng.uniform(1e4, 1e5, maxdim) * np.arange(1, maxdim + 1)  # G
        r[maxdim + 1:maxdim + 7] = rng.normal(size=6) * [1e8, 10, 1e8, 10, 1e8, 10]
        for c in range(3):
            r[(c + 1) * maxdim + 7:(c + 2) * maxdim + 7] = rng.normal(size=maxdim) * 1e-6 / np.arange(1, maxdim + 1) ** 2
        kq = maxdim - 2
        r[4 * maxdim + 7] = kq + 1
        r[4 * maxdim + 8:4 * maxdim + 11] = [kq, kq - 1, kq - 2]
    return recs

for typ, body, maxdim in [(1, 1001, 15), (21, 1021, 20)]:
    n = 230
    ep = np.linspace(T0, T1, n + 1)[1:]
    recs = mda_records(n, maxdim, ep)
    p_recs, a_recs = DA(recs.flatten()); p_ep, a_ep = DA(ep)
    segid = f"type{typ}".ljust(8)
    if typ == 1:
        L.spkw01_(I(h), I(body), I(10), S("J2000"), D(T0), D(T1), S(segid), I(n), p_recs, p_ep,
                  ctypes.c_int(5), ctypes.c_int(len(segid)))
    else:
        L.spkw21_(I(h), I(body), I(10), S("J2000"), D(T0), D(T1), S(segid), I(n), I(4 * maxdim + 11), p_recs, p_ep,
                  ctypes.c_int(5), ctypes.c_int(len(segid)))
    check_err()

# Type 19: two intervals with different subtypes, and a single-interval one with select-first
def type19(body, sellst):
    bnds = np.array([T0, -1e7, 4e7, T1])
    subs = [0, 1, 2]; degs = [7, 9, 3]
    pkts_all, eps_all, npk = [], [], []
    for i in range(3):
        m = 57 + 13 * i
        ep = np.sort(rng.uniform(bnds[i], bnds[i + 1], m)); ep[0] = bnds[i]; ep[-1] = bnds[i + 1]
        st = kepler_states(ep)
        if subs[i] == 0:
            pk = np.hstack([st, st[:, 3:] + 1e-4, rng.normal(size=(m, 3)) * 1e-9])
        else:
            pk = st
        pkts_all.extend(pk.flatten()); eps_all.extend(ep); npk.append(m)
    segid = "type19".ljust(8)
    p_npk, a1 = IA(npk); p_sub, a2 = IA(subs); p_deg, a3 = IA(degs)
    p_pk, a4 = DA(pkts_all); p_ep, a5 = DA(eps_all); p_iv, a6 = DA(bnds)
    L.spkw19_(I(h), I(body), I(10), S("J2000"), D(T0), D(T1), S(segid), I(3), p_npk, p_sub, p_deg,
              p_pk, p_ep, p_iv, ctypes.byref(c_int(1 if sellst else 0)), ctypes.c_int(5), ctypes.c_int(len(segid)))
    check_err()
type19(1190, True)
type19(1191, False)

sp.spkcls(h)
check_err()

# --------------------------------------------------------------------------------------------
# Binary PCK: types 2 and 20 for class 3000 (ITRF93) relative to ECLIPJ2000, and type 2 for a
# kernel-defined body-fixed frame.
# --------------------------------------------------------------------------------------------
pck = fresh(f"{OUT}/synthetic.bpc")
h = sp.pckopn(pck, "synthetic pck", 0)
n = 30; deg = 8; intlen = (T1 - T0) / n
co = rng.normal(size=(n, 3, deg + 1)) * np.array([1.0 / (k + 1) ** 2 for k in range(deg + 1)])
co[:, 2, 1] += 7.29e-5 * intlen / 2  # spin
sp.pckw02(h, 3000, "ECLIPJ2000", T0, T1, "itrf93-like", intlen, n, deg, co.flatten(), T0)
# PCK type 20 for class 3001 relative to J2000
nt = 7; cdata = []
for r in range(n):
    for c in range(3):
        cdata.extend(list(rng.normal(size=nt) * 1e-6)); cdata.append(rng.normal())
segid = "pck20".ljust(8)
p_c, a_c = DA(cdata)
L.pckw20_(I(h), I(3001), S("J2000"), D(T0), D(T1), S(segid), D(intlen / 86400.0), I(n), I(nt - 1), p_c,
          D(1.0), D(1.0), D(2451545.0 - 1158.0), D(T0 / 86400.0 + 1158.0), ctypes.c_int(5), ctypes.c_int(len(segid)))
check_err()
sp.pckcls(h)
check_err()

# PCK type 3: copy the raw generic-segment data of the SPK type 14 segment into a PCK segment.
hs = sp.dafopr(spk); sp.dafbfs(hs); raw = None
while sp.daffna():
    s = sp.dafgs(); dc, ic = sp.dafus(s, 2, 6)
    if ic[3] == 14:
        raw = sp.dafgda(hs, int(ic[4]), int(ic[5])); t14 = (dc[0], dc[1])
sp.dafcls(hs)
pck3 = fresh(f"{OUT}/synthetic3.bpc")
h3 = sp.pckopn(pck3, "pck3", 0)
# Build the summary manually: dafbna takes packed summary doubles
pk = np.zeros(2 + 3, dtype=np.float64)
pk[0], pk[1] = t14
ints = np.array([3002, 1, 3, 0, 0, 0], dtype=np.int32)
packed = np.concatenate([pk[:2], np.frombuffer(ints.tobytes(), dtype=np.float64)])
p_sum, a_sum = DA(packed)
name = "pck type 3".ljust(40)
L.dafbna_(I(h3), p_sum, S(name), ctypes.c_int(len(name)))
p_raw, a_raw = DA(raw)
L.dafada_(p_raw, I(len(raw)))
L.dafena_()
check_err()
sp.pckcls(h3)

# --------------------------------------------------------------------------------------------
# Text kernels: frames, names, IAU model, and a meta-kernel.
# --------------------------------------------------------------------------------------------
open(f"{OUT}/frames.tf", "w").write(r"""
Synthetic frame kernel
\begindata
FRAME_MYTK_MAT      = 1400001
FRAME_1400001_NAME  = 'MYTK_MAT'
FRAME_1400001_CLASS = 4
FRAME_1400001_CLASS_ID = 1400001
FRAME_1400001_CENTER = 399
TKFRAME_1400001_RELATIVE = 'ITRF93'
TKFRAME_1400001_SPEC = 'MATRIX'
TKFRAME_1400001_MATRIX = ( 0.4 0.6 0.0 -0.6 0.4 0.0 0.0 0.0 1.0 )

FRAME_MYTK_ANG      = 1400002
FRAME_1400002_NAME  = 'MYTK_ANG'
FRAME_1400002_CLASS = 4
FRAME_1400002_CLASS_ID = 1400002
FRAME_1400002_CENTER = 'EARTH'
TKFRAME_MYTK_ANG_RELATIVE = 'MYTK_MAT'
TKFRAME_MYTK_ANG_SPEC = 'ANGLES'
TKFRAME_MYTK_ANG_ANGLES = ( -71.3, -42.1, 180.0 )
TKFRAME_MYTK_ANG_AXES = ( 3, 2, 3 )
TKFRAME_MYTK_ANG_UNITS = 'DEGREES'

FRAME_MYTK_Q      = 1400003
FRAME_1400003_NAME  = 'MYTK_Q'
FRAME_1400003_CLASS = 4
FRAME_1400003_CLASS_ID = 1400003
FRAME_1400003_CENTER = 0
TKFRAME_1400003_RELATIVE = 'GALACTIC'
TKFRAME_1400003_SPEC = 'QUATERNION'
TKFRAME_1400003_Q = ( 0.5 0.5 -0.5 0.5 )

FRAME_MYPCK20   = 1400004
FRAME_1400004_NAME = 'MYPCK20'
FRAME_1400004_CLASS = 2
FRAME_1400004_CLASS_ID = 3001
FRAME_1400004_CENTER = 399

FRAME_MYPCK3   = 1400005
FRAME_1400005_NAME = 'MYPCK3'
FRAME_1400005_CLASS = 2
FRAME_1400005_CLASS_ID = 3002
FRAME_1400005_CENTER = 399

NAIF_BODY_NAME += ( 'MY ROCK', 'OTHER_ROCK' )
NAIF_BODY_CODE += ( 1002, 1003 )
\begintext
""")

open(f"{OUT}/iau.tpc", "w").write(r"""
\begindata
BODY499_POLE_RA  = ( 317.269202  -0.10927547   0.0 )
BODY499_POLE_DEC = (  54.432516  -0.05827105   0.0 )
BODY499_PM       = ( 176.049863 +350.891982443297  0. )
BODY499_NUT_PREC_RA  = ( 0 0 0 0 0 0 0 0 0 0   0.000068  0.000238  0.000052  0.000009  0.419057 )
BODY499_NUT_PREC_DEC = ( 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0   0.000051  0.000141  0.000031  0.000005  1.591274 )
BODY499_NUT_PREC_PM  = ( 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0   0.000145  0.000157  0.000040  0.000001  0.000001  0.584542 )
BODY4_NUT_PREC_ANGLES  = (
   190.72646643  15917.10818695  21.46892470  31834.27934054  332.86082793  19139.89694742
   394.93256437  38280.79631835  189.63271560  41215158.18420050  121.46893664  660.22803474
   231.05028581  660.99123540  251.37314025  1320.50145245  217.98635955  38279.96125550
   196.19729402  19139.83628608  198.991226  19139.4819985  226.292679  38280.8511281
   249.663391  57420.7251593  266.183510  76560.6367950  79.398797  0.5042615
   122.433576  19139.9407476  43.058401  38280.8753272  57.663379  57420.7517205
   79.476401  76560.6495004  166.325722  0.5042615  129.071773  19140.0328244
   36.352167  38281.0473591  56.668646  57420.9295360  67.364003  76560.2552215
   104.792680  95700.4387578  95.391654  0.5042615 )
BODY399_GM = 398600.435436
BODY10_GM  = 132712440041.939400
\begintext
""")

# Meta-kernel loading the synthetic files with a path symbol
open(f"{OUT}/meta.tm", "w").write(f"""
\\begindata
PATH_VALUES  = ( '{OUT}' )
PATH_SYMBOLS = ( 'K' )
KERNELS_TO_LOAD = ( '$K/synthetic.bsp', '$K/synth+'
                    'etic.bpc', '$K/synthetic3.bpc', '$K/frames.tf', '$K/iau.tpc' )
\\begintext
""")
print("generated")
