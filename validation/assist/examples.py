"""Reproduce ASSIST's published examples with spacerocks and ASSIST side by side.

    pip install assist matplotlib
    cd py_spacerocks && maturin develop --release
    SPACEROCKS_KERNELS=/path/to/kernels python validation/assist/examples.py [outdir]

`SPACEROCKS_KERNELS` must hold `de440.bsp` and the full `sb441-n16.bsp` (the same files ASSIST's
examples use). The JPL reference data (`apophis_sb_daily_v2.txt` and the `*_sb.txt` files) come
from ASSIST's repository (jupyter_examples/); set `ASSIST_EXAMPLES` to that directory.

Every example runs in both codes with the same kernels, initial conditions and force model
(ASSIST's defaults, plus Einstein-Infeld-Hoffmann GR from all 11 major bodies where the ASSIST
example turns it on). Results are written to `results.json` in the output directory; `plots.py`
turns them into figures.
"""
import json
import os
import sys
import time

import numpy as np
import rebound
import assist

from spacerocks import SpaceRock, RockCollection
from spacerocks.spice import SpiceKernel
from spacerocks.time import Time
from spacerocks.assist import SpiceSimulation, Force

AU_M = 149597870700.0
AU_KM = AU_M / 1000.0

KDIR = os.environ["SPACEROCKS_KERNELS"]
EXDIR = os.environ.get("ASSIST_EXAMPLES", os.path.join(os.path.dirname(__file__), "data"))
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(__file__), "out")
os.makedirs(OUT, exist_ok=True)

DE = os.path.join(KDIR, "de440.bsp")
SB = os.path.join(KDIR, "sb441-n16.bsp")
ephem = assist.Ephem(DE, SB)
JREF = ephem.jd_ref
kernel = SpiceKernel()
kernel.load(DE)
kernel.load(SB)

# Initial conditions from ASSIST's notebooks and unit tests (barycentric unless noted).
HOLMAN_T0 = 2458849.5
HOLMAN = [3.338875348598862E+00, -9.176518412197102E-01, -5.038590741719294E-01,
          2.805663364339457E-03, 7.550408665778840E-03, 2.980028207875623E-03]
FARNOCCHIA_T0 = 2458849.5
FARNOCCHIA = [9.572786624350362E-01, -2.101947659454845E+00, -7.418444809938682E-01,
              9.715385994801106E-03, 6.153629531316274E-03, 1.549521077070705E-03]
APOPHIS_NG = [4.999999873689E-13, -2.901085508711E-14, 0.0]


def timed(f, repeat=3):
    """Run f() `repeat` times; return (last result, best wall time in s)."""
    best = np.inf
    for _ in range(repeat):
        t = time.perf_counter()
        r = f()
        best = min(best, time.perf_counter() - t)
    return r, best


def sun_state(jd):
    s = ephem.get_particle("sun", jd - JREF)
    return np.array([s.x, s.y, s.z, s.vx, s.vy, s.vz])


# ---------------------------------------------------------------------------------------------
# Thin wrappers so both codes are driven the same way.

class AssistRun:
    def __init__(self, t0, states, gr_sources=1, nongrav=None, min_dt=None, dt=None, variations=()):
        self.sim = rebound.Simulation()
        for s in states:
            self.sim.add(x=s[0], y=s[1], z=s[2], vx=s[3], vy=s[4], vz=s[5])
        self.vars = []
        for (tp, comp) in variations:
            v = self.sim.add_variation(testparticle=tp)
            setattr(v.particles[0], comp, 1.0)
            self.vars.append(v)
        self.sim.t = t0 - JREF
        if dt is not None:
            self.sim.dt = dt
        if min_dt is not None:
            self.sim.ri_ias15.min_dt = min_dt
        self.ex = assist.Extras(self.sim, ephem)
        self.ex.gr_eih_sources = gr_sources
        if nongrav is not None:
            n = self.sim.N
            params = np.zeros(3 * n)
            params[:3] = nongrav
            self.ex.particle_params = params

    def at(self, jd):
        self.ex.integrate_or_interpolate(jd - JREF)
        return np.array([list(p.xyz) + list(p.vxyz) for p in self.sim.particles[: self.sim.N - self.sim.N_var]])

    def var_states(self):
        return np.array([list(v.particles[0].xyz) + list(v.particles[0].vxyz) for v in self.vars])


class SpacerocksRun:
    mode = None  # step-size criterion; None keeps spacerocks' default (PRS23)
    summation = os.environ.get("SUMMATION")  # None keeps spacerocks' default (Kahan)

    def __init__(self, t0, states, gr_sources=1, nongrav=None, min_dt=None, variations=()):
        self.t0 = Time(t0, "tdb", "jd")
        self.sim = SpiceSimulation.horizons(self.t0, kernel)
        if self.mode is not None:
            self.sim.adaptive_mode = self.mode
        if self.summation is not None:
            self.sim.summation = self.summation
        if gr_sources != 1:
            self.sim.set_forces([Force.nongravitational(), Force.earth_harmonics(kernel), Force.solar_j2(kernel),
                                 Force.gr_eih(kernel, sources=gr_sources), Force.newtonian_gravity()])
        if min_dt is not None:
            self.sim.min_timestep = min_dt
        for i, s in enumerate(states):
            r = SpaceRock.from_xyz(f"p{i}", *s, self.t0, "J2000", "SSB")
            if nongrav is not None and i == 0:
                r.set_nongrav(*nongrav)
            self.sim.add(r)
        for (tp, comp) in variations:
            self.sim.add_variation(comp, f"p{tp}")

    def at(self, jd):
        self.sim.integrate(Time(jd, "tdb", "jd"), kernel)
        return np.array([[p.x, p.y, p.z, p.vx, p.vy, p.vz] for p in self.sim.particles])

    def var_states(self):
        return self.sim.variational_particles


class SpacerocksGlobalRun(SpacerocksRun):
    mode = "global"


# ASSIST (its default: REBOUND's "global" step criterion), spacerocks with its default
# (PRS23) criterion, and spacerocks with ASSIST's criterion.
CODES = {"ASSIST": AssistRun, "spacerocks": SpacerocksRun, "spacerocks (global)": SpacerocksGlobalRun}
results = {}


def record(name, **kw):
    results[name] = kw
    print(f"== {name}")
    for k, v in kw.items():
        if isinstance(v, dict):
            print(f"   {k}: " + ", ".join(f"{a}={b:.4g}" if isinstance(b, float) else f"{a}={b}" for a, b in v.items() if not isinstance(b, (list, dict))))
        elif not isinstance(v, list):
            print(f"   {k}: {v}")


# ---------------------------------------------------------------------------------------------
# 1. Holman, 30 days, against JPL Horizons (ASSIST unit test holman_spk; 10 cm tolerance).
def ex_holman_30d():
    t0, t1 = JREF + 8416.5, JREF + 8446.5
    x0 = [-2.724183384883979E+00, -3.523994546329214E-02, 9.036596202793466E-02,
          -1.374545432301129E-04, -1.027075301472321E-02, -4.195690627695180E-03]
    jpl = np.array([-2.710320457933958E+00, -3.424507930535848E-01, -3.582442972611413E-02,
                    1.059255302926290E-03, -1.018748422976772E-02, -4.207712906489264E-03])
    out = {}
    for code, cls in CODES.items():
        s, dt = timed(lambda: cls(t0, [x0]).at(t1)[0], repeat=5)
        out[code] = {"err_m": float(np.linalg.norm(s[:3] - jpl[:3]) * AU_M), "time_s": dt}
    record("holman_30d", **out)


# 2. ASSIST "Getting started": (84100) Farnocchia for 10,000 days.
def ex_getting_started():
    t0 = FARNOCCHIA_T0
    out, final = {}, {}
    for code, cls in CODES.items():
        s, dt = timed(lambda: cls(t0, [FARNOCCHIA]).at(t0 + 10_000)[0])
        final[code] = s
        out[code] = {"time_s": dt}
    for c in CODES:
        out[c]["difference_from_assist_m"] = float(np.linalg.norm(final["ASSIST"][:3] - final[c][:3]) * AU_M)
    record("getting_started", **out)


# 3. Apophis 2029 close approach against JPL's orbit solution (ASSIST Apophis notebook / paper Fig. 3).
APOPHIS_T0 = 2.4621385359989386E+06
APOPHIS_T1 = 2.4625030372426095E+06
APOPHIS_HELIO = [-5.5946538550488512E-01, 8.5647564757574512E-01, 3.0415066217102493E-01,
                 -1.3818324735921638E-02, -6.0088275597939191E-03, -2.5805044631309632E-03]
APOPHIS_FINAL_HELIO = [1.7028330901729331E-02, 1.2193934090901304E+00, 4.7823589236374386E-01,
                       -1.3536187639388663E-02, 5.3200999989786943E-04, -1.6648346717629861E-05]


def ex_apophis():
    x0 = np.array(APOPHIS_HELIO) + sun_state(APOPHIS_T0)
    jpl = np.array(APOPHIS_FINAL_HELIO) + sun_state(APOPHIS_T1)
    times = np.linspace(APOPHIS_T0, APOPHIS_T1, 10000)
    out, tracks = {}, {}
    for code, cls in CODES.items():
        def run():
            r = cls(APOPHIS_T0, [x0], gr_sources=11, nongrav=APOPHIS_NG, min_dt=0.001)
            return np.array([r.at(t)[0] for t in times])
        track, dt = timed(run)
        tracks[code] = track
        out[code] = {"err_m": float(np.linalg.norm(track[-1, :3] - jpl[:3]) * AU_M), "time_s": dt}
    earth = np.array([ephem.get_particle("earth", t - JREF).xyz for t in times])
    geo = {c: ((tracks[c][:, :3] - earth) * AU_KM).tolist() for c in CODES}
    diff = {c: (np.linalg.norm(tracks["ASSIST"][:, :3] - tracks[c][:, :3], axis=1) * AU_M).tolist() for c in CODES if c != "ASSIST"}
    i = int(np.argmin(np.linalg.norm(tracks["ASSIST"][:, :3] - earth, axis=1)))
    out["closest_approach"] = {"jd": float(times[i]), "geocentric_km": float(np.linalg.norm(tracks["ASSIST"][i, :3] - earth[i]) * AU_KM)}
    out["series"] = {"t": (times - APOPHIS_T0).tolist(), "diff_from_assist_m": diff,
                     "geo_km": {c: g[::10] for c, g in geo.items()}}
    record("apophis_2029", **out)

    # Step sizes through the encounter (paper Fig. 4).
    steps = {}
    r = AssistRun(APOPHIS_T0, [x0], gr_sources=11, nongrav=APOPHIS_NG, min_dt=0.001, dt=1.0)
    ts, dts = [], []
    while r.sim.t < APOPHIS_T1 - JREF:
        r.sim.step()
        ts.append(r.sim.t + JREF - APOPHIS_T0)
        dts.append(r.sim.dt_last_done)
    steps["ASSIST"] = [ts, dts]
    for c in ["spacerocks", "spacerocks (global)"]:
        r = CODES[c](APOPHIS_T0, [x0], gr_sources=11, nongrav=APOPHIS_NG, min_dt=0.001)
        ts, dts = [], []
        while r.sim.step_epoch < APOPHIS_T1:
            r.sim.step(kernel)
            ts.append(r.sim.step_epoch - APOPHIS_T0)
            dts.append(r.sim.last_timestep)
        steps[c] = [ts, dts]
    results["apophis_steps"] = steps
    print("   steps: " + ", ".join(f"{c} {len(v[0])}" for c, v in steps.items()))


# 4. Apophis daily positions against JPL's small-body integrator (paper Fig. 5).
def ex_apophis_daily():
    t0 = 2462137.5
    helio = [-5.450770368702937E-01, 8.625884531011220E-01, 3.067841508125688E-01,
             -1.395794822904029E-02, -5.791529394125580E-03, -2.503280483635055E-03]
    x0 = np.array(helio) + sun_state(t0)
    rows = []
    with open(os.path.join(EXDIR, "apophis_sb_daily_v2.txt")) as f:
        for _ in range(3):
            f.readline()
        for line in f:
            parts = line.split()
            rows.append((float(parts[3]), np.array(parts[-3:], dtype=float)))
    times = np.array([r[0] for r in rows])
    ref = np.array([v / AU_KM + sun_state(t)[:3] for t, v in rows])
    out, series = {}, {"t": (times - t0).tolist()}
    for code, cls in CODES.items():
        def run():
            r = cls(t0, [x0], gr_sources=11, nongrav=APOPHIS_NG, min_dt=0.001)
            return np.array([r.at(t)[0, :3] for t in times])
        pos, dt = timed(run)
        err = np.linalg.norm(pos - ref, axis=1) * AU_M
        series[code] = err.tolist()
        out[code] = {"max_err_m": float(err.max()), "final_err_m": float(err[-1]), "time_s": dt}
    out["series"] = series
    out["n_epochs"] = len(times)
    record("apophis_daily", **out)


# 5. JPL small-body integrator positions at 10 to 100,000 days (the *_sb.txt files).
def ex_sb_long():
    cases = {
        "Holman": (HOLMAN_T0, HOLMAN, None, "holman_sb.txt"),
        "Farnocchia": (FARNOCCHIA_T0, FARNOCCHIA, None, "farnocchia_sb.txt"),
        "Apophis": (2462137.5, list(np.array([-5.450770368702937E-01, 8.625884531011220E-01, 3.067841508125688E-01,
                                              -1.395794822904029E-02, -5.791529394125580E-03, -2.503280483635055E-03])
                                    + sun_state(2462137.5)), APOPHIS_NG, "apophis_sb.txt"),
    }
    out = {}
    for name, (t0, x0, ng, fname) in cases.items():
        rows = []
        with open(os.path.join(EXDIR, fname)) as f:
            for line in f:
                p = line.split()
                if len(p) == 7 and p[1] == "d:":
                    rows.append((float(p[0]), float(p[3]), np.array(p[4:7], dtype=float)))
        entry = {"span_days": [r[0] for r in rows]}
        for code, cls in CODES.items():
            def run():
                r = cls(t0, [x0], gr_sources=11, nongrav=ng)
                return [r.at(jd)[0, :3] for (_, jd, _) in rows]
            pos, dt = timed(run, repeat=1)
            errs = [float(np.linalg.norm(p - (v / AU_KM + sun_state(jd)[:3])) * AU_M) for p, (_, jd, v) in zip(pos, rows)]
            entry[code] = {"err_m": errs, "time_s": dt}
        out[name] = entry
        print(f"== sb {name}")
        for c in CODES:
            print(f"   {c:22s} " + ", ".join(f"{int(sp)} d: {e:.3g} m" for sp, e in zip(entry["span_days"], entry[c]["err_m"])))
    results["sb_long"] = out


# 6. Round trip: integrate Holman forward and back (paper Fig. 2).
def ex_round_trip():
    spans = np.logspace(0, 5, 21)
    out = {"span_days": spans.tolist()}
    for code in CODES:
        errs = []
        t_start = time.perf_counter()
        for T in spans:
            if code == "ASSIST":
                r = AssistRun(HOLMAN_T0, [HOLMAN], gr_sources=11, dt=20.0)
                r.sim.integrate(HOLMAN_T0 + T - JREF)
                r.sim.dt *= -1
                r.sim.integrate(HOLMAN_T0 - JREF)
                back = np.array(r.sim.particles[0].xyz)
            else:
                r = CODES[code](HOLMAN_T0, [HOLMAN], gr_sources=11)
                r.at(HOLMAN_T0 + T)
                back = r.at(HOLMAN_T0)[0, :3]
            errs.append(float(np.linalg.norm(back - np.array(HOLMAN[:3])) * AU_M))
        out[code] = {"err_m": errs, "time_s": time.perf_counter() - t_start}
    record("round_trip", **out)


# 7. Variational equations vs a shadow particle (ASSIST VariationalEquations notebook / Fig. 1).
def ex_variational():
    scale = 1e-8
    shadow = list(np.array(HOLMAN) + np.array([scale, 0, 0, 0, 0, 0]))
    times = np.linspace(HOLMAN_T0, HOLMAN_T0 + 10_000, 500)
    out = {"t": (times - HOLMAN_T0).tolist()}
    for code, cls in CODES.items():
        def run():
            r = cls(HOLMAN_T0, [HOLMAN, shadow], variations=[(0, "x")])
            d, v = [], []
            for t in times:
                s = r.at(t)
                d.append(s[1, :3] - s[0, :3])
                v.append(r.var_states()[0][:3] * scale)
            return np.array(d), np.array(v)
        (d, v), dt = timed(run, repeat=1)
        out[code] = {"shadow_km": (d[:, 0] * AU_KM).tolist(),
                     "residual_mm": (np.linalg.norm(d - v, axis=1) * AU_M * 1e3).tolist(), "time_s": dt}
    record("variational", **out)


# 8. (5303) Parijskij's encounter with Ceres (paper Fig. 6).
def ex_ceres_5303():
    t0 = 2449718.5
    x0 = [-2.232847879711731E+00, 1.574146331186095E+00, 8.329414259670296E-01,
          -6.247432571575564E-03, -7.431073424167182E-03, -3.231725223736132E-03]
    times = np.linspace(t0, t0 + 3653.0, 3653)
    ceres = np.array([ephem.get_particle("ceres", t - JREF).xyz for t in times])
    out = {"t": (times - t0).tolist()}
    for code, cls in CODES.items():
        def run():
            r = cls(t0, [x0])
            return np.array([r.at(t)[0, :3] for t in times])
        pos, dt = timed(run)
        sep = np.linalg.norm(pos - ceres, axis=1)
        i = int(np.argmin(sep))
        out[code] = {"closest_jd": float(times[i]), "closest_km": float(sep[i] * AU_KM), "time_s": dt,
                     "separation_au": sep.tolist(), "pos": pos.tolist()}
    for c in CODES:
        out[c]["max_diff_from_assist_m"] = float(np.max(np.linalg.norm(np.array(out["ASSIST"]["pos"]) - np.array(out[c]["pos"]), axis=1)) * AU_M)
    for code in CODES:
        del out[code]["pos"]
    record("ceres_5303", **out)


# 9. Convergence with the tolerance epsilon through two encounters (Newtonian gravity only, no
#    step floor): (5303) with Ceres and Apophis with the Earth. Positions at the end are compared
#    with the tightest converged run.
def ex_convergence():
    out = {}
    cases = {
        "ceres_5303": (2449718.5, 2449718.5 + 3653.0,
                       [-2.232847879711731E+00, 1.574146331186095E+00, 8.329414259670296E-01,
                        -6.247432571575564E-03, -7.431073424167182E-03, -3.231725223736132E-03]),
        "apophis_2029": (APOPHIS_T0, APOPHIS_T1, list(np.array(APOPHIS_HELIO) + sun_state(APOPHIS_T0))),
    }
    eps_list = [1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 1e-13]
    for name, (t0, t1, x0) in cases.items():
        entry = {"eps": eps_list}
        finals = {}
        for code in CODES:
            finals[code] = []
            for eps in eps_list:
                if code != "spacerocks" and eps < 1e-11:
                    # The last-term ("global") criterion is round-off limited below ~1e-11: steps
                    # shrink without bound (in ASSIST too), so those runs are skipped.
                    finals[code].append(None)
                    continue
                r = CODES[code](t0, [x0])
                if code == "ASSIST":
                    r.ex.forces = ["SUN", "PLANETS", "ASTEROIDS"]
                    r.sim.ri_ias15.epsilon = eps
                else:
                    r.sim.newtonian_only()
                    r.sim.epsilon = eps
                tt = time.perf_counter()
                finals[code].append(r.at(t1)[0, :3])
                print(f"      {name} {code} eps={eps:g}: {time.perf_counter() - tt:.2f} s", flush=True)
        # Reference: the run that all converged configurations approach (ASSIST's tightest
        # epsilon for Ceres; for Apophis, spacerocks PRS23 at 1e-13, which is converged to cm).
        ref = finals["ASSIST"][3] if name == "ceres_5303" else finals["spacerocks"][-1]
        for code in CODES:
            entry[code] = [None if f is None else float(np.linalg.norm(f - ref) * AU_M) for f in finals[code]]
        out[name] = entry
        print(f"== convergence {name}")
        for code in CODES:
            print(f"   {code:22s} " + ", ".join("-" if e is None else f"{e:.3g}" for e in entry[code]))
    results["convergence"] = out


# 10. Throughput: ASSIST's benchmark (N copies of Holman, 10 years, no outputs), N = 1 .. 1000.
def ex_benchmark():
    t0 = JREF + 8416.5
    x0 = [-2.724183384883979E+00, -3.523994546329214E-02, 9.036596202793466E-02,
          -1.374545432301129E-04, -1.027075301472321E-02, -4.195690627695180E-03]
    out = {"n": [1, 10, 100, 1000]}
    for code in ["ASSIST", "spacerocks", "spacerocks (global)", "spacerocks batch"]:
        out[code] = []
    for n in out["n"]:
        states = [[x0[0] + i * 1e-10] + x0[1:] for i in range(n)]
        _, ta = timed(lambda: AssistRun(t0, states).at(t0 + 3652.5), repeat=3 if n < 1000 else 1)
        _, ts = timed(lambda: SpacerocksRun(t0, states).at(t0 + 3652.5), repeat=3 if n < 1000 else 1)
        _, tg = timed(lambda: SpacerocksGlobalRun(t0, states).at(t0 + 3652.5), repeat=3 if n < 1000 else 1)

        def batch():
            rc = RockCollection()
            for i, s in enumerate(states):
                rc.add(SpaceRock.from_xyz(f"p{i}", *s, Time(t0, "tdb", "jd"), "J2000", "SSB"))
            rc.propagate(Time(t0 + 3652.5, "tdb", "jd"), kernel)
        _, tb = timed(batch, repeat=3 if n < 1000 else 1)
        out["ASSIST"].append(ta)
        out["spacerocks"].append(ts)
        out["spacerocks (global)"].append(tg)
        out["spacerocks batch"].append(tb)
        print(f"   N={n}: ASSIST {ta*1e3:.1f} ms, spacerocks {ts*1e3:.1f} ms, spacerocks (global) {tg*1e3:.1f} ms, spacerocks batch {tb*1e3:.1f} ms")
    results["benchmark"] = out


# 12. Both codes with the PRS23 step rule (REBOUND's default, available in ASSIST as
#     sim.ri_ias15.adaptive_mode = "prs23"): (5303)-Ceres encounter (Newtonian, final position
#     vs ASSIST's global rule at epsilon 1e-11) and the Apophis flyby (full model, vs JPL).
def ex_prs23_crosscheck():
    out = {"eps": [1e-9, 1e-10, 1e-11, 1e-12, 1e-13]}
    t0 = 2449718.5
    x0 = [-2.232847879711731E+00, 1.574146331186095E+00, 8.329414259670296E-01,
          -6.247432571575564E-03, -7.431073424167182E-03, -3.231725223736132E-03]
    ref = AssistRun(t0, [x0])
    ref.sim.ri_ias15.epsilon = 1e-11
    ref.ex.forces = ["SUN", "PLANETS", "ASTEROIDS"]
    pref = ref.at(t0 + 3653.0)[0, :3]
    xa = np.array(APOPHIS_HELIO) + sun_state(APOPHIS_T0)
    jpl = (np.array(APOPHIS_FINAL_HELIO) + sun_state(APOPHIS_T1))[:3]
    for name in ["ceres_5303", "apophis_2029"]:
        entry = {"ASSIST": [], "spacerocks": [], "ASSIST vs spacerocks": []}
        for eps in out["eps"]:
            fin = {}
            for code in ["ASSIST", "spacerocks"]:
                if name == "ceres_5303":
                    r = CODES[code](t0, [x0])
                    if code == "ASSIST":
                        r.ex.forces = ["SUN", "PLANETS", "ASTEROIDS"]
                    else:
                        r.sim.newtonian_only()
                    t1 = t0 + 3653.0
                else:
                    r = CODES[code](APOPHIS_T0, [xa], gr_sources=11, nongrav=APOPHIS_NG, min_dt=0.001)
                    t1 = APOPHIS_T1
                if code == "ASSIST":
                    r.sim.ri_ias15.adaptive_mode = "prs23"
                    r.sim.ri_ias15.epsilon = eps
                else:
                    r.sim.epsilon = eps
                fin[code] = r.at(t1)[0, :3]
            target = pref if name == "ceres_5303" else jpl
            for code in ["ASSIST", "spacerocks"]:
                entry[code].append(float(np.linalg.norm(fin[code] - target) * AU_M))
            entry["ASSIST vs spacerocks"].append(float(np.linalg.norm(fin["ASSIST"] - fin["spacerocks"]) * AU_M))
        out[name] = entry
        print(f"== prs23 {name} (" + ("vs converged global" if name == "ceres_5303" else "vs JPL") + ")")
        for k, v in entry.items():
            print(f"   {k:22s} " + ", ".join(f"{x:.3g}" for x in v))
    # Against JPL's small-body integrator at epsilon 1e-9, as in ex_sb_long.
    for name, t0, x0, fname in [("Holman", HOLMAN_T0, HOLMAN, "holman_sb.txt"),
                                ("Farnocchia", FARNOCCHIA_T0, FARNOCCHIA, "farnocchia_sb.txt")]:
        rows = []
        with open(os.path.join(EXDIR, fname)) as f:
            for line in f:
                p = line.split()
                if len(p) == 7 and p[1] == "d:":
                    rows.append((float(p[0]), float(p[3]), np.array(p[4:7], dtype=float)))
        entry = {"span_days": [r[0] for r in rows]}
        for code in ["ASSIST", "spacerocks"]:
            r = CODES[code](t0, [x0], gr_sources=11)
            if code == "ASSIST":
                r.sim.ri_ias15.adaptive_mode = "prs23"
            entry[code] = [float(np.linalg.norm(r.at(jd)[0, :3] - (v / AU_KM + sun_state(jd)[:3])) * AU_M) for (_, jd, v) in rows]
        out["sb_" + name] = entry
        print(f"== prs23 sb {name}: " + "; ".join(f"{c} " + ", ".join(f"{x:.3g}" for x in entry[c]) for c in ["ASSIST", "spacerocks"]))
    results["prs23_crosscheck"] = out


# 11. Compensated summation (spacerocks only): round trips of 64 Holman clones (offset by ~1e-9
#     AU so each rounds differently) with plain sums, Kahan-compensated positions and velocities
#     (the default) and REBOUND's full scheme; and the cost per step (1000 clones, 10 years).
def ex_summation():
    rng = np.random.default_rng(1)
    states = [list(np.array(HOLMAN) + np.r_[rng.normal(0, 1e-9, 3), np.zeros(3)]) for _ in range(64)]
    x0 = np.array(states)[:, :3]
    modes = ["none", "kahan", "full"]
    spans = [100.0, 1000.0, 10000.0, 100000.0]
    out = {"spans": spans, "modes": modes, "configs": {}}
    for rule, eps in [("global", 1e-9), ("global", 1e-11), ("prs23", 1e-13)]:
        cfg = {}
        for m in modes:
            med, p90 = [], []
            for T in spans:
                r = SpacerocksRun(HOLMAN_T0, states, gr_sources=11)
                r.sim.adaptive_mode = rule
                r.sim.summation = m
                r.sim.epsilon = eps
                r.at(HOLMAN_T0 + T)
                e = np.linalg.norm(r.at(HOLMAN_T0)[:, :3] - x0, axis=1) * AU_M
                med.append(float(np.median(e)))
                p90.append(float(np.percentile(e, 90)))
            cfg[m] = {"median_m": med, "p90_m": p90}
            print(f"   {rule} eps={eps:g} {m:6s} " + ", ".join(f"{int(T)} d: {a:.2e} m" for T, a in zip(spans, med)), flush=True)
        out["configs"][f"{rule}, epsilon {eps:g}"] = cfg
    t0 = JREF + 8416.5
    xh = [-2.724183384883979E+00, -3.523994546329214E-02, 9.036596202793466E-02,
          -1.374545432301129E-04, -1.027075301472321E-02, -4.195690627695180E-03]
    clones = [[xh[0] + i * 1e-10] + xh[1:] for i in range(1000)]
    out["time_s"] = {}
    for m in modes:
        def run():
            r = SpacerocksRun(t0, clones)
            r.sim.summation = m
            r.at(t0 + 3652.5)
        _, out["time_s"][m] = timed(run, repeat=3)
    print("   1000 clones, 10 years: " + ", ".join(f"{m} {t:.3f} s" for m, t in out["time_s"].items()))
    results["summation"] = out


if __name__ == "__main__":
    only = os.environ.get("ONLY")
    if only and os.path.exists(os.path.join(OUT, "results.json")):
        results.update(json.load(open(os.path.join(OUT, "results.json"))))
    for f in [ex_holman_30d, ex_getting_started, ex_apophis, ex_apophis_daily, ex_sb_long,
              ex_round_trip, ex_variational, ex_ceres_5303, ex_convergence, ex_benchmark, ex_summation, ex_prs23_crosscheck]:
        if only and f.__name__ not in only.split(","):
            continue
        f()
        with open(os.path.join(OUT, "results.json"), "w") as fh:
            json.dump(results, fh)
