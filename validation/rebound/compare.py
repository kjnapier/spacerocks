"""Compare spacerocks' WisdomHolman, Trace and IAS15 integrators with REBOUND's.

Each scenario builds initial conditions with REBOUND, integrates them with both codes (the
spacerocks side through `examples/nbody_compare.rs`) and with IAS15 in both codes as the
reference, and prints JSON with final-state errors, energy errors and wall times.

    cargo build --release --example nbody_compare
    python3 validation/rebound/compare.py > validation/rebound/results.json
"""
import json
import math
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import rebound

G = 2.9591220828559104e-4  # spacerocks GRAVITATIONAL_CONSTANT (AU^3 / Msun / day^2)
BIN = Path(__file__).resolve().parents[2] / "target/release/examples/nbody_compare"
REPEATS = 3

GIANTS = [  # mass (Msun), a (AU), e, inc (deg), Omega, omega, M (deg); approximate J2000 values
    ("jupiter", 9.547919e-4, 5.2026, 0.0484, 1.304, 100.5, 273.9, 20.0),
    ("saturn", 2.858860e-4, 9.5549, 0.0539, 2.486, 113.7, 339.4, 317.0),
    ("uranus", 4.366244e-5, 19.218, 0.0473, 0.773, 74.0, 96.9, 142.2),
    ("neptune", 5.151389e-5, 30.110, 0.0086, 1.770, 131.8, 273.2, 256.2),
]


def new_sim():
    sim = rebound.Simulation()
    sim.G = G
    return sim


def add_giants(sim, names=("jupiter", "saturn", "uranus", "neptune")):
    for name, m, a, e, inc, Om, om, M in GIANTS:
        if name in names:
            sim.add(m=m, a=a, e=e, inc=math.radians(inc), Omega=math.radians(Om), omega=math.radians(om), M=math.radians(M), primary=sim.particles[0])


def export(sim, names):
    sim.move_to_com()
    return [(n, p.m, p.x, p.y, p.z, p.vx, p.vy, p.vz) for n, p in zip(names, sim.particles)]


def giants():
    sim = new_sim()
    sim.add(m=1.0)
    add_giants(sim)
    return export(sim, ["sun"] + [g[0] for g in GIANTS])


def giants_with_kbos(n, seed=1):
    rng = np.random.default_rng(seed)
    sim = new_sim()
    sim.add(m=1.0)
    add_giants(sim)
    names = ["sun"] + [g[0] for g in GIANTS]
    for k in range(n):
        sim.add(m=0.0, a=rng.uniform(38, 50), e=rng.uniform(0, 0.2), inc=math.radians(rng.uniform(0, 20)),
                Omega=rng.uniform(0, 2 * math.pi), omega=rng.uniform(0, 2 * math.pi), M=rng.uniform(0, 2 * math.pi), primary=sim.particles[0])
        names.append(f"kbo{k}")
    return export(sim, names)


def jupiter_encounter():
    sim = new_sim()
    sim.add(m=1.0)
    vj = math.sqrt(G / 5.2)
    sim.add(m=9.547919e-4, x=5.2, vy=vj)
    sim.add(m=0.0, x=5.2 - 0.3, y=-0.25, z=0.01, vy=vj + 1.2e-3)
    return export(sim, ["sun", "jupiter", "tp"])


def sungrazer():
    sim = new_sim()
    sim.add(m=1.0)
    vj = math.sqrt(G / 5.2)
    sim.add(m=9.547919e-4, y=5.2, vx=-vj)
    q, e = 0.05, 0.98
    ra = q / (1 - e) * (1 + e)
    sim.add(m=0.0, x=-ra, vy=-math.sqrt(G * (1 - e) / ra))
    return export(sim, ["sun", "jupiter", "comet"])


def centaurs(n, seed=2):
    """Test particles crossing Uranus' and Neptune's orbits: frequent close encounters."""
    rng = np.random.default_rng(seed)
    sim = new_sim()
    sim.add(m=1.0)
    add_giants(sim)
    names = ["sun"] + [g[0] for g in GIANTS]
    for k in range(n):
        sim.add(m=0.0, a=rng.uniform(20, 30), e=rng.uniform(0.1, 0.3), inc=math.radians(rng.uniform(0, 10)),
                Omega=rng.uniform(0, 2 * math.pi), omega=rng.uniform(0, 2 * math.pi), M=rng.uniform(0, 2 * math.pi), primary=sim.particles[0])
        names.append(f"cen{k}")
    return export(sim, names)


def run_spacerocks(integrator, parts, dt, n_steps, n_out):
    text = f"{integrator} {dt!r} {n_steps} {n_out}\n" + "\n".join(" ".join([p[0]] + [repr(v) for v in p[1:]]) for p in parts)
    out = subprocess.run([str(BIN)], input=text, capture_output=True, text=True, check=True).stdout
    states, energy, seconds = {}, [], None
    for line in out.splitlines():
        f = line.split()
        if f[0] == "state":
            states[f[1]] = np.array([float(x) for x in f[2:]])
        elif f[0] == "energy":
            energy.append((float(f[1]), float(f[2])))
        elif f[0] == "seconds":
            seconds = float(f[1])
    return {"states": states, "energy": energy, "seconds": seconds}


def run_rebound(integrator, parts, dt, n_steps, n_out, n_active=False, exact_finish_time=1, **opts):
    sim = new_sim()
    for n, m, x, y, z, vx, vy, vz in parts:
        sim.add(m=m, x=x, y=y, z=z, vx=vx, vy=vy, vz=vz)
    if n_active:
        # Massive bodies come first; skips test particle-test particle gravity (spacerocks
        # skips it too) but also stops TRACE checking test particles for pericenter passages.
        sim.N_active = sum(1 for p in parts if p[1] > 0)
    sim.integrator = integrator
    for k, v in opts.items():
        setattr(sim.integrator, k, v)
    sim.dt = dt
    energy = [(0.0, sim.energy())]
    per_out = max(n_steps // max(n_out, 1), 1)
    taken, seconds = 0, 0.0
    while taken < n_steps:
        chunk = min(per_out, n_steps - taken)
        start = time.perf_counter()
        if integrator == "ias15":
            sim.integrate((taken + chunk) * dt, exact_finish_time=exact_finish_time)
        else:
            sim.steps(chunk)
        seconds += time.perf_counter() - start
        taken += chunk
        if n_out:
            energy.append((sim.t, sim.energy()))
    states = {n: np.array([p.x, p.y, p.z, p.vx, p.vy, p.vz]) for (n, *_), p in zip(parts, sim.particles)}
    return {"states": states, "energy": energy, "seconds": seconds}


def timed(fn, *args, **kw):
    """Run `fn` REPEATS times (no energy outputs) and keep the fastest wall time."""
    best = None
    for _ in range(REPEATS):
        r = fn(*args, **kw)
        if best is None or r["seconds"] < best["seconds"]:
            best = r
    return best


def pos_err(a, b, names):
    return max(float(np.linalg.norm(a["states"][n][:3] - b["states"][n][:3])) for n in names)


def max_de(r):
    e0 = r["energy"][0][1]
    return max(abs((e - e0) / e0) for _, e in r["energy"])


def scenario(label, parts, dt, n_steps, compare, n_out=200):
    """Integrate with all six integrators; errors of the `compare` particles against REBOUND IAS15."""
    res = {"label": label, "n_particles": len(parts), "dt_days": dt, "n_steps": n_steps, "span_days": dt * n_steps}
    runs = {
        "spacerocks_wh": run_spacerocks("whfast", parts, dt, n_steps, n_out),
        "spacerocks_trace": run_spacerocks("trace", parts, dt, n_steps, n_out),
        "spacerocks_ias15": run_spacerocks("ias15", parts, dt, n_steps, n_out),
        "spacerocks_ias15_interp": run_spacerocks("ias15i", parts, dt, n_steps, n_out),
        "rebound_whfast_dh": run_rebound("whfast", parts, dt, n_steps, n_out, coordinates="democraticheliocentric"),
        "rebound_trace": run_rebound("trace", parts, dt, n_steps, n_out),
        "rebound_ias15": run_rebound("ias15", parts, dt, n_steps, n_out),
    }
    ref = runs["rebound_ias15"]
    res["pos_err_vs_rebound_ias15_au"] = {k: pos_err(r, ref, compare) for k, r in runs.items()}
    res["max_rel_energy_error"] = {k: max_de(r) for k, r in runs.items()}
    res["wh_vs_whfast_dh_au"] = pos_err(runs["spacerocks_wh"], runs["rebound_whfast_dh"], [p[0] for p in parts])
    res["trace_vs_rebound_trace_au"] = pos_err(runs["spacerocks_trace"], runs["rebound_trace"], compare)
    res["ias15_vs_ias15_au"] = pos_err(runs["spacerocks_ias15"], runs["rebound_ias15"], compare)
    res["ias15_interp_vs_ias15_au"] = pos_err(runs["spacerocks_ias15_interp"], runs["rebound_ias15"], compare)
    return res


def ias15_timing(label, parts, span, outputs, compare):
    """IAS15 wall time and accuracy against output cadence. spacerocks lands on each output with
    `integrate` or interpolates with `integrate_or_interpolate`; REBOUND lands on each output
    (exact_finish_time=1) or overshoots it (exact_finish_time=0, no interpolation)."""
    rows = []
    ref = run_rebound("ias15", parts, span, 1, 1)
    for n in outputs:
        dt = span / n
        runs = {
            "spacerocks_integrate": timed(run_spacerocks, "ias15", parts, dt, n, n),
            "spacerocks_integrate_or_interpolate": timed(run_spacerocks, "ias15i", parts, dt, n, n),
            "rebound_exact_finish": timed(run_rebound, "ias15", parts, dt, n, n, n_active=True),
            "rebound_overshoot": timed(run_rebound, "ias15", parts, dt, n, n, n_active=True, exact_finish_time=0),
        }
        rows.append({"n_outputs": n,
                     "ms": {k: 1e3 * r["seconds"] for k, r in runs.items()},
                     # An overshooting run doesn't end at `span`, so it has no comparable final state.
                     "final_pos_err_vs_rebound_one_output_au": {k: pos_err(r, ref, compare) for k, r in runs.items() if k != "rebound_overshoot"}})
    return {"label": label, "span_days": span, "rows": rows}


def timing(label, parts, dt, n_steps):
    t = {
        "spacerocks_wh": timed(run_spacerocks, "whfast", parts, dt, n_steps, 1)["seconds"],
        "spacerocks_trace": timed(run_spacerocks, "trace", parts, dt, n_steps, 1)["seconds"],
        "rebound_whfast_dh": timed(run_rebound, "whfast", parts, dt, n_steps, 1, n_active=True, coordinates="democraticheliocentric")["seconds"],
        "rebound_whfast_jacobi": timed(run_rebound, "whfast", parts, dt, n_steps, 1, n_active=True)["seconds"],
        "rebound_trace": timed(run_rebound, "trace", parts, dt, n_steps, 1, n_active=True)["seconds"],
        "rebound_trace_all_checked": timed(run_rebound, "trace", parts, dt, n_steps, 1)["seconds"],
    }
    return {"label": label, "n_particles": len(parts), "n_steps": n_steps,
            "us_per_step": {k: 1e6 * v / n_steps for k, v in t.items()}}


def main():
    giant_names = [g[0] for g in GIANTS]
    results = {"rebound_version": rebound.__version__, "accuracy": [], "timing": []}
    acc = results["accuracy"]
    acc.append(scenario("Sun + 4 giants, 1000 yr", giants(), 120.0, 3044, giant_names))
    acc.append(scenario("Sun + 4 giants, 100 kyr (energy)", giants(), 120.0, 304400, giant_names, n_out=400))
    acc.append(scenario("Test particle inside Jupiter's Hill sphere", jupiter_encounter(), 20.0, 150, ["tp"]))
    acc.append(scenario("Sungrazing comet (q = 0.05 AU), 5 passages", sungrazer(), 5.0, 1600, ["comet"]))
    acc.append(scenario("Giants + 100 KBOs, 10 kyr", giants_with_kbos(100), 120.0, 30440, [f"kbo{k}" for k in range(100)]))
    acc.append(scenario("Giants + 50 Centaurs, 1 kyr", centaurs(50), 30.0, 12176, [f"cen{k}" for k in range(50)]))
    for r in acc:
        print(r["label"], file=sys.stderr)
    tim = results["timing"]
    tim.append(timing("Sun + 4 giants", giants(), 120.0, 100000))
    for n in (10, 100, 1000):
        tim.append(timing(f"Giants + {n} KBOs", giants_with_kbos(n), 120.0, max(100000 // n, 200)))
    tim.append(timing("Giants + 50 Centaurs", centaurs(50), 30.0, 5000))
    results["ias15_timing"] = [
        ias15_timing("Sungrazing comet, 8000 days", sungrazer(), 8000.0, (1, 160, 1600, 16000), ["comet"]),
        ias15_timing("Giants + 100 KBOs, 1 kyr", giants_with_kbos(100), 365250.0, (1, 100, 1000, 10000), [f"kbo{k}" for k in range(100)]),
    ]
    for r in results["ias15_timing"]:
        print(r["label"], file=sys.stderr)
    json.dump(results, sys.stdout, indent=1)


if __name__ == "__main__":
    main()
