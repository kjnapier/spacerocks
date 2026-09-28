"""Cross-check spacerocks' force models and variational equations against ASSIST.

    pip install assist
    cargo build --release --example assist_compare
    SPACEROCKS_KERNELS=/path/to/kernels python validation/assist/compare.py

Integrates five test particles (a close Earth approacher, an Apophis-like orbit with A2, a
near-Sun orbit with A1-A3, a main-belt object and a TNO) for 400 days with each force switched
on in turn, in ASSIST and in spacerocks, and compares the effect of every force, the full
default model, and the state transition matrix and non-gravitational partials.

ASSIST 1.2.3 misreads the Jupiter GM (GM5) from NAIF's de440s.bsp comments because it falls on a
comment-record boundary, and then integrates with Jupiter massless. This script gives ASSIST a
temporary copy of de440s.bsp with that boundary patched.
"""
import os, sys, tempfile
import assist, rebound, numpy as np, subprocess

KDIR = os.environ["SPACEROCKS_KERNELS"]
DRIVER = os.path.join(os.path.dirname(__file__), "..", "..", "target", "release", "examples", "assist_compare")

def patched_de440s():
    data = bytearray(open(os.path.join(KDIR, "de440s.bsp"), "rb").read())
    start = data.find(b"Initial conditions and constants")
    i = data.find(b"GM5 ", start)
    if i > 0 and i % 1024 == 0 and data[i - 1] == 0:
        data[i - 1] = ord("\n")
    f = tempfile.NamedTemporaryFile(suffix=".bsp", delete=False)
    f.write(data)
    f.close()
    return f.name

ephem = assist.Ephem(patched_de440s(), os.path.join(KDIR, "sb441-n16.bsp"))
JREF = ephem.jd_ref
T0 = 2460300.5; T1 = T0 + 400.0
AU_KM = 1.495978707e8
e = ephem.get_particle("earth", T0 - JREF)
sun = ephem.get_particle("sun", T0 - JREF)
# test objects: close Earth approacher, Apophis-like with nongrav, near-Sun object, main-belt, TNO
objs = []
# 1. NEO 0.003 AU from Earth moving slowly relative to it
objs.append([e.x+0.002, e.y-0.001, e.z+0.001, e.vx+0.001, e.vy-0.0005, e.vz+0.0003, 0,0,0])
# 2. Apophis-like elliptical with nongrav A2
objs.append([0.6, 0.55, 0.2, -0.014, 0.012, 0.004, 0, -2.9e-14, 0])
# 3. comet-ish / near-Sun: perihelion ~0.1 AU with A1,A2,A3
objs.append([0.3, 0.1, 0.05, -0.01, 0.028, 0.004, 1e-8, -2e-9, 5e-10])
# 4. main belt
objs.append([2.5, -1.2, 0.3, 0.005, 0.009, 0.001, 0,0,0])
# 5. TNO
objs.append([30.0, 20.0, 5.0, -0.0015, 0.0022, 0.0003, 0,0,0])

def run_assist(obj, forces, eps=1e-9, var=False):
    sim = rebound.Simulation()
    sim.t = T0 - JREF
    sim.add(x=obj[0], y=obj[1], z=obj[2], vx=obj[3], vy=obj[4], vz=obj[5])
    params = [obj[6], obj[7], obj[8]]
    if var:
        for i in range(6):
            vp = sim.add_variation(testparticle=0)
            setattr(vp.particles[0], "xyz"[i] if i < 3 else "v"+"xyz"[i-3], 1.0)
            params += [0,0,0]
        for k in range(3):
            vp = sim.add_variation(testparticle=0)
            dp = [0,0,0]; dp[k] = 1.0
            params += dp
    ex = assist.Extras(sim, ephem)
    ex.forces = forces
    ex.particle_params = np.array(params)
    sim.ri_ias15.epsilon = eps
    sim.integrate(T1 - JREF)
    out = list(sim.particles[0].xyz) + list(sim.particles[0].vxyz)
    if var:
        for i in range(1, sim.N):
            p = sim.particles[i]
            out += list(p.xyz) + list(p.vxyz)
    return np.array(out)

def run_ours(obj, cfg, eps=1e-9, var=False):
    inp = f"{T0} {T1} {cfg} {eps} {'var' if var else 'novar'}\n" + " ".join(repr(x) for x in obj) + "\n"
    out = subprocess.run([DRIVER], input=inp.encode(), capture_output=True, check=True).stdout.decode()
    return np.array([float(x) for x in out.split()])

base_a = ["SUN", "PLANETS", "ASTEROIDS"]
configs = [
    ("newtonian", [], ""),
    ("+ nongrav", ["NON_GRAVITATIONAL"], "ng"),
    ("+ earth J2-J4", ["EARTH_HARMONICS"], "earth"),
    ("+ solar J2", ["SUN_HARMONICS"], "sunj2"),
    ("+ GR EIH (Sun)", ["GR_EIH"], "eih"),
    ("+ GR simple", ["GR_SIMPLE"], "grsimple"),
    ("+ GR potential", ["GR_POTENTIAL"], "grpot"),
]
names = ["close-Earth NEO", "Apophis-like+NG", "near-Sun+NG", "main belt", "TNO"]
print("400-day integrations. For each force: size of its effect (km) and |ASSIST effect - ours| (km); then full-model difference.")
for oi, obj in enumerate(objs):
    a0 = run_assist(obj, base_a); o0 = run_ours(obj, "newton")
    print(f"\n{names[oi]}: newtonian-only |ASSIST - ours| = {np.linalg.norm(a0[:3]-o0[:3])*AU_KM:.4f} km")
    for label, af, of in configs[1:]:
        a = run_assist(obj, base_a + af); o = run_ours(obj, of + ",newton")
        ea = (a - a0)[:3]; eo = (o - o0)[:3]
        print(f"   {label:16s} effect {np.linalg.norm(ea)*AU_KM:12.4f} km   diff {np.linalg.norm(ea-eo)*AU_KM:.2e} km")
    full_a = run_assist(obj, base_a + ["NON_GRAVITATIONAL","EARTH_HARMONICS","SUN_HARMONICS","GR_EIH"])
    full_o = run_ours(obj, "ng,earth,sunj2,eih,newton")
    print(f"   full ASSIST default model: |ASSIST - ours| = {np.linalg.norm(full_a[:3]-full_o[:3])*AU_KM:.4f} km")


print("\nVariational particles (full default model): max relative difference from ASSIST per variation")
full_a = ["SUN", "PLANETS", "ASTEROIDS", "NON_GRAVITATIONAL","EARTH_HARMONICS","SUN_HARMONICS","GR_EIH"]
labels = ["x","y","z","vx","vy","vz","A1","A2","A3"]
for oi, obj in enumerate(objs):
    a = run_assist(obj, full_a, var=True).reshape(-1, 6)
    o = run_ours(obj, "ng,earth,sunj2,eih,newton", var=True).reshape(-1, 6)
    worst = []
    for j in range(1, 10):
        den = np.abs(a[j]).max(); rel = np.abs(a[j] - o[j]).max() / den if den > 0 else np.abs(o[j]).max()
        worst.append(f"{labels[j-1]} {rel:.1e}")
    print(f"{names[oi]:16s}", ", ".join(worst))
