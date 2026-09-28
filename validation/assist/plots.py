"""Figures for `examples.py`: `python validation/assist/plots.py [outdir]` reads results.json there
and writes PNGs next to it."""
import json
import os
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(__file__), "out")
R = json.load(open(os.path.join(OUT, "results.json")))

# Categorical slots in fixed order (validated: adjacent CVD dE >= 9.1); markers and line styles
# carry identity as well, since two of the hues sit below 3:1 contrast on white.
STYLE = {
    "ASSIST": dict(color="#2a78d6", marker="o", ls="-"),
    "spacerocks": dict(color="#eb6834", marker="s", ls="--"),
    "spacerocks (global)": dict(color="#1baf7a", marker="^", ls="-."),
    "spacerocks batch": dict(color="#eda100", marker="D", ls=":"),
}
LABEL = {
    "ASSIST": "ASSIST (global step rule)",
    "spacerocks": "spacerocks (PRS23 step rule, default)",
    "spacerocks (global)": "spacerocks (global step rule)",
    "spacerocks batch": "spacerocks RockCollection.propagate",
}
CODES = ["ASSIST", "spacerocks", "spacerocks (global)"]

plt.rcParams.update({
    "figure.dpi": 150, "font.size": 9, "axes.titlesize": 10, "axes.spines.top": False,
    "axes.spines.right": False, "axes.grid": True, "grid.color": "#e4e3df", "grid.linewidth": 0.6,
    "axes.edgecolor": "#8a8984", "axes.labelcolor": "#2b2a27", "xtick.color": "#52514e",
    "ytick.color": "#52514e", "legend.frameon": False, "lines.linewidth": 1.6, "lines.markersize": 4,
})


def line(ax, x, y, code, every=None, **kw):
    st = dict(STYLE[code])
    if every:
        st["markevery"] = every
    ax.plot(x, y, label=LABEL[code], **st, **kw)


def save(fig, name):
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, name), facecolor="white")
    plt.close(fig)
    print("wrote", name)


# 1. Against JPL's small-body integrator, 10 to 100,000 days.
fig, axs = plt.subplots(1, 3, figsize=(10, 3.2), sharey=True)
for ax, (name, e) in zip(axs, R["sb_long"].items()):
    for c in CODES:
        line(ax, e["span_days"], e[c]["err_m"], c)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_title(name); ax.set_xlabel("integration time (days)")
axs[0].set_ylabel("position error vs JPL (m)")
axs[0].legend(loc="upper left", fontsize=7)
save(fig, "fig_jpl_long.png")

# 2. Apophis daily positions against JPL's small-body integrator.
e = R["apophis_daily"]
fig, ax = plt.subplots(figsize=(6.5, 3.2))
t = np.array(e["series"]["t"])
for c in CODES:
    line(ax, t, e["series"][c], c, every=50)
ax.set_yscale("log"); ax.set_xlabel("days from 2029-01-01"); ax.set_ylabel("position error vs JPL (m)")
ax.axvline(103.5, color="#8a8984", lw=0.8, ls=":")
ax.text(106, ax.get_ylim()[0] * 3, "Earth flyby\n2029-04-13", fontsize=7, color="#52514e")
ax.set_title("Apophis, daily positions vs JPL small-body integration")
ax.legend(fontsize=7)
save(fig, "fig_apophis_daily.png")

# 3. Apophis 2029 encounter: step sizes, and difference from ASSIST.
fig, axs = plt.subplots(1, 2, figsize=(10, 3.2))
for c, (ts, dts) in R["apophis_steps"].items():
    line(axs[0], ts, np.abs(dts), c, every=max(1, len(ts) // 25))
axs[0].set_yscale("log"); axs[0].set_xlabel("days from 2029-01-01"); axs[0].set_ylabel("step size (days)")
axs[0].set_title("Step size through the Apophis flyby (min_dt = 0.001 d)")
axs[0].legend(fontsize=7)
s = R["apophis_2029"]["series"]
for c, d in s["diff_from_assist_m"].items():
    line(axs[1], s["t"], np.maximum(d, 1e-6), c, every=500)
axs[1].set_yscale("log"); axs[1].set_xlabel("days from 2029-01-01"); axs[1].set_ylabel("difference from ASSIST (m)")
axs[1].set_title("Position difference from ASSIST")
save(fig, "fig_apophis_2029.png")

# 4. Round trip.
e = R["round_trip"]
fig, ax = plt.subplots(figsize=(5, 3.2))
for c in CODES:
    line(ax, e["span_days"], np.maximum(e[c]["err_m"], 1e-6), c)
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("outgoing time (days)"); ax.set_ylabel("round-trip error (m)")
ax.set_title("Holman forward and back")
ax.legend(fontsize=7)
save(fig, "fig_round_trip.png")

# 5. Variational equations vs shadow particle.
e = R["variational"]
fig, axs = plt.subplots(2, 1, figsize=(6.5, 4.2), sharex=True)
for c in CODES:
    line(axs[0], e["t"], e[c]["shadow_km"], c, every=40)
    line(axs[1], e["t"], np.maximum(e[c]["residual_mm"], 1e-6), c, every=40)
axs[0].set_ylabel("x displacement (km)")
axs[0].set_title("Holman: shadow particle displaced 1e-8 AU in x")
axs[1].set_yscale("log"); axs[1].set_ylabel("|shadow - variational| (mm)"); axs[1].set_xlabel("days")
axs[0].legend(fontsize=7)
save(fig, "fig_variational.png")

# 6. (5303)-Ceres encounter and convergence studies.
e = R["ceres_5303"]
fig, axs = plt.subplots(1, 3, figsize=(11, 3.2))
for c in CODES:
    line(axs[0], e["t"], e[c]["separation_au"], c, every=300)
axs[0].set_yscale("log"); axs[0].set_xlabel("days from 1995-01-01"); axs[0].set_ylabel("separation from Ceres (AU)")
axs[0].set_title("(5303) Parijskij - Ceres")
axs[0].legend(fontsize=7)
for ax, (name, title) in zip(axs[1:], [("ceres_5303", "Ceres encounter: final position"),
                                       ("apophis_2029", "Apophis flyby: final position")]):
    cv = R["convergence"][name]
    for c in CODES:
        eps = [x for x, y in zip(cv["eps"], cv[c]) if y is not None]
        err = [max(y, 1e-4) for y in cv[c] if y is not None]
        line(ax, eps, err, c)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.invert_xaxis()
    ax.set_xlabel("IAS15 epsilon"); ax.set_ylabel("distance from converged answer (m)")
    ax.set_title(title + " vs epsilon")
save(fig, "fig_ceres_convergence.png")

# 7. Speed, as a speed-up over ASSIST (ASSIST's time / spacerocks' time) on a log2 axis, so each
#    gridline is a factor of 2.
from matplotlib.ticker import FixedLocator, NullLocator


def log2_axis(ax, lo, hi):
    ax.set_yscale("log", base=2)
    ticks = [2.0 ** k for k in range(int(np.floor(np.log2(lo))), int(np.ceil(np.log2(hi))) + 1)]
    ax.yaxis.set_major_locator(FixedLocator(ticks))
    ax.yaxis.set_minor_locator(NullLocator())
    ax.set_yticklabels([f"{t:g}×" for t in ticks])
    ax.set_ylim(ticks[0], ticks[-1])
    ax.axhline(1.0, color="#2a78d6", lw=1.2)


names = [("holman_30d", "Holman\n30 d"), ("getting_started", "Getting\nstarted"), ("apophis_2029", "Apophis\n2029"),
         ("apophis_daily", "Apophis\ndaily"), ("variational", "Variational\n10,000 d"), ("ceres_5303", "5303-Ceres\n10 yr"),
         ("round_trip", "Round trip\n(21 runs)")]
SR = ["spacerocks", "spacerocks (global)"]
fig, axs = plt.subplots(1, 2, figsize=(11, 3.4))
w = 0.36
x = np.arange(len(names))
allv = []
for i, c in enumerate(SR):
    vals = [R[k]["ASSIST"]["time_s"] / R[k][c]["time_s"] for k, _ in names]
    allv += vals
    bars = axs[0].bar(x + (i - 0.5) * w, vals, w - 0.03, color=STYLE[c]["color"], label=LABEL[c])
    for b_, v in zip(bars, vals):
        axs[0].text(b_.get_x() + b_.get_width() / 2, v * 1.06, f"{v:.1f}", ha="center", va="bottom", fontsize=6.5, color="#52514e")
log2_axis(axs[0], min(1.0, min(allv)), max(allv) * 1.3)
axs[0].set_xticks(x, [n for _, n in names], fontsize=7)
axs[0].set_ylabel("speed-up over ASSIST"); axs[0].set_title("Each example, called from Python (ASSIST = 1×)")
axs[0].legend(fontsize=7, loc="upper left")
b = R["benchmark"]
allv = []
for c in ["spacerocks", "spacerocks (global)", "spacerocks batch"]:
    v = np.array(b["ASSIST"]) / np.array(b[c])
    allv += list(v)
    line(axs[1], b["n"], v, c)
log2_axis(axs[1], min(1.0, min(allv)), max(allv) * 1.3)
axs[1].set_xscale("log")
axs[1].set_xlabel("particles"); axs[1].set_ylabel("speed-up over ASSIST")
axs[1].set_title("ASSIST benchmark: N particles, 10 years (ASSIST = 1×)")
axs[1].legend(fontsize=7, loc="upper left")
save(fig, "fig_speed.png")

# 8. Compensated summation: round-trip error of 64 Holman clones, and cost per step.
if "summation" in R:
    e = R["summation"]
    SUM_STYLE = {"none": dict(color="#8a8984", marker="o", ls=":"),
                 "kahan": dict(color="#eb6834", marker="s", ls="--"),
                 "full": dict(color="#2a78d6", marker="^", ls="-")}
    SUM_LABEL = {"none": "plain sums", "kahan": "Kahan on x, v (default)", "full": "Kahan on x, v and corrector (REBOUND)"}
    cfgs = list(e["configs"].items())
    fig, axs = plt.subplots(1, len(cfgs) + 1, figsize=(3.2 * (len(cfgs) + 1), 3.2))
    for ax, (name, cfg) in zip(axs, cfgs):
        for m in e["modes"]:
            ax.plot(e["spans"], np.maximum(cfg[m]["median_m"], 1e-6), label=SUM_LABEL[m], **SUM_STYLE[m])
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel("outgoing time (days)"); ax.set_title(name)
    axs[0].set_ylabel("median round-trip error (m)")
    axs[0].legend(fontsize=7, loc="upper left")
    ax = axs[-1]
    base = e["time_s"]["none"]
    vals = [e["time_s"][m] / base for m in e["modes"]]
    ax.bar(range(len(vals)), vals, color=[SUM_STYLE[m]["color"] for m in e["modes"]])
    ax.set_xticks(range(len(vals)), ["plain", "Kahan x, v", "full"], fontsize=8)
    ax.set_ylim(0.9, max(vals) * 1.03); ax.set_ylabel("time relative to plain sums")
    ax.set_title("1000 clones, 10 years")
    save(fig, "fig_summation.png")

# 9. PRS23 in both codes.
if "prs23_crosscheck" in R:
    e = R["prs23_crosscheck"]
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.2))
    for ax, name, title, ylab in [(axs[0], "ceres_5303", "(5303)-Ceres, PRS23 step rule", "distance from converged answer (m)"),
                                  (axs[1], "apophis_2029", "Apophis flyby, PRS23 step rule", "distance from JPL (m)")]:
        for c in ["ASSIST", "spacerocks"]:
            st = dict(STYLE[c]); lab = c + " (PRS23)"
            ax.plot(e["eps"], e[name][c], label=lab, **st)
        ax.plot(e["eps"], np.maximum(e[name]["ASSIST vs spacerocks"], 1e-4), label="ASSIST vs spacerocks", color="#1baf7a", marker="^", ls="-.")
        ax.set_xscale("log"); ax.set_yscale("log"); ax.invert_xaxis()
        ax.set_xlabel("IAS15 epsilon"); ax.set_ylabel(ylab); ax.set_title(title)
    axs[0].legend(fontsize=7)
    save(fig, "fig_prs23.png")
