"""STEP 6.6 — Compare the measured outline area with ellipse- and echo-style area estimates.

For the aortic annulus, the strict LVOT section at 2 mm and the chord-closed section at 4 mm:
  true             polygon (shoelace) area of the outline                      [reference]
  ellipse_moment   pi/4 * major * minor   (moment-equivalent ellipse)
  ellipse_feret    pi/4 * Feret max * Feret min   (caliper max/min)
  ellipse_2view    pi/4 * D_plax * D_perp  (two orthogonal echo-style diameters, e.g. biplane)
  echo_circle      pi/4 * D_plax^2         (standard echo: PLAX diameter, circular LVOT assumed)
D_plax = edge-to-edge line through the outline centroid towards the mitral centre (septum -> aorto-mitral
side), i.e. the direction of the parasternal long-axis cut; D_perp = the orthogonal line.
Outputs: outputs/tables/area_methods_summary.csv, area_methods_by_mesh.csv,
         outputs/figures/analysis/fig_area_methods.png
Caveat: meshes are end-diastolic; clinical LVOT diameters are measured in mid-systole.
Usage: python scripts/06_area_methods.py [--in lvot_all]
"""
import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src.config import load_config  # noqa: E402
from src.provenance import log_run  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--in", dest="inp", default="lvot_all")
a = ap.parse_args()
cfg = load_config()
COL = cfg["analysis"]["colours"]
LAB = {"F": "Women", "M": "Men"}
T = ROOT / "outputs/tables"
FIG = ROOT / "outputs/figures/analysis"
d = pd.read_csv(T / f"{a.inp}.csv")
need = ["ao_d_plax", "lvot2_d_plax", "lvot4cc_d_plax"]
if any(c not in d for c in need):
    sys.exit(f"{a.inp}.csv lacks echo diameters ({need}); rerun 01_measure.py with the current code")

LEVELS = [  # name, prefix for area/diameters
    ("Aortic annulus", "ao_"),
    ("LVOT 2 mm (closed)", "lvot2_"),
    ("LVOT 4 mm (chord-closed)", "lvot4cc_"),
]
METHODS = [
    ("ellipse_moment", "Ellipse (moment major × minor)"),
    ("ellipse_feret", "Ellipse (Feret max × min)"),
    ("ellipse_2view", "Ellipse (PLAX × orthogonal)"),
    ("echo_circle", "Circle from PLAX diameter (echo)"),
]
q = np.pi / 4
rows = []
for lev, p in LEVELS:
    g = pd.DataFrame({
        "mesh_id": d.mesh_id, "sex": d.sex, "age": d.age, "bmi": d.bmi, "level": lev,
        "true": d[p + "area"],
        "ellipse_moment": q * d[p + "dmax"] * d[p + "dmin"],
        "ellipse_feret": q * d[p + "feret_max"] * d[p + "feret_min"],
        "ellipse_2view": q * d[p + "d_plax"] * d[p + "d_perp"],
        "echo_circle": q * d[p + "d_plax"] ** 2,
        "d_plax": d[p + "d_plax"], "d_perp": d[p + "d_perp"], "d_area": d[p + "d_area"],
        "plax_vs_minor_deg": d[p + "plax_vs_minor_deg"],
    })
    rows.append(g)
L = pd.concat(rows).dropna(subset=["true"])
for m, _ in METHODS:
    L[m + "_pct"] = 100 * (L[m] - L["true"]) / L["true"]
L.round(3).to_csv(T / "area_methods_by_mesh.csv", index=False)

summ = []
for (lev, sx), g in L.groupby(["level", "sex"], sort=False):
    r = {"level": lev, "sex": sx, "n": len(g), "true_area_mm2": g["true"].mean(),
         "d_area_mm": g["d_area"].mean(), "d_plax_mm": g["d_plax"].mean(), "d_perp_mm": g["d_perp"].mean(),
         "plax_vs_minor_deg_median": g["plax_vs_minor_deg"].median()}
    for m, _ in METHODS:
        r[f"{m}_mm2"] = g[m].mean()
        r[f"{m}_pct_mean"] = g[m + "_pct"].mean()
        r[f"{m}_pct_sd"] = g[m + "_pct"].std()
        r[f"{m}_pct_min"], r[f"{m}_pct_max"] = g[m + "_pct"].min(), g[m + "_pct"].max()
    summ.append(r)
S = pd.DataFrame(summ)
S.round(2).to_csv(T / "area_methods_summary.csv", index=False)

# figure: row 1 = % difference from the true area by method and sex; row 2 = echo circle vs true
fig, axs = plt.subplots(2, 3, figsize=(16, 9))
for j, (lev, _) in enumerate(LEVELS):
    ax = axs[0, j]
    g = L[L.level == lev]
    for i, (m, ml) in enumerate(METHODS):
        for k, sx in enumerate("FM"):
            v = g.loc[g.sex == sx, m + "_pct"].dropna()
            x = i + (k - 0.5) * 0.3
            ax.scatter(x + np.random.default_rng(0).uniform(-.07, .07, len(v)), v, s=8, color=COL[sx], alpha=.45, lw=0,
                       label=LAB[sx] if i == 0 else None)
            ax.plot([x - .1, x + .1], [v.mean()] * 2, color="k", lw=2)
    ax.axhline(0, color="0.4", lw=1)
    ax.set_xticks(range(len(METHODS)))
    ax.set_xticklabels(["Ellipse\nmoment axes", "Ellipse\nFeret max×min", "Ellipse\nPLAX × orthogonal", "Echo circle\nfrom PLAX D"], fontsize=8)
    ax.set_ylabel("% difference from measured outline area")
    ax.set_title(lev, fontsize=10)
    ax.grid(axis="y", alpha=.25)
    if j == 0:
        ax.legend(fontsize=8, frameon=False)
    ax = axs[1, j]
    for sx in "FM":
        h = g[g.sex == sx]
        ax.scatter(h["true"], h["echo_circle"], s=10, color=COL[sx], alpha=.6, lw=0, label=LAB[sx])
    lim = [g[["true", "echo_circle"]].min().min() * 0.9, g[["true", "echo_circle"]].max().max() * 1.05]
    ax.plot(lim, lim, "k--", lw=1, label="identity")
    ax.set_xlim(lim); ax.set_ylim(lim); ax.set_aspect("equal")
    ax.set_xlabel("Measured outline area (mm²)"); ax.set_ylabel("Echo circular area π/4·D_PLAX² (mm²)")
    mb = g["echo_circle_pct"].mean()
    ax.set_title(f"{lev}: echo circle vs measured (mean {mb:+.0f}%)", fontsize=10)
    ax.grid(alpha=.25)
    if j == 0:
        ax.legend(fontsize=8, frameon=False)
fig.suptitle("Area estimates vs the measured outline area (end-diastolic meshes). Black bar = mean.", fontsize=11)
fig.tight_layout(rect=(0, 0, 1, 0.96))
fig.savefig(FIG / "fig_area_methods.png", dpi=120)
plt.close(fig)

pd.set_option("display.width", 250)
show = ["level", "sex", "n", "true_area_mm2", "d_area_mm", "d_plax_mm", "d_perp_mm", "plax_vs_minor_deg_median"] + \
       [f"{m}_pct_mean" for m, _ in METHODS]
print(S[show].round(1).to_string(index=False))
log_run(__file__, input=a.inp, n_meshes=int(d.shape[0]),
        outputs=["outputs/tables/area_methods_summary.csv", "outputs/tables/area_methods_by_mesh.csv",
                 "outputs/figures/analysis/fig_area_methods.png"])
