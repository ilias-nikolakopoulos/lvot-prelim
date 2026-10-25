"""STEP 6.8 — Ellipticity (major / minor diameter) by method, level and sex: median and IQR.

Methods (same diameters as the area estimates):
  moment              moment-ellipse major / minor axis
  feret               widest / narrowest calliper width
  minor-major axis    longest diameter through the centre / the diameter at 90 deg to it
Levels: aortic annulus, LVOT 2 mm (strict, closed), LVOT 4 mm (chord-closed). All meshes (end-diastole).
Outputs: outputs/tables/ellipticity_methods.csv, outputs/figures/analysis/fig_ellipticity_methods.png
Needs lvot_all.csv (STEP 1-5.1) and area_methods_by_mesh.csv (STEP 6.6, for the minor-major axis diameters).
Usage: python scripts/10_ellipticity_methods.py
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import font_manager
from matplotlib.patches import Patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src.provenance import log_run  # noqa: E402

for f_ in sorted((ROOT / "assets/fonts").glob("Roboto-*.ttf")):
    font_manager.fontManager.addfont(str(f_))
plt.rcParams["font.family"] = "Roboto"
NAME_COL = "#102D40"
LEVEL_COL = {"Aortic annulus": "#E6C229", "LVOT 2 mm (closed)": "#009E73", "LVOT 4 mm (chord-closed)": "#7B4FA6"}
LEVEL_TITLE = {"Aortic annulus": "Aortic annulus", "LVOT 2 mm (closed)": "2 mm from annulus",
               "LVOT 4 mm (chord-closed)": "4 mm from annulus"}
PREFIX = {"Aortic annulus": "ao_", "LVOT 2 mm (closed)": "lvot2_", "LVOT 4 mm (chord-closed)": "lvot4cc_"}
METHODS = [("moment", "Moment\nellipse"), ("feret", "Feret\nellipse"), ("minor_major", "Minor–major\naxis ellipse")]

T = ROOT / "outputs/tables"
d = pd.read_csv(T / "lvot_all.csv")
A = pd.read_csv(T / "area_methods_by_mesh.csv")
rows = []
for lev, p in PREFIX.items():
    ct = A[A.level == lev].set_index("mesh_id")[["ct_long", "ct_short"]]
    g = d[["mesh_id", "sex", "age", "bmi"]].assign(
        level=lev,
        moment=d[p + "dmax"] / d[p + "dmin"],
        feret=d[p + "feret_max"] / d[p + "feret_min"],
    ).set_index("mesh_id")
    g["minor_major"] = ct["ct_long"] / ct["ct_short"]
    rows.append(g.reset_index())
E = pd.concat(rows).dropna(subset=["moment"])

summ = []
for lev in PREFIX:
    for m, _ in METHODS:
        for sx in "FM":
            v = E[(E.level == lev) & (E.sex == sx)][m].dropna()
            q1, med, q3 = v.quantile([.25, .5, .75])
            summ.append({"level": lev, "method": m, "sex": sx, "n": len(v), "median": med, "q1": q1, "q3": q3,
                         "iqr": q3 - q1, "p5": v.quantile(.05), "p95": v.quantile(.95)})
S = pd.DataFrame(summ)
S.round(3).to_csv(T / "ellipticity_methods.csv", index=False)
W = S.pivot_table(index=["level", "method"], columns="sex", values=["median", "q1", "q3"], sort=False)
pd.set_option("display.width", 200)
print(S.assign(txt=S["median"].map("{:.2f}".format) + " (" + S.q1.map("{:.2f}".format) + "–" + S.q3.map("{:.2f}".format) + ")")
      .pivot_table(index=["level", "method"], columns="sex", values="txt", aggfunc="first", sort=False).to_string())

# figure: one panel per level; per method a women (lighter) and a men (darker) box; box = IQR, line = median,
# whiskers = 5th-95th percentile, no outlier points
fig, axs = plt.subplots(1, 3, figsize=(16, 6), sharey=True)
for ax, lev in zip(axs, PREFIX):
    col = LEVEL_COL[lev]
    for i, (m, _) in enumerate(METHODS):
        for k, sx in enumerate("FM"):
            v = E[(E.level == lev) & (E.sex == sx)][m].dropna().values
            pos = i + (k - .5) * .36
            bp = ax.boxplot([v], positions=[pos], widths=.3, whis=(5, 95), showfliers=False, patch_artist=True,
                            medianprops=dict(color=NAME_COL, lw=2), whiskerprops=dict(color="0.4"),
                            capprops=dict(color="0.4"), boxprops=dict(edgecolor=col, lw=1.2))
            bp["boxes"][0].set_facecolor(col); bp["boxes"][0].set_alpha({"F": .35, "M": .75}[sx])
            ax.text(pos, np.percentile(v, 95) + .01, f"{np.median(v):.2f}", ha="center", va="bottom", fontsize=9,
                    color=NAME_COL, fontweight="bold")
    ax.axhline(1, color="0.5", lw=.8, ls=":")
    ax.set_xticks(range(len(METHODS)))
    ax.set_xticklabels([f"{i}\n{l}" for i, (_, l) in enumerate(METHODS, 1)], fontsize=9.5, fontweight="bold", color=NAME_COL)
    ax.set_title(LEVEL_TITLE[lev], fontsize=13, fontweight="bold", pad=12)
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)
    ax.legend(handles=[Patch(fc=col, alpha=.35, ec=col, label="Women"), Patch(fc=col, alpha=.75, ec=col, label="Men")],
              fontsize=9, frameon=False, loc="lower right")
    ax.set_ylim(0.98, 1.75)
axs[0].set_ylabel("Ellipticity (major / minor diameter)")
fig.suptitle("LVOT ellipticity by method", fontsize=16, fontweight="bold", y=0.99)
fig.text(0.5, 0.005, "box = interquartile range, line and number = median, whiskers = 5th–95th percentile; all 1,423 meshes (end-diastole)",
         ha="center", fontsize=9, color="0.35")
fig.tight_layout(rect=(0, 0.03, 1, 0.95))
out = ROOT / "outputs/figures/analysis/fig_ellipticity_methods.png"
fig.savefig(out, dpi=125)
print("wrote", out)
log_run(__file__, n_meshes=int(E.mesh_id.nunique()), outputs=["outputs/tables/ellipticity_methods.csv", str(out.relative_to(ROOT))])
