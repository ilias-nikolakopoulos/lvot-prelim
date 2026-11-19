"""STEP 7g — Figure: echo LVOT area and VTI err in opposite directions, but do not cancel case by case.

Echo stroke volume = LVOT area (circle from one diameter) x VTI (PW Doppler in the LVOT core).
  area error = echo circle area / true lumen area   (< 1: the elliptical LVOT is underestimated)
  VTI error  = core (gated) velocity / true mean velocity  (> 1: the core is faster than the cross-section mean)
  SV error   = area error x VTI error
On a log scale the two errors add, so the bars in panel A sum to the SV error.
  A  per case: area error (left), VTI error (right) and the net SV error (diamond), sorted by SV error
  B  area error vs VTI error with lines of equal SV error; the errors cancel only on the SV = 1 line
Data: CT-based CFD twins (ImageCAS, real LVOT lumen with the LV attached, steady laminar flow, h = 1 mm),
  imagecas_cfd/cfd/cfd_results_lv_lam_h1.csv (rows twin == "real"); depth set by DEPTH below.
Outputs: outputs/figures/analysis/fig_opposing_errors.png
Usage: python scripts/29_opposing_errors_figure.py [--depth 5.0]
"""
import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import font_manager

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src.provenance import log_run  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--depth", type=float, default=5.0, help="mm below the LV-aorta label boundary (5, 7.5 or 10)")
a = ap.parse_args()

for f_ in sorted((ROOT / "assets/fonts").glob("Roboto-*.ttf")):
    font_manager.fontManager.addfont(str(f_))
if any("Roboto" in f.name for f in font_manager.fontManager.ttflist):
    plt.rcParams["font.family"] = "Roboto"
FG, MUTED, GRID = "#102D40", "#5A6B77", "#D9E0E5"
C_AREA, C_VTI, C_SV = "#009E73", "#CC79A7", FG  # Okabe-Ito green / pink (not the sex colours)

SRC = ROOT / "imagecas_cfd/cfd/cfd_results_lv_lam_h1.csv"
d = pd.read_csv(SRC)
d = d[(d.twin == "real") & (d.depth_mm == a.depth)].copy()
d["area"], d["vti"], d["sv"] = d.area_ratio, d.profile_ratio, d.SV_ratio
assert np.allclose(d.area * d.vti, d.sv, rtol=1e-6)
d = d.sort_values("sv").reset_index(drop=True)
n = len(d)

# STEP 7g.1 — panel A: log-scale diverging bars, ticks labelled in %
fig, (axA, axB) = plt.subplots(1, 2, figsize=(13, 5.6), gridspec_kw=dict(width_ratios=[1.15, 1]))
y = np.arange(n)
la, lv, ls = np.log(d.area), np.log(d.vti), np.log(d.sv)
axA.barh(y + .17, la, height=.32, color=C_AREA, label="LVOT area (circle from one diameter)")
axA.barh(y - .17, lv, height=.32, color=C_VTI, label="VTI (velocity in the LVOT core)")
axA.scatter(ls, y, marker="D", s=46, color=C_SV, zorder=3, label="Stroke volume = area x VTI")
axA.axvline(0, color=FG, lw=1)
pct = np.array([-40, -30, -20, -10, 0, 10, 20, 30, 40])
axA.set_xticks(np.log(1 + pct / 100), [f"{p:+d}%" if p else "0" for p in pct])
lim = np.log([0.6, 1.85])
axA.set_xlim(*lim)
axA.set_yticks(y, [f"Case {i} ({arm[:3]}.)" for i, arm in zip(d.id, d.arm)], fontsize=8.5)
for yi, s in zip(y, d.sv):
    axA.text(lim[1] - .01, yi, f"SV {100 * (s - 1):+.0f}%", ha="right", va="center", fontsize=8.5,
             color=FG, fontweight="bold" if abs(s - 1) > .1 else "normal")
axA.set_xlabel("Echo error (echo / true)", color=FG)
axA.set_title("A   Every case: the two errors pull in opposite directions,\n      but by different amounts",
              loc="left", fontsize=11.5, color=FG)
axA.legend(loc="upper center", bbox_to_anchor=(.5, -.13), ncol=3, fontsize=8.5, frameon=False, handlelength=1.4,
           columnspacing=1.2)
axA.grid(axis="x", color=GRID, lw=.6)
axA.set_axisbelow(True)

# STEP 7g.2 — panel B: area error vs VTI error, lines of equal SV error (VTI = k / area)
xs = np.linspace(.55, 1.05, 200)
for k in (.6, .7, .8, .9, 1.0, 1.1, 1.2):
    yy = k / xs
    ok = (yy > .95) & (yy < 1.5)
    axB.plot(xs[ok], yy[ok], color=FG if k == 1 else GRID, lw=1.6 if k == 1 else 1, zorder=1)
    i = np.where(ok)[0][-1]  # label where the line leaves the plot (bottom or right edge), away from the points
    lab = "SV exact" if k == 1 else f"SV {100 * (k - 1):+.0f}%"
    if k == 1:  # the exact line is labelled at its top, clear of the axes
        axB.text(.705, 1.49, lab, fontsize=8.5, color=FG, fontweight="bold", va="top")
    elif yy[i] < .97:
        axB.text(xs[i] - .006, .955, lab, fontsize=8, color=FG if k == 1 else MUTED, ha="right", va="bottom",
                 fontweight="bold" if k == 1 else "normal")
    else:
        axB.text(1.045, yy[i] + .008, lab, fontsize=8, color=MUTED, ha="right", va="bottom")
axB.scatter(d.area, d.vti, s=48, color=C_SV, zorder=3)
OFF = {978: (-6, 7), 663: (-16, -11), 167: (5, 1), 467: (-18, 4)}  # untangle three overlapping points
for _, r in d.iterrows():
    axB.annotate(str(r.id), (r.area, r.vti), xytext=OFF.get(r.id, (4, 3)), textcoords="offset points",
                 fontsize=7.5, color=MUTED)
gx, gy = np.exp(la.mean()), np.exp(lv.mean())
axB.scatter([gx], [gy], marker="x", s=90, color=FG, lw=2, zorder=4)
axB.annotate(f"average of {n}:\narea {100 * (gx - 1):+.0f}%, VTI {100 * (gy - 1):+.0f}%\nSV {100 * (gx * gy - 1):+.0f}%",
             (gx, gy), xytext=(.565, 1.36), textcoords="data", fontsize=8.5, color=FG,
             arrowprops=dict(arrowstyle="-", color=FG, lw=.7))
axB.axvline(1, color=MUTED, lw=.7, ls=":")
axB.axhline(1, color=MUTED, lw=.7, ls=":")
axB.set_xlim(.55, 1.05)
axB.set_ylim(.95, 1.5)
axB.set_xticks([.6, .7, .8, .9, 1.0], ["-40%", "-30%", "-20%", "-10%", "0"])
axB.set_yticks([1.0, 1.1, 1.2, 1.3, 1.4, 1.5], ["0", "+10%", "+20%", "+30%", "+40%", "+50%"])
axB.set_xlabel("LVOT area error (echo circle vs true lumen)", color=FG)
axB.set_ylabel("VTI error (core velocity vs true mean)", color=FG)
axB.set_title("B   They cancel only on the black line;\n      the spread around it is the per-patient SV error",
              loc="left", fontsize=11.5, color=FG)

for ax in (axA, axB):
    ax.tick_params(colors=FG, labelsize=8.5)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(MUTED)
sv_pct = 100 * (d.sv - 1)
fig.text(.01, .01, f"CT-based CFD twins (ImageCAS, n = {n}; real lumen with LV attached, steady laminar flow), sample "
         f"{a.depth:g} mm below the LV-aorta label boundary. SV error range {sv_pct.min():+.0f}% to "
         f"{sv_pct.max():+.0f}%. Log scale: the area and VTI bars add up to the SV diamond.",
         fontsize=8, color=MUTED)
fig.tight_layout(rect=(0, .03, 1, 1))
OUT = ROOT / "outputs/figures/analysis/fig_opposing_errors.png"
fig.savefig(OUT, dpi=200, facecolor="white")
plt.close(fig)
print(d[["id", "pair", "arm", "area", "vti", "sv"]].round(3).to_string(index=False))
print(f"geometric mean: area {gx:.3f}, VTI {gy:.3f}, SV {gx * gy:.3f}; SV range {d.sv.min():.3f}-{d.sv.max():.3f}")
log_run(__file__, depth_mm=a.depth, n_cases=n, input=str(SRC.relative_to(ROOT)),
        outputs=[str(OUT.relative_to(ROOT))])
