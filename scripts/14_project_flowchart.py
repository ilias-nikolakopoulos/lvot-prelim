"""STEP 7c — Project flowchart, PRISMA style (top-down, counts in the main column, exclusions to the right).

Phases (left bars): Identification, Processing, Eligibility, Measurement, Analysis, Results.
All counts are recomputed from outputs/tables/lvot_all.csv so the chart stays correct after a rerun.
House style: Roboto (assets/fonts), navy #102D40, accent #B94F00, result colours as in the figures.
Output: outputs/figures/project_flowchart.png (portrait)
Usage: python scripts/14_project_flowchart.py
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import font_manager
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src.provenance import log_run  # noqa: E402

for f_ in sorted((ROOT / "assets/fonts").glob("Roboto-*.ttf")):
    font_manager.fontManager.addfont(str(f_))
plt.rcParams["font.family"] = "Roboto"
NAVY, ORANGE, TINT, INK = "#102D40", "#B94F00", "#EEF3F6", "#1E2A33"

# ---- counts from the data ----
d = pd.read_csv(ROOT / "outputs/tables/lvot_all.csv")
err_file = ROOT / "outputs/tables/lvot_all_errors.csv"
n_err = len(pd.read_csv(err_file)) if err_file.exists() else 0
N, nF, nM = len(d), int((d.sex == "F").sum()), int((d.sex == "M").sum())
miss2 = d[d.lvot2_d_area.isna()]
n2, n4 = int(d.lvot2_d_area.notna().sum()), int(d.lvot4cc_d_area.notna().sum())
tube3 = round(100 * (d.lvot_tube_len == 3).mean())
F = d[d.sex == "F"].set_index(["age", "bmi"]); M = d[d.sex == "M"].set_index(["age", "bmi"])
npair = len(F.index.intersection(M.index))
fmt = "{:,}".format

fig = plt.figure(figsize=(10, 13.5))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 10); ax.set_ylim(0, 13.5); ax.axis("off")
ax.text(5, 13.0, "LVOT size and ellipticity in UK Biobank hearts", fontsize=20, fontweight="bold", color=NAVY, ha="center")
ax.text(5, 12.58, "Study flow", fontsize=13, color="0.35", ha="center")

MX0, MX1 = 1.25, 6.35          # main column
EX0, EX1 = 6.85, 9.7           # exclusion column
MC = (MX0 + MX1) / 2
LC, RC = MX0 + (MC - MX0) / 2 - .05, MC + (MX1 - MC) / 2 + .05


def box(x0, y0, x1, y1, lines, ec=NAVY, fc="white", lw=1.6, size=11.5, align="center"):
    ax.add_patch(FancyBboxPatch((x0, y0), x1 - x0, y1 - y0, boxstyle="round,pad=0,rounding_size=0.08",
                                fc=fc, ec=ec, lw=lw))
    n = len(lines)
    step = min(.34, (y1 - y0 - .25) / max(n - 1, 1)) if n > 1 else 0
    top = (y0 + y1) / 2 + step * (n - 1) / 2
    xt = (x0 + x1) / 2 if align == "center" else x0 + .2
    for i, ln in enumerate(lines):
        ax.text(xt, top - i * step, ln, ha=align, va="center", fontsize=size + .5 if i == 0 else size,
                fontweight="bold" if i == 0 else "normal", color=NAVY if i == 0 else INK)


def arrow(x0, y0, x1, y1):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=15, color=NAVY, lw=1.4,
                                 shrinkA=0, shrinkB=0))


def line(x0, y0, x1, y1):
    ax.plot([x0, x1], [y0, y1], color=NAVY, lw=1.4)


def phase(y0, y1, label):
    ax.add_patch(FancyBboxPatch((0.3, y0), 0.6, y1 - y0, boxstyle="round,pad=0,rounding_size=0.08", fc=NAVY, ec="none"))
    ax.text(0.6, (y0 + y1) / 2, label, rotation=90, ha="center", va="center", fontsize=11, fontweight="bold", color="white")


def exclusion(y0, y1, lines, from_y):
    box(EX0, y0, EX1, y1, lines, ec=ORANGE, fc=TINT, lw=1.4, size=10.5, align="left")
    arrow(MC, from_y, EX0, from_y)


def split(y_from, y_to):
    """Main stream splits into two columns (left/right half of the main column)."""
    ys = y_from - .3
    line(MC, y_from, MC, ys); line(LC, ys, RC, ys)
    arrow(LC, ys, LC, y_to); arrow(RC, ys, RC, y_to)
    return ys


Y = {"id": (11.0, 12.2), "proc": (9.45, 10.45), "sec": (7.95, 8.95), "elig": (6.3, 7.3),
     "meas": (4.7, 5.95), "ana": (3.1, 4.2), "res": (0.45, 2.5)}

# Identification
phase(*Y["id"], "Identification")
box(MX0, Y["id"][0], MX1, Y["id"][1], ["Representative meshes, UK Biobank", "Ugurlu et al., PLoS One 2025",
                                       f"n = {fmt(N)}  ({fmt(nF)} women, {fmt(nM)} men)",
                                       "averaged by sex, age (46–83) and BMI; end-diastole"])
arrow(MC, Y["id"][0], MC, Y["proc"][1])

# Processing
phase(Y["sec"][0], Y["proc"][1], "Processing")
box(MX0, Y["proc"][0], MX1, Y["proc"][1], ["Mesh geometry", "LV cavity, aortic and mitral rims, axes",
                                           f"n = {fmt(N - n_err)}  (3 closed surfaces in every mesh)"])
gap_y = (Y["proc"][0] + Y["sec"][1]) / 2  # connector between geometry and sections
exclusion(gap_y - .32, gap_y + .32, ["Processing errors", f"n = {n_err}"], gap_y)
arrow(MC, Y["proc"][0], MC, Y["sec"][1])
box(MX0, Y["sec"][0], MX1, Y["sec"][1], ["LVOT sections", "planes perpendicular to the aortic axis, every 1 mm",
                                         f"closed tube 2–5 mm long ({tube3}% exactly 3 mm)"])

# Eligibility: 2 mm and 4 mm
phase(*Y["elig"], "Eligibility")
ys = split(Y["sec"][0], Y["elig"][1])
box(MX0, Y["elig"][0], MC - .1, Y["elig"][1], ["Primary: 2 mm", "strict closed section", f"n = {fmt(n2)}"], size=10.5)
box(MC + .1, Y["elig"][0], MX1, Y["elig"][1], ["Secondary: 4 mm", "gap bridged by a chord", f"n = {fmt(n4)}"], size=10.5)
exclusion(Y["elig"][1] - .55, Y["elig"][1] + .5,
          ["Excluded from 2 mm", f"n = {len(miss2)}  ({int((miss2.sex == 'F').sum())} women, "
           f"{int((miss2.sex == 'M').sum())} men)", "closed tube only to 1 mm"], ys)

# Measurement
phase(*Y["meas"], "Measurement")
arrow(LC, Y["elig"][0], LC, Y["meas"][1]); arrow(RC, Y["elig"][0], RC, Y["meas"][1])
box(MX0, Y["meas"][0], MX1, Y["meas"][1], ["Measured per mesh", "area (shoelace), diameters,",
                                           "ellipticity (3 ellipse methods), echo PLAX circle,", "LV volume and length"], size=10.5)

# Analysis: matched pairs only (regression removed, DECISIONS #42)
phase(*Y["ana"], "Analysis")
arrow(MC, Y["meas"][0], MC, Y["ana"][1])
gap2 = (Y["meas"][0] + Y["ana"][1]) / 2
box(MX0, Y["ana"][0], MX1, Y["ana"][1], ["Men vs women, age and BMI held fixed",
                                         f"{npair} same-age, same-BMI pairs ({fmt(2 * npair)} meshes)",
                                         "paired difference, overall and by age group"], size=10.5)
exclusion(gap2 - .55, gap2 + .55,
          ["No same-age, same-BMI partner", f"n = {len(F) + len(M) - 2 * npair}  ({len(F) - npair} women, {len(M) - npair} men)",
           "not in the sex comparison"], gap2)

# Results
phase(*Y["res"], "Results")
arrow(MC, Y["ana"][0], MC, Y["res"][1])
ax.add_patch(FancyBboxPatch((MX0, Y["res"][0]), EX1 - MX0, Y["res"][1] - Y["res"][0],
                            boxstyle="round,pad=0,rounding_size=0.08", fc="white", ec=NAVY, lw=1.6))
RES = [("~3 mm", "men's annulus and LVOT larger at every age (15–16%);\n~6% after indexing to LV size", "#0072B2"),
       ("1.14 vs 1.4–1.5", "round annulus, oval LVOT just below it", "#009E73"),
       ("+0.07 to +0.10", "women's LVOT more oval than men's, at every age", "#7B4FA6"),
       ("−30%", "echo circular area underestimates LVOT area\n(women −33%, men −29%)", ORANGE)]
for i, (big, txt, col) in enumerate(RES):
    y = Y["res"][1] - .3 - i * .5
    ax.text(MX0 + .25, y, big, fontsize=14, fontweight="bold", color=col, va="center")
    ax.text(MX0 + 2.3, y, txt, fontsize=10.5, color=INK, va="center", linespacing=1.15)

out = ROOT / "outputs/figures/project_flowchart.png"
fig.savefig(out, dpi=150)
print("wrote", out)
log_run(__file__, output=str(out.relative_to(ROOT)), n=N, n2=n2, n4=n4, pairs=npair)
