"""STEP 7c — Presentation flowchart of the whole project (16:9 PNG, house style).

Six stages left to right (data -> geometry -> LVOT sections -> measures -> QC -> analysis), the four key results
below, and a reproducibility strip. Roboto (assets/fonts); names #102D40; level/result colours as in the figures.
Output: outputs/figures/project_flowchart.png
Usage: python scripts/14_project_flowchart.py
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src.provenance import log_run  # noqa: E402

for f_ in sorted((ROOT / "assets/fonts").glob("Roboto-*.ttf")):
    font_manager.fontManager.addfont(str(f_))
plt.rcParams["font.family"] = "Roboto"
NAVY, ORANGE, TINT, INK = "#102D40", "#B94F00", "#EEF3F6", "#1E2A33"

STAGES = [
    ("Data", ["1,423 averaged meshes", "UK Biobank (Ugurlu 2025)", "averaged by sex, age 46–83", "and BMI; end-diastole"]),
    ("Geometry", ["LV cavity = LV wall +", "mitral & aortic caps", "aortic & mitral rims,", "plane fits, LV long axis"]),
    ("LVOT sections", ["planes perpendicular to", "the aortic axis, 1 mm apart", "exact cut gives outline", "closed tube only 2–5 mm"]),
    ("Measures", ["area (shoelace formula)", "diameters, ellipticity", "(3 ellipse methods)", "echo PLAX circle"]),
    ("Quality control", ["1,423 / 1,423, 0 errors", "2 mm section: 1,413", "4 mm chord-closed: all", "renders + outlier review"]),
    ("Analysis", ["sex × age, BMI held fixed", "622 matched pairs", "+ regression", "indexed to LV size"]),
]
RESULTS = [
    ("~3 mm", "men's annulus and LVOT\nlarger at every age (15–16%)", "#0072B2"),
    ("1.14 vs 1.4–1.5", "round annulus,\noval LVOT just below", "#009E73"),
    ("+0.07 to +0.10", "women's LVOT\nmore oval than men's", "#7B4FA6"),
    ("−30%", "echo circle underestimates\nLVOT area (women −33%)", ORANGE),
]

fig = plt.figure(figsize=(16, 9))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 16); ax.set_ylim(0, 9); ax.axis("off")
ax.text(0.6, 8.35, "LVOT size and ellipticity in UK Biobank hearts: what we did", fontsize=24, fontweight="bold",
        color=NAVY, va="center")

# stage cards
n, x0, gap, y0, h = len(STAGES), 0.6, 0.35, 4.55, 3.0
w = (16 - 2 * x0 - (n - 1) * gap) / n
for i, (title, lines) in enumerate(STAGES):
    x = x0 + i * (w + gap)
    ax.add_patch(FancyBboxPatch((x, y0), w, h, boxstyle="round,pad=0,rounding_size=0.15", fc=TINT, ec="none"))
    ax.add_patch(FancyBboxPatch((x, y0 + h - 0.75), w, 0.75, boxstyle="round,pad=0,rounding_size=0.15", fc=NAVY, ec="none"))
    ax.add_patch(plt.Rectangle((x, y0 + h - 0.75), w, 0.3, fc=NAVY, ec="none"))  # square off the header's lower corners
    ax.text(x + 0.18, y0 + h - 0.375, f"{i + 1}", fontsize=16, fontweight="bold", color=ORANGE, va="center")
    ax.text(x + 0.48, y0 + h - 0.375, title, fontsize=15, fontweight="bold", color="white", va="center")
    for j, line in enumerate(lines):
        ax.text(x + 0.18, y0 + h - 1.1 - j * 0.42, line, fontsize=11.5, color=INK, va="center")
    if i < n - 1:
        ax.add_patch(FancyArrowPatch((x + w + 0.03, y0 + h / 2), (x + w + gap - 0.03, y0 + h / 2),
                                     arrowstyle="-|>", mutation_scale=16, color=NAVY, lw=1.5))

# arrow from the stages down to the results
ax.add_patch(FancyArrowPatch((8, y0 - 0.05), (8, 3.75), arrowstyle="-|>", mutation_scale=20, color=NAVY, lw=2))
ax.text(8.2, 4.1, "key results", fontsize=12, color=NAVY, fontweight="bold", va="center")

# results
rgap, ry, rh = 0.35, 1.75, 1.95
rw = (16 - 2 * x0 - 3 * rgap) / 4
rx0 = x0
for i, (big, label, col) in enumerate(RESULTS):
    x = rx0 + i * (rw + rgap)
    ax.add_patch(FancyBboxPatch((x, ry), rw, rh, boxstyle="round,pad=0,rounding_size=0.15", fc="white", ec=col, lw=2.2))
    ax.text(x + rw / 2, ry + rh - 0.55, big, fontsize=24, fontweight="bold", color=col, ha="center", va="center")
    ax.text(x + rw / 2, ry + 0.62, label, fontsize=12, color=INK, ha="center", va="center", linespacing=1.3)

# reproducibility strip
ax.add_patch(FancyBboxPatch((0.6, 0.45), 14.8, 0.85, boxstyle="round,pad=0,rounding_size=0.15", fc=TINT, ec="none"))
ax.text(0.85, 0.875, "Reproducible:", fontsize=12.5, fontweight="bold", color=NAVY, va="center")
ax.text(2.2, 0.875, "one command rebuilds everything  ·  all parameters in one config file  ·  every step documented with "
        "its code  ·  decisions log  ·  run log  ·  git history", fontsize=12, color=INK, va="center")

out = ROOT / "outputs/figures/project_flowchart.png"
fig.savefig(out, dpi=150)
print("wrote", out)
log_run(__file__, output=str(out.relative_to(ROOT)))
