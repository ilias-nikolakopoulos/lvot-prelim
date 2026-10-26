"""STEP 7e — Slide figure: LVOT ellipticity by depth (x-axis) and sex, in two build steps.

One method (minor-major axis ellipse: longest diameter through the centre / the diameter at 90 deg, as in a CT
reading), matched same-age, same-BMI pairs. Two messages, one per build:
  build 1  grey boxes + arrow: the LVOT gets more elliptical below the annulus
  build 2  women vs men coloured + gap brackets: women's LVOT is more elliptical, more so deeper
Box = IQR, line = median, whiskers = 5th-95th percentile.
Themes: dark (transparent background, light text, for a dark slide) and light (white background).
Outputs: outputs/figures/slides/ellipticity_depth_{dark,light}_build{1,2}.png
Needs: outputs/tables/area_methods_by_mesh.csv (STEP 6.6).
Usage: python scripts/19_ellipticity_depth_slide.py
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import font_manager
from matplotlib.patches import FancyArrowPatch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src.config import load_config  # noqa: E402
from src.provenance import log_run  # noqa: E402

for f_ in sorted((ROOT / "assets/fonts").glob("Roboto-*.ttf")):
    font_manager.fontManager.addfont(str(f_))
plt.rcParams["font.family"] = "Roboto"
SEXC = load_config()["analysis"]["colours"]  # women #D55E00, men #0072B2 (CVD-validated pair)
LEVELS = [("Aortic annulus", "Aortic\nannulus"), ("LVOT 2 mm (closed)", "2 mm below\nannulus"),
          ("LVOT 4 mm (chord-closed)", "4 mm below\nannulus")]
THEMES = {
    "dark": dict(fg="#EEF3F6", muted="#AEBBC4", grey="#8C9AA4", grey_edge="#C9D3DA", med="#FFFFFF", bg=None),
    "light": dict(fg="#102D40", muted="#5A6B77", grey="#B7C1C8", grey_edge="#7C8B95", med="#102D40", bg="white"),
}

A = pd.read_csv(ROOT / "outputs/tables/area_methods_by_mesh.csv")
A = A[A.paired & A.level.isin([lv for lv, _ in LEVELS])].copy()
A["ellip"] = A.ct_long / A.ct_short
V = {(lv, sx): A[(A.level == lv) & (A.sex == sx)].ellip.dropna().values for lv, _ in LEVELS for sx in "FM"}
MED = {k: float(np.median(v)) for k, v in V.items()}
for lv, _ in LEVELS:
    print(f"{lv:26s} women {MED[(lv, 'F')]:.2f}  men {MED[(lv, 'M')]:.2f}  gap {MED[(lv, 'F')] - MED[(lv, 'M')]:+.2f}  "
          f"n={len(V[(lv, 'F')])}")

ODIR = ROOT / "outputs/figures/slides"
ODIR.mkdir(parents=True, exist_ok=True)
OFF, W = .2, .32


def draw(theme, build):
    t = THEMES[theme]
    fig, ax = plt.subplots(figsize=(13.33, 6.2))
    if t["bg"] is None:
        fig.patch.set_alpha(0); ax.patch.set_alpha(0)
    else:
        fig.patch.set_facecolor(t["bg"]); ax.set_facecolor(t["bg"])
    x = np.arange(len(LEVELS))
    for i, (lv, _) in enumerate(LEVELS):
        for sx, dx in (("F", -OFF), ("M", OFF)):
            col = t["grey"] if build == 1 else SEXC[sx]
            edge = t["grey_edge"] if build == 1 else SEXC[sx]
            bp = ax.boxplot([V[(lv, sx)]], positions=[i + dx], widths=W, whis=(5, 95), showfliers=False,
                            patch_artist=True, medianprops=dict(color=t["med"], lw=2.5),
                            whiskerprops=dict(color=t["muted"], lw=1.4), capprops=dict(color=t["muted"], lw=1.4),
                            boxprops=dict(edgecolor=edge, lw=1.6))
            bp["boxes"][0].set_facecolor(col)
            bp["boxes"][0].set_alpha(.55 if build == 1 else .8)
            ax.text(i + dx, np.percentile(V[(lv, sx)], 95) + .015, f"{MED[(lv, sx)]:.2f}", ha="center", va="bottom",
                    fontsize=15, fontweight="bold", color=t["fg"])
    # depth arrow (both builds; quieter in build 2)
    y0, y1 = 1.30, 1.70
    ax.add_patch(FancyArrowPatch((x[0] + .05, y0), (x[2] - .25, y1), connectionstyle="arc3,rad=-0.25",
                                 arrowstyle="-|>", mutation_scale=28, lw=3.2, color=t["fg"],
                                 alpha=1 if build == 1 else .35, zorder=1))
    ax.text(-0.05, 1.50, "more elliptical\nwith depth", fontsize=19, fontweight="bold", color=t["fg"],
            alpha=1 if build == 1 else .35, ha="center", va="center", linespacing=1.05)
    if build == 2:
        # direct labels on the first pair instead of a legend
        for sx, dx, lab in (("F", -OFF, "Women"), ("M", OFF, "Men")):
            ax.text(x[0] + dx, np.percentile(V[(LEVELS[0][0], sx)], 5) - .03, lab, ha="center", va="top",
                    fontsize=16, fontweight="bold", color=SEXC[sx])
        # gap brackets: women - men median at each level
        for i, (lv, _) in enumerate(LEVELS):
            mf, mm = MED[(lv, "F")], MED[(lv, "M")]
            gap = mf - mm
            xb = i + OFF + W / 2 + .07
            if abs(gap) < .005:
                ax.text(xb, mf, "no gap", ha="left", va="center", fontsize=15, color=t["fg"], fontweight="bold")
                continue
            ax.plot([xb, xb + .05, xb + .05, xb], [mf, mf, mm, mm], color=t["fg"], lw=2)
            ax.text(xb + .09, (mf + mm) / 2, f"women\n+{gap:.2f}", ha="left", va="center", fontsize=15,
                    fontweight="bold", color=SEXC["F"], linespacing=1.05)
    ax.axhline(1.0, color=t["muted"], lw=1, ls=":")
    ax.text(2.62, 1.0, "1.0 = circle", ha="right", va="bottom", fontsize=13, color=t["muted"])
    ax.set_xticks(x)
    ax.set_xticklabels([lab for _, lab in LEVELS], fontsize=17, fontweight="bold", color=t["fg"])
    ax.tick_params(axis="y", labelsize=15, colors=t["fg"])
    ax.tick_params(axis="x", length=0, pad=10)
    ax.set_ylabel("Ellipticity (major / minor diameter)", fontsize=16, color=t["fg"])
    ax.set_ylim(.97, 1.78); ax.set_xlim(-.6, 2.75)
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)
    for s_ in ("left", "bottom"):
        ax.spines[s_].set_color(t["muted"])
    ax.grid(axis="y", color=t["muted"], alpha=.2, lw=.8)
    fig.tight_layout()
    out = ODIR / f"ellipticity_depth_{theme}_build{build}.png"
    fig.savefig(out, dpi=200, transparent=t["bg"] is None)
    plt.close(fig)
    return out


outs = [draw(th, b) for th in THEMES for b in (1, 2)]
for o in outs:
    print("wrote", o.relative_to(ROOT))
log_run(__file__, outputs=[str(o.relative_to(ROOT)) for o in outs],
        medians={f"{lv}|{sx}": round(m, 3) for (lv, sx), m in MED.items()})
