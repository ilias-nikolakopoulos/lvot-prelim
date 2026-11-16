"""STEP 7e — Slide figure: LVOT ellipticity by depth and sex, two panels, in two build steps.

One method (minor-major axis ellipse: longest diameter through the centre / the diameter at 90 deg, as in a CT
reading), matched same-age, same-BMI pairs (one woman's and one man's mesh per pair).
  A  ellipticity vs depth below the annulus (annulus, 2 mm, 4 mm chord-closed): box plots, line through the medians
  B  paired difference women - men at each depth: box plots, line through the medians; mean and 95% CI (t) overlaid
Box = IQR, line in box = median, whiskers = 5th-95th percentile.
Builds (same layout, so the slide does not jump):
  build 1  A only, both sexes pooled (grey): the LVOT gets more elliptical below the annulus
  build 2  A by sex (coloured) + B: women's LVOT is more elliptical, and the gap grows with depth
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
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src.config import load_config  # noqa: E402
from src.provenance import log_run  # noqa: E402

for f_ in sorted((ROOT / "assets/fonts").glob("Roboto-*.ttf")):
    font_manager.fontManager.addfont(str(f_))
plt.rcParams["font.family"] = "Roboto"
SEXC = load_config()["analysis"]["colours"]  # women #D55E00, men #0072B2 (CVD-validated pair)
LEVELS = [("Aortic annulus", 0, "Annulus"), ("LVOT 2 mm (closed)", 2, "2 mm"), ("LVOT 4 mm (chord-closed)", 4, "4 mm")]
THEMES = {
    "dark": dict(fg="#EEF3F6", muted="#AEBBC4", grey="#C9D3DA", bg=None),
    "light": dict(fg="#102D40", muted="#5A6B77", grey="#7C8B95", bg="white"),
}

A = pd.read_csv(ROOT / "outputs/tables/area_methods_by_mesh.csv")
A = A[A.paired & A.level.isin([lv for lv, _, _ in LEVELS])].copy()
A["ellip"] = A.ct_long / A.ct_short
X = np.array([dp for _, dp, _ in LEVELS], float)


# panel A: values by level, pooled and by sex
V = {g: [A[(A.level == lv) & ((A.sex == g) if g != "all" else True)].ellip.dropna().values for lv, _, _ in LEVELS]
     for g in ("all", "F", "M")}
# panel B: paired women - men (same age and BMI) at each level
D, DIFF = [], []
for lv, dp, _ in LEVELS:
    s_ = A[A.level == lv]
    p = s_[s_.sex == "F"].merge(s_[s_.sex == "M"], on=["age", "bmi"], suffixes=("_F", "_M"))[["ellip_F", "ellip_M"]].dropna()
    dd = (p.ellip_F - p.ellip_M).values
    half = stats.t.ppf(.975, len(dd) - 1) * dd.std(ddof=1) / np.sqrt(len(dd))
    D.append(dd)
    DIFF.append(dict(level=lv, depth=dp, n=len(dd), median=np.median(dd), mean=dd.mean(), lo=dd.mean() - half,
                     hi=dd.mean() + half))
DIFF = pd.DataFrame(DIFF)
for (lv, _, _), f, m in zip(LEVELS, V["F"], V["M"]):
    print(f"{lv:26s} median women {np.median(f):.2f}  men {np.median(m):.2f}")
print(DIFF.round(3).to_string(index=False))


def boxes(ax, data, xs, col, t, w):
    bp = ax.boxplot(data, positions=xs, widths=w, whis=(5, 95), showfliers=False, patch_artist=True,
                    medianprops=dict(color=t["fg"], lw=2.2), whiskerprops=dict(color=col, lw=1.4),
                    capprops=dict(color=col, lw=1.4), boxprops=dict(edgecolor=col, lw=1.6), manage_ticks=False)
    for bx in bp["boxes"]:
        bx.set_facecolor(col); bx.set_alpha(.55)
    ax.plot(xs, [np.median(v) for v in data], "-", color=col, lw=2.4, zorder=3)


ODIR = ROOT / "outputs/figures/slides"
ODIR.mkdir(parents=True, exist_ok=True)


def style(ax, t):
    ax.set_xticks(X, [lab for _, _, lab in LEVELS])
    ax.set_xlim(-.5, 4.5)
    ax.tick_params(labelsize=14, colors=t["fg"])
    ax.set_xlabel("Level below the aortic annulus", fontsize=15, color=t["fg"])
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)
    for s_ in ("left", "bottom"):
        ax.spines[s_].set_color(t["muted"])
    ax.grid(axis="y", color=t["muted"], alpha=.2, lw=.8)
    if t["bg"] is None:
        ax.patch.set_alpha(0)


def draw(theme, build):
    t = THEMES[theme]
    fig, (a, b) = plt.subplots(1, 2, figsize=(13.33, 6.2), gridspec_kw=dict(width_ratios=[1.25, 1], wspace=.32))
    if t["bg"] is None:
        fig.patch.set_alpha(0)
    else:
        fig.patch.set_facecolor(t["bg"])

    # A: ellipticity by level, box plots + line through the medians
    groups = [("all", t["grey"], "Women and\nmen", 0)] if build == 1 else \
        [("F", SEXC["F"], "Women", -.38), ("M", SEXC["M"], "Men", .38)]
    for k, (g, col, lab, dx) in enumerate(groups):
        boxes(a, V[g], X + dx, col, t, .62)
        a.text(-.6, 1.66 - .065 * k, lab.replace("\n", " "), color=col, fontsize=15, fontweight="bold", va="center")
    a.axhline(1.0, color=t["muted"], lw=1, ls=":")
    a.text(-.75, 1.005, "1.0 = circle", fontsize=12, color=t["muted"], va="bottom")
    a.set_ylim(.97, 1.72)
    a.set_ylabel("Ellipticity (major / minor diameter)", fontsize=15, color=t["fg"])
    a.set_title("A  More elliptical below the annulus", loc="left", fontsize=17, fontweight="bold", color=t["fg"])
    style(a, t)
    a.set_xlim(-.8, 4.8)

    # B: paired difference women - men, box plots + line through the medians, mean and 95% CI overlaid
    b.set_title("B  Women − men, same age and BMI", loc="left", fontsize=17, fontweight="bold", color=t["fg"])
    b.axhline(0, color=t["muted"], lw=1.2)
    b.set_ylim(-.09, .22)
    b.set_ylabel("Difference in ellipticity", fontsize=15, color=t["fg"])
    style(b, t)
    b.set_xlim(-.7, 5.1)
    if build == 1:
        b.set_visible(False)
    else:
        boxes(b, D, X, SEXC["F"], t, .62)
        b.errorbar(X, DIFF["mean"], yerr=[DIFF["mean"] - DIFF.lo, DIFF.hi - DIFF["mean"]], fmt="D", ms=6,
                   mfc="white", mec=t["fg"], color=t["fg"], capsize=3, lw=1.4, zorder=4)
        for x_, r in zip(X, DIFF.itertuples()):
            b.text(x_ + .36, r.median, f"{r.median:+.2f}", ha="left", va="center", fontsize=14, fontweight="bold",
                   color=t["fg"])
        b.text(3.0, -.006, "0 = no difference", fontsize=12, color=t["muted"], ha="center", va="top")

    note = ("Minor–major axis ellipse; box = IQR, whiskers = 5th–95th percentile, line = medians" +
            ("; B: paired differences, white diamond = mean and 95% CI, numbers = medians" if build == 2 else "") + ".\n" +
            f"{DIFF.n.min()}–{DIFF.n.max()} same-age, same-BMI pairs; 4 mm section closed with a chord; end-diastole.")
    fig.text(.01, .01, note, fontsize=10.5, color=t["muted"], ha="left", va="bottom")
    fig.subplots_adjust(left=.07, right=.98, top=.9, bottom=.22)
    out = ODIR / f"ellipticity_depth_{theme}_build{build}.png"
    fig.savefig(out, dpi=200, transparent=t["bg"] is None)
    plt.close(fig)
    return out


outs = [draw(th, bd) for th in THEMES for bd in (1, 2)]
for o in outs:
    print("wrote", o.relative_to(ROOT))
log_run(__file__, outputs=[str(o.relative_to(ROOT)) for o in outs],
        diff={r.level: [round(r.mean, 3), round(r.lo, 3), round(r.hi, 3), int(r.n)] for r in DIFF.itertuples()})
