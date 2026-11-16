"""STEP 7e — Slide figure: LVOT ellipticity by depth and sex, two panels, in two build steps.

One method (minor-major axis ellipse: longest diameter through the centre / the diameter at 90 deg, as in a CT
reading), matched same-age, same-BMI pairs (one woman's and one man's mesh per pair).
  A  ellipticity vs depth below the annulus (annulus, 2 mm, 4 mm chord-closed): median line and IQR band
  B  paired difference women - men at each depth: mean and 95% CI (t), band = IQR of the pair differences
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


def q(v):
    return np.percentile(v, [25, 50, 75])


# panel A: median and IQR by depth, pooled and by sex
PROF = {g: np.array([q(A[(A.level == lv) & ((A.sex == g) if g != "all" else True)].ellip.dropna()) for lv, _, _ in LEVELS])
        for g in ("all", "F", "M")}
# panel B: paired women - men (same age and BMI) at each depth
DIFF = []
for lv, dp, _ in LEVELS:
    s = A[A.level == lv]
    p = s[s.sex == "F"].merge(s[s.sex == "M"], on=["age", "bmi"], suffixes=("_F", "_M"))[["ellip_F", "ellip_M"]].dropna()
    dd = p.ellip_F - p.ellip_M
    half = stats.t.ppf(.975, len(dd) - 1) * dd.std(ddof=1) / np.sqrt(len(dd))
    DIFF.append(dict(level=lv, depth=dp, n=len(dd), mean=dd.mean(), lo=dd.mean() - half, hi=dd.mean() + half,
                     q25=dd.quantile(.25), q75=dd.quantile(.75)))
DIFF = pd.DataFrame(DIFF)
for (lv, _, _), f, m in zip(LEVELS, PROF["F"], PROF["M"]):
    print(f"{lv:26s} median women {f[1]:.2f}  men {m[1]:.2f}")
print(DIFF.round(3).to_string(index=False))

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

    # A: ellipticity profile
    groups = [("all", t["grey"], "Women and\nmen")] if build == 1 else [("F", SEXC["F"], "Women"), ("M", SEXC["M"], "Men")]
    for g, col, lab in groups:
        P = PROF[g]
        a.fill_between(X, P[:, 0], P[:, 2], color=col, alpha=.18, lw=0)
        a.plot(X, P[:, 1], "-o", color=col, lw=2.6, ms=8, zorder=3)
        a.text(X[-1] + .18, P[-1, 1], lab, color=col, fontsize=15, fontweight="bold", va="center")
    a.axhline(1.0, color=t["muted"], lw=1, ls=":")
    a.text(-.4, 1.005, "1.0 = circle", fontsize=12, color=t["muted"], va="bottom")
    a.set_ylim(.97, 1.68)
    a.set_ylabel("Ellipticity (major / minor diameter)", fontsize=15, color=t["fg"])
    a.set_title("A  More elliptical below the annulus", loc="left", fontsize=17, fontweight="bold", color=t["fg"])
    style(a, t)
    a.set_xlim(-.5, 5.3)

    # B: paired difference women - men
    b.set_title("B  Women − men, same age and BMI", loc="left", fontsize=17, fontweight="bold", color=t["fg"])
    b.axhline(0, color=t["muted"], lw=1.2)
    b.set_ylim(-.04, .14)
    b.set_ylabel("Difference in ellipticity", fontsize=15, color=t["fg"])
    style(b, t)
    b.set_xlim(-.5, 4.9)
    if build == 1:
        b.set_visible(False)
    else:
        b.fill_between(X, DIFF.q25, DIFF.q75, color=SEXC["F"], alpha=.15, lw=0)
        b.plot(X, DIFF["mean"], "-", color=SEXC["F"], lw=2.2)
        b.errorbar(X, DIFF["mean"], yerr=[DIFF["mean"] - DIFF.lo, DIFF.hi - DIFF["mean"]], fmt="o", ms=8,
                   color=SEXC["F"], capsize=5, lw=2, zorder=3)
        for x_, r in zip(X, DIFF.itertuples()):
            b.text(x_ + .2, r.mean, f"{r.mean:+.2f}", ha="left", va="center", fontsize=14, fontweight="bold",
                   color=t["fg"])
        b.text(4.45, -.004, "0 = no difference", fontsize=12, color=t["muted"], ha="right", va="top")

    note = ("Minor–major axis ellipse; A: median and IQR" +
            ("; B: mean and 95% CI of the paired differences, band = IQR" if build == 2 else "") +
            f". {DIFF.n.min()}–{DIFF.n.max()} same-age, same-BMI pairs; 4 mm section closed with a chord; end-diastole.")
    fig.text(.01, .01, note, fontsize=10.5, color=t["muted"], ha="left", va="bottom")
    fig.subplots_adjust(left=.07, right=.98, top=.9, bottom=.2)
    out = ODIR / f"ellipticity_depth_{theme}_build{build}.png"
    fig.savefig(out, dpi=200, transparent=t["bg"] is None)
    plt.close(fig)
    return out


outs = [draw(th, bd) for th in THEMES for bd in (1, 2)]
for o in outs:
    print("wrote", o.relative_to(ROOT))
log_run(__file__, outputs=[str(o.relative_to(ROOT)) for o in outs],
        diff={r.level: [round(r.mean, 3), round(r.lo, 3), round(r.hi, 3), int(r.n)] for r in DIFF.itertuples()})
