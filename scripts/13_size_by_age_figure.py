"""STEP 6.10 — Annulus and LVOT (2 mm) diameter by age block, women vs men, in the house figure style.

Same data as deck slide 6: mean area-derived diameter of the common-support pairs (same age and BMI in
both sexes) per age block. Error bars = 95% CI of each sex's mean (t-based, over the pairs in the block, from
lvot_all.csv); labels = paired men - women difference with its 95% CI (sexdiff_by_age_block.csv).
Style: Roboto (assets/fonts), level colours (annulus yellow, 2 mm green), women lighter, men darker.
Output: outputs/figures/analysis/fig_size_by_age.png
Usage: python scripts/13_size_by_age_figure.py
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from matplotlib import font_manager
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src.provenance import log_run  # noqa: E402

for f_ in sorted((ROOT / "assets/fonts").glob("Roboto-*.ttf")):
    font_manager.fontManager.addfont(str(f_))
plt.rcParams["font.family"] = "Roboto"
NAME_COL = "#102D40"
PANELS = [("ao_d_area", "Aortic annulus", "#E6C229"), ("lvot2_d_area", "2 mm from annulus", "#009E73")]
ALPHA = {"F": .45, "M": 1.0}

B = pd.read_csv(ROOT / "outputs/tables/sexdiff_by_age_block.csv")
L = pd.read_csv(ROOT / "outputs/tables/lvot_all.csv")
P = L[L.sex == "F"].merge(L[L.sex == "M"], on=["age", "bmi"], suffixes=("_F", "_M"))
BLOCKS = [tuple(map(int, b.split("-"))) for b in B[B.measure == "ao_d_area"].age_block]


def block_ci(m, sx):
    """Mean and 95% CI half-width per age block, over the common-support pairs with both values present."""
    out = []
    for lo, hi in BLOCKS:
        q = P[(P.age >= lo) & (P.age <= hi)][[f"{m}_F", f"{m}_M"]].dropna()[f"{m}_{sx}"]
        out.append((q.mean(), stats.t.ppf(.975, len(q) - 1) * q.std(ddof=1) / np.sqrt(len(q))))
    return np.array(out)
fig, axs = plt.subplots(1, 2, figsize=(13, 5.6), sharey=True)
for ax, (m, title, col) in zip(axs, PANELS):
    d = B[B.measure == m].reset_index(drop=True)
    x = np.arange(len(d))
    for sx, lab in (("F", "Women"), ("M", "Men")):
        mc = block_ci(m, sx)
        y, ci = mc[:, 0], mc[:, 1]
        assert np.allclose(y, d[f"mean_{sx}"].values, atol=1e-3), "block means differ from the analysis table"
        ax.errorbar(x, y, yerr=ci, fmt="none", ecolor=col, alpha=min(1, ALPHA[sx] + .2), elinewidth=1.6, capsize=5,
                    capthick=1.6, zorder=2)
        ax.plot(x, y, "-o", color=col, alpha=ALPHA[sx], lw=3, ms=9, mec=col, mfc=col, label=lab, zorder=3)
        for xi, yi, ci_ in zip(x, y, ci):
            ax.text(xi, yi + (ci_ + .1 if sx == "M" else -ci_ - .1), f"{yi:.1f}", ha="center",
                    va="bottom" if sx == "M" else "top", fontsize=9, color=NAME_COL)
    for xi, (lo, hi, df, cl, ch) in enumerate(zip(d.mean_F, d.mean_M, d.diff_M_minus_F, d.lo, d.hi)):
        ax.text(xi, (lo + hi) / 2, f"+{df:.1f}\n({cl:.1f}–{ch:.1f})", ha="center", va="center", fontsize=9.5,
                color=NAME_COL, fontweight="bold", linespacing=1.2,
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=col, lw=1))
    ax.text(len(d) - .55, (d.mean_F.iloc[-1] + d.mean_M.iloc[-1]) / 2, "men − women (95% CI)", ha="left", va="center",
            fontsize=9, color="0.4", rotation=90)
    ax.set_xticks(x)
    ax.set_xticklabels([b.replace("-", "–") for b in d.age_block], fontsize=10.5, color=NAME_COL, fontweight="bold")
    ax.set_xlabel("Age, years", fontsize=11, color=NAME_COL)
    ax.set_xlim(-.4, len(d) - .3)
    ax.set_title(title, fontsize=14, fontweight="bold", color="black", pad=12)
    ax.grid(axis="y", alpha=.25, lw=.6)
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)
    ax.legend(handles=[Line2D([], [], color=col, alpha=ALPHA["F"], lw=3, marker="o", ms=8, label="Women"),
                       Line2D([], [], color=col, alpha=ALPHA["M"], lw=3, marker="o", ms=8, label="Men")],
              fontsize=10, frameon=False, loc="upper right", ncol=2)
axs[0].set_ylabel("Area-derived diameter, mm", fontsize=11, color=NAME_COL)
axs[0].set_ylim(16.8, 22.6)
fig.suptitle("Men's outflow is ~3 mm wider at every age", fontsize=16, fontweight="bold", y=0.99)
fig.text(0.5, 0.01, "Means of 622 same-age, same-BMI pairs per age block; error bars and values in parentheses = 95% confidence intervals",
         ha="center", fontsize=10, color="0.4")
fig.tight_layout(rect=(0, 0.04, 1, 0.95), w_pad=3)
out = ROOT / "outputs/figures/analysis/fig_size_by_age.png"
fig.savefig(out, dpi=150)
print("wrote", out)
log_run(__file__, output=str(out.relative_to(ROOT)))
