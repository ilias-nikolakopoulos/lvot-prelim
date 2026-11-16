"""STEP 6.10 — Annulus and LVOT (2 mm) diameter by age block, women vs men, in the house figure style.

Common-support pairs (same age and BMI in both sexes) per age block, from lvot_all.csv:
points = median area-derived diameter per sex, error bars = interquartile range (25th-75th percentile);
labels = median (IQR) of the paired men - women differences.
Style: Roboto (assets/fonts), level colours (annulus yellow, 2 mm green), women lighter, men darker.
Outputs: outputs/figures/analysis/fig_size_by_age.png (area-derived diameter, mm; default)
         outputs/figures/analysis/fig_area_by_age.png (outline area, mm2; with --area)
Usage: python scripts/13_size_by_age_figure.py [--area]
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import font_manager
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src.provenance import log_run  # noqa: E402

for f_ in sorted((ROOT / "assets/fonts").glob("Roboto-*.ttf")):
    font_manager.fontManager.addfont(str(f_))
plt.rcParams["font.family"] = "Roboto"
NAME_COL = "#102D40"
AREA = "--area" in sys.argv
SUF, UNIT, NDP = ("area", "mm²", 0) if AREA else ("d_area", "mm", 1)
PANELS = [(f"ao_{SUF}", "Aortic annulus", "#E6C229"), (f"lvot2_{SUF}", "2 mm from annulus", "#009E73")]
ALPHA = {"F": .45, "M": 1.0}

B = pd.read_csv(ROOT / "outputs/tables/sexdiff_by_age_block.csv")
L = pd.read_csv(ROOT / "outputs/tables/lvot_all.csv")
P = L[L.sex == "F"].merge(L[L.sex == "M"], on=["age", "bmi"], suffixes=("_F", "_M"))
BLOCKS = [tuple(map(int, b.split("-"))) for b in B[B.measure == "ao_d_area"].age_block]


def block_stats(m):
    """Per age block: median and IQR of each sex, and median (IQR) of the paired men - women difference."""
    rows = []
    for lo, hi in BLOCKS:
        q = P[(P.age >= lo) & (P.age <= hi)][[f"{m}_F", f"{m}_M"]].dropna()
        dif = q[f"{m}_M"] - q[f"{m}_F"]
        r = {"n": len(q)}
        for sx in "FM":
            r[f"q1_{sx}"], r[f"med_{sx}"], r[f"q3_{sx}"] = q[f"{m}_{sx}"].quantile([.25, .5, .75])
        r["d_q1"], r["d_med"], r["d_q3"] = dif.quantile([.25, .5, .75])
        rows.append(r)
    return pd.DataFrame(rows)


STATS = {m: block_stats(m) for m, _, _ in PANELS}
lo_ = min(st[["q1_F", "q1_M"]].min().min() for st in STATS.values())
hi_ = max(st[["q3_F", "q3_M"]].max().max() for st in STATS.values())
pad = (hi_ - lo_) * .12
OFF = (hi_ - lo_) * .02  # value-label offset
fig, axs = plt.subplots(1, 2, figsize=(13, 5.6), sharey=True)
for ax, (m, title, col) in zip(axs, PANELS):
    d = B[B.measure == "ao_d_area"].reset_index(drop=True)  # age-block labels
    st = STATS[m]
    x = np.arange(len(d))
    for sx, lab in (("F", "Women"), ("M", "Men")):
        y, q1, q3 = st[f"med_{sx}"].values, st[f"q1_{sx}"].values, st[f"q3_{sx}"].values
        ax.errorbar(x, y, yerr=[y - q1, q3 - y], fmt="none", ecolor=col, alpha=min(1, ALPHA[sx] + .2), elinewidth=1.6,
                    capsize=5, capthick=1.6, zorder=2)
        ax.plot(x, y, "-o", color=col, alpha=ALPHA[sx], lw=3, ms=9, mec=col, mfc=col, label=lab, zorder=3)
        for xi, yi, a1, a3 in zip(x, y, q1, q3):
            ax.text(xi, a3 + OFF if sx == "M" else a1 - OFF, f"{yi:.{NDP}f}", ha="center",
                    va="bottom" if sx == "M" else "top", fontsize=9, color=NAME_COL)
    for xi, r in st.iterrows():
        ax.text(xi, (r.q1_F + r.q3_M) / 2 if r.q3_F < r.q1_M else (r.med_F + r.med_M) / 2,
                f"+{r.d_med:.{NDP}f}\n({r.d_q1:.{NDP}f}–{r.d_q3:.{NDP}f})", ha="center", va="center", fontsize=9.5,
                color=NAME_COL, fontweight="bold", linespacing=1.2,
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=col, lw=1))
    ax.text(len(d) - .55, (st.med_F.iloc[-1] + st.med_M.iloc[-1]) / 2, "men − women, median (IQR)", ha="left", va="center",
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
axs[0].set_ylabel(f"{'Outline area' if AREA else 'Area-derived diameter'}, {UNIT} (median, IQR)", fontsize=11, color=NAME_COL)
axs[0].set_ylim(lo_ - pad * 1.0, hi_ + pad * 1.4)
dmed = np.median(np.concatenate([STATS[m].d_med.values for m, _, _ in PANELS]))
TITLE = (f"Men's outflow area is ~{round(dmed, -1):.0f} mm² larger at every age" if AREA
         else f"Men's outflow is ~{dmed:.0f} mm wider at every age")
fig.suptitle(TITLE, fontsize=16, fontweight="bold", y=0.99)
fig.text(0.5, 0.01, "Medians of the same-age, same-BMI pairs per age block; error bars and values in parentheses = interquartile range",
         ha="center", fontsize=10, color="0.4")
fig.tight_layout(rect=(0, 0.04, 1, 0.95), w_pad=3)
out = ROOT / "outputs/figures/analysis" / ("fig_area_by_age.png" if AREA else "fig_size_by_age.png")
fig.savefig(out, dpi=150)
print("wrote", out)
log_run(__file__, output=str(out.relative_to(ROOT)))
