"""STEP 6.11 — Continuity-equation AVA with a circular vs elliptical LVOT area, for the same VTIs.

AVA = LVOT area x VTI_LVOT / VTI_AV. With both VTIs fixed, AVA scales with the LVOT area used, so
  AVA_elliptical / AVA_circular = k = A_elliptical / A_circular,  A_circular = pi/4 * D_PLAX^2.
Elliptical area: the measured outline area (equals the moment ellipse within ~1%); the minor-major axis ellipse
(longest diameter x the one at 90 deg) is reported as a CT-style alternative.
Simulation: each heart's k (matched same-age, same-BMI pairs) is applied to a range of circular AVAs (0.5-1.5 cm2);
k is unitless, so this does not depend on the meshes' (end-diastolic, smaller) absolute LVOT size. Reported: the
corrected AVA and the share of hearts reclassified out of severe (< 1.0 cm2), and the circular AVA that corresponds to
a corrected AVA of 1.0 cm2 (= 1/k).
Caveats: end-diastolic averaged meshes, no leaflets; echo LVOT is measured in mid-systole; AVA cut-offs were
validated with the circular method. Illustrative only, not a clinical correction.
Outputs: outputs/tables/ava_correction_factors.csv, ava_simulation_reclassification.csv,
         outputs/figures/analysis/fig_ava_simulation.png
Usage: python scripts/18_ava_simulation.py
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
from matplotlib.patches import Patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src.provenance import log_run  # noqa: E402

for f_ in sorted((ROOT / "assets/fonts").glob("Roboto-*.ttf")):
    font_manager.fontManager.addfont(str(f_))
plt.rcParams["font.family"] = "Roboto"
NAVY = "#102D40"
LEVEL_COL = {"Aortic annulus": "#E6C229", "LVOT 2 mm (closed)": "#009E73", "LVOT 4 mm (chord-closed)": "#7B4FA6"}
LEVEL_TITLE = {"Aortic annulus": "Aortic annulus", "LVOT 2 mm (closed)": "2 mm from annulus",
               "LVOT 4 mm (chord-closed)": "4 mm from annulus"}
ALPHA = {"F": .4, "M": .85}
T = ROOT / "outputs/tables"

L = pd.read_csv(T / "area_methods_by_mesh.csv")
L = L[L.paired].copy()
L["A_circle"] = L.echo_circle / 100  # mm2 -> cm2
L["A_true"] = L["true"] / 100
L["A_mm"] = L.ellipse_ct / 100
L["k_true"] = L.A_true / L.A_circle
L["k_mm"] = L.A_mm / L.A_circle

# correction factors
rows = []
for (lev, sx), g in L.groupby(["level", "sex"], sort=False):
    r = {"level": lev, "sex": sx, "n": len(g)}
    for k in ("k_true", "k_mm"):
        q = g[k].quantile([.05, .25, .5, .75, .95]).values
        r.update({f"{k}_p5": q[0], f"{k}_q1": q[1], f"{k}_median": q[2], f"{k}_q3": q[3], f"{k}_p95": q[4]})
    r["circ_AVA_equiv_to_1cm2"] = 1 / r["k_true_median"]
    r["circ_AVA_equiv_to_1p5cm2"] = 1.5 / r["k_true_median"]
    rows.append(r)
K = pd.DataFrame(rows)
K.round(3).to_csv(T / "ava_correction_factors.csv", index=False)

# reclassification sweep: a given circular AVA multiplied by each heart's k (unitless, so independent of LVOT size)
sim = []
for (lev, sx), g in L.groupby(["level", "sex"], sort=False):
    for v in np.round(np.arange(0.50, 1.501, 0.05), 2):
        c = v * g.k_true
        sim.append({"level": lev, "sex": sx, "ava_circular": v, "ava_elliptical_median": c.median(),
                    "ava_elliptical_q1": c.quantile(.25), "ava_elliptical_q3": c.quantile(.75),
                    "pct_no_longer_severe": 100 * (c >= 1.0).mean() if v < 1.0 else np.nan,
                    "pct_no_longer_moderate_plus": 100 * (c >= 1.5).mean() if v < 1.5 else np.nan})
SIM = pd.DataFrame(sim)
SIM.round(3).to_csv(T / "ava_simulation_reclassification.csv", index=False)

pd.set_option("display.width", 220)
print(K[["level", "sex", "n", "k_true_q1", "k_true_median", "k_true_q3", "k_mm_median", "circ_AVA_equiv_to_1cm2"]].round(2).to_string(index=False))
ex = SIM[SIM.ava_circular.isin([0.6, 0.7, 0.8, 0.9])]
print(ex[["level", "sex", "ava_circular", "ava_elliptical_median", "pct_no_longer_severe"]].round(2).to_string(index=False))

# ---------------- figure ----------------
fig, axs = plt.subplots(1, 3, figsize=(18, 6.4), gridspec_kw={"width_ratios": [1.05, 1, 1.15]})

# (a) correction factor k by level and sex
ax = axs[0]
for i, lev in enumerate(LEVEL_COL):
    for j, sx in enumerate("FM"):
        v = L[(L.level == lev) & (L.sex == sx)].k_true.values
        pos = i + (j - .5) * .36
        bp = ax.boxplot([v], positions=[pos], widths=.3, whis=(5, 95), showfliers=False, patch_artist=True,
                        medianprops=dict(color=NAVY, lw=2), whiskerprops=dict(color="0.4"), capprops=dict(color="0.4"),
                        boxprops=dict(edgecolor=LEVEL_COL[lev], lw=1.2))
        bp["boxes"][0].set_facecolor(LEVEL_COL[lev]); bp["boxes"][0].set_alpha(ALPHA[sx])
        ax.text(pos, np.percentile(v, 95) + .015, f"×{np.median(v):.2f}", ha="center", va="bottom", fontsize=10,
                color=NAVY, fontweight="bold")
ax.axhline(1, color="0.5", lw=.8, ls=":")
ax.set_xticks(range(3)); ax.set_xticklabels([LEVEL_TITLE[l].replace(" from", "\nfrom") for l in LEVEL_COL], fontsize=10.5,
                                             fontweight="bold", color=NAVY)
ax.set_ylim(.95, 1.85); ax.set_ylabel("AVA elliptical / AVA circular  (= area ratio)", fontsize=11, color=NAVY)
ax.set_title("1  Correction factor for the same VTIs", fontsize=13.5, fontweight="bold", color=NAVY, loc="left")
ax.legend(handles=[Patch(fc="0.55", alpha=ALPHA["F"], label="Women"), Patch(fc="0.55", alpha=ALPHA["M"], label="Men")],
          fontsize=10, frameon=False, loc="upper left")

# (b) circular vs elliptical AVA, 2 mm, median k with IQR band
ax = axs[1]
lev, col = "LVOT 2 mm (closed)", LEVEL_COL["LVOT 2 mm (closed)"]
x = np.linspace(0.4, 2.0, 100)
ax.axhspan(0, 1.0, color="#B94F00", alpha=.07, lw=0)
for sx in "FM":
    r = K[(K.level == lev) & (K.sex == sx)].iloc[0]
    ax.fill_between(x, x * r.k_true_q1, x * r.k_true_q3, color=col, alpha=ALPHA[sx] * .45, lw=0)
    ax.plot(x, x * r.k_true_median, color=col, alpha=ALPHA[sx] + .1, lw=3,
            label=f"{'Women' if sx == 'F' else 'Men'}: ×{r.k_true_median:.2f}")
ax.plot(x, x, "--", color="0.45", lw=1.2, label="no correction")
for yv in (1.0, 1.5):
    ax.axhline(yv, color="0.6", lw=.8); ax.axvline(yv, color="0.6", lw=.8)
rF = K[(K.level == lev) & (K.sex == "F")].iloc[0]; rM = K[(K.level == lev) & (K.sex == "M")].iloc[0]
ax.annotate(f"circular 0.80 cm²\nbecomes {0.8 * rF.k_true_median:.2f} (women), {0.8 * rM.k_true_median:.2f} (men)",
            xy=(0.8, 0.8 * rF.k_true_median), xytext=(1.05, 0.62), fontsize=9.5, color=NAVY,
            arrowprops=dict(arrowstyle="->", color=NAVY))
ax.plot([0.8, 0.8], [0.8, 0.8 * rF.k_true_median], color=NAVY, lw=1.2)
ax.text(1.98, 0.95, "corrected AVA < 1.0 cm² (severe)", fontsize=9, color="#B94F00", va="top", ha="right")
ax.set_xlim(.4, 2.0); ax.set_ylim(.4, 3.0)
ax.set_xlabel("AVA, circular LVOT (cm²)", fontsize=11, color=NAVY)
ax.set_ylabel("AVA, elliptical LVOT (cm²)", fontsize=11, color=NAVY)
ax.set_title("2  Same VTIs, 2 mm from annulus", fontsize=13.5, fontweight="bold", color=NAVY, loc="left")
ax.legend(fontsize=10, frameon=False, loc="upper left")
ax.grid(alpha=.2)

# (c) reclassification: for a given circular AVA, share of hearts whose corrected AVA is >= 1.0 cm2
ax = axs[2]
xa = np.round(np.arange(0.50, 1.001, 0.01), 2)
for lev_, ls in (("LVOT 2 mm (closed)", "-"), ("Aortic annulus", ":")):
    for sx in "FM":
        k = L[(L.level == lev_) & (L.sex == sx)].k_true.values
        pct = [100 * np.mean(v * k >= 1.0) for v in xa]
        ax.plot(xa, pct, ls, color=LEVEL_COL[lev_], alpha=ALPHA[sx] + .1, lw=3)
ax.axvline(0.8, color="0.6", lw=.8)
ax.set_xlabel("AVA measured with a circular LVOT (cm²)", fontsize=11, color=NAVY)
ax.set_ylabel("No longer severe after correction (%)", fontsize=11, color=NAVY)
ax.set_title("3  Who stops being 'severe' (< 1.0 cm²)?", fontsize=13.5, fontweight="bold", color=NAVY, loc="left")
ax.legend(handles=[Line2D([], [], color=LEVEL_COL["LVOT 2 mm (closed)"], lw=3, label="LVOT 2 mm"),
                   Line2D([], [], color=LEVEL_COL["Aortic annulus"], ls=":", lw=3, label="annulus"),
                   Line2D([], [], color="0.55", alpha=ALPHA["F"] + .1, lw=3, label="women (lighter)"),
                   Line2D([], [], color="0.55", alpha=ALPHA["M"] + .1, lw=3, label="men (darker)")],
          fontsize=10, frameon=False, loc="upper left")
ax.set_xlim(.5, 1.0); ax.set_ylim(-2, 102); ax.grid(alpha=.2)
for a_ in axs:
    for s_ in ("top", "right"):
        a_.spines[s_].set_visible(False)

fig.suptitle("Aortic valve area with a circular vs elliptical LVOT", fontsize=17, fontweight="bold", color=NAVY, y=0.995)
fig.tight_layout(rect=(0, 0, 1, 0.95), w_pad=3)
out = ROOT / "outputs/figures/analysis/fig_ava_simulation.png"
fig.savefig(out, dpi=140)
print("wrote", out)
log_run(__file__, n=len(L), outputs=["outputs/tables/ava_correction_factors.csv", "outputs/tables/ava_simulation_dvi.csv",
                                     str(out.relative_to(ROOT))])
