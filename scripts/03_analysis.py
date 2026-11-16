"""STEP 6 — Age x sex analysis of LVOT / annulus measures.

Input: outputs/tables/<name>.csv from 01_measure.py (default lvot_all).
Outputs (outputs/tables, outputs/figures/analysis):
  exclusions.csv            mesh_id, measure group, reason
  sexdiff_by_age_block.csv  common-support paired M-F differences by age block
  fig_raw.png, fig_common.png, fig_sexdiff.png
Sex comparison = matched pairs only (same age and BMI).
Meshes are averages of many participants with unknown counts: CIs / p-values are descriptive.
Usage: python scripts/03_analysis.py [--in lvot_all]
"""
import argparse
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

sys_path_root = Path(__file__).resolve().parents[1]
import sys  # noqa: E402
sys.path.insert(0, str(sys_path_root))
from src.config import load_config  # noqa: E402
from src.provenance import log_run  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
ap = argparse.ArgumentParser()
ap.add_argument("--in", dest="inp", default="lvot_all")
a = ap.parse_args()
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", message=".*(fragmented|lsuffix).*")
T = ROOT / "outputs/tables"
FIG = ROOT / "outputs/figures/analysis"
FIG.mkdir(parents=True, exist_ok=True)
A_ = load_config()["analysis"]
COL = A_["colours"]  # Okabe-Ito vermillion / blue (validated CVD-safe)
LAB = {"F": "Women", "M": "Men"}
BLOCKS = [tuple(b) for b in A_["age_blocks"]]

d = pd.read_csv(T / f"{a.inp}.csv")
# STEP 6.1 — exclusions (processing errors, implausible rings, failed sections), logged with reason
excl = []
err = T / f"{a.inp}_errors.csv"
if err.exists():
    for _, r in pd.read_csv(err).iterrows():
        excl.append((r.mesh_id, "all", f"processing error: {r.error}"))
flags = d["qc_flags"].astype("string").fillna("")
ao_cols = [c for c in d if c.startswith(("ao_", "lvot"))]
bad = flags.str.contains("ao_size_implausible")
for m in d.mesh_id[bad]:
    excl.append((m, "aortic+LVOT", "aortic ring size implausible"))
d.loc[bad, ao_cols] = np.nan
bad = flags.str.contains("mv_size_implausible")
for m in d.mesh_id[bad]:
    excl.append((m, "mitral", "mitral ring size implausible"))
d.loc[bad, [c for c in d if c.startswith("mv_")] + ["ao_mv_angle", "ao_mv_centre_dist"]] = np.nan
for m in d.mesh_id[d["lvot2_d_area"].isna() & d["ao_d_area"].notna()]:
    r = d.loc[d.mesh_id == m, "section_reasons"].iloc[0]
    excl.append((m, "LVOT 2 mm (strict)", "section at 2 mm not a closed tube: " + str(r)[:120]))
for m in d.mesh_id[d["lvot4cc_d_area"].isna() & d["ao_d_area"].notna()]:
    g, r = d.loc[d.mesh_id == m, ["lvot4cc_gap_deg", "section_reasons"]].iloc[0]
    why = f"gap {g:.0f} deg > 120" if g > 120 else "section at 4 mm not chord-closed: " + str(r)[:120]
    excl.append((m, "LVOT 4 mm (chord-closed)", why))
if flags.str.contains("surface_regions").any():
    for m in d.mesh_id[flags.str.contains("surface_regions")]:
        excl.append((m, "note (kept)", "unexpected number of surface components"))

# STEP 6.2 — indexed size: diameter / LV cavity volume^(1/3) and % of LV length
d = d.copy()
d["lv_cav_cuberoot"] = d["lv_cavity_ml"] ** (1 / 3)
for k in ("ao_d_area", "lvot2_d_area", "lvot4cc_d_area", "mv_d_area"):
    d[k + "_per_edv13"] = d[k] / d["lv_cav_cuberoot"]
    d[k + "_pct_lvlen"] = 100 * d[k] / d["lv_length"]

MEAS = [  # key, label
    ("ao_d_area", "Aortic annulus diameter (area-derived, mm)"),
    ("ao_dmax", "Aortic annulus major diameter (mm)"),
    ("ao_dmin", "Aortic annulus minor diameter (mm)"),
    ("ao_ellip_ratio", "Aortic annulus ellipticity (major/minor)"),
    ("lvot2_d_area", "LVOT 2 mm diameter (area-derived, mm)"),
    ("lvot2_dmax", "LVOT 2 mm major diameter (mm)"),
    ("lvot2_dmin", "LVOT 2 mm minor diameter (mm)"),
    ("lvot2_ellip_ratio", "LVOT 2 mm ellipticity (major/minor)"),
    ("lvot4cc_d_area", "LVOT 4 mm diameter, chord-closed (mm)"),
    ("lvot4cc_ellip_ratio", "LVOT 4 mm ellipticity, chord-closed"),
    ("ao_d_area_per_edv13", "Annulus diam / LV cavity vol^(1/3)"),
    ("lvot2_d_area_per_edv13", "LVOT 2 mm diam / LV cavity vol^(1/3)"),
    ("lvot4cc_d_area_per_edv13", "LVOT 4 mm diam / LV cavity vol^(1/3)"),
    ("ao_d_area_pct_lvlen", "Annulus diam, % of LV length"),
    ("lvot2_d_area_pct_lvlen", "LVOT 2 mm diam, % of LV length"),
    ("mv_d_area", "Mitral ring diameter (area-derived, mm)"),
    ("ao_mv_angle", "Aorto-mitral plane angle (deg)"),
    ("ao_axis_vs_lv_long_angle", "Aortic axis vs LV long axis (deg)"),
    ("lv_cavity_ml", "LV cavity volume (mL)"),
    ("lv_length", "LV length, mitral centre-apex (mm)"),
]
MEAS = [(k, l) for k, l in MEAS if k in d and d[k].notna().sum() > 3]

# STEP 6.3 — common support: F/M pairs with the same age and BMI (holds BMI fixed)
F = d[d.sex == "F"].set_index(["age", "bmi"])
M = d[d.sex == "M"].set_index(["age", "bmi"])
pairs = F.join(M, lsuffix="_F", rsuffix="_M", how="inner").reset_index().copy()
print(f"{len(d)} meshes, {len(pairs)} common-support F/M pairs")


tq = lambda dof: stats.t.ppf(0.975, dof)
block_rows = []
# STEP 6.4 — per measure: paired M-F differences by age block (matched pairs)
for k, lab in MEAS:
    # common-support blocks
    for lo, hi in BLOCKS:
        p = pairs[(pairs.age >= lo) & (pairs.age <= hi)][[k + "_F", k + "_M"]].dropna()
        if len(p) < 2:
            continue
        diff = p[k + "_M"] - p[k + "_F"]
        se = diff.std(ddof=1) / np.sqrt(len(p))
        block_rows.append({"measure": k, "age_block": f"{lo}-{hi}", "n_pairs": len(p),
                           "mean_F": p[k + "_F"].mean(), "mean_M": p[k + "_M"].mean(),
                           "diff_M_minus_F": diff.mean(), "lo": diff.mean() - tq(len(p) - 1) * se,
                           "hi": diff.mean() + tq(len(p) - 1) * se,
                           "pct_diff": 100 * diff.mean() / p[k + "_F"].mean()})
blocks = pd.DataFrame(block_rows, columns=["measure", "age_block", "n_pairs", "mean_F", "mean_M", "diff_M_minus_F", "lo", "hi", "pct_diff"])
blocks.round(4).to_csv(T / "sexdiff_by_age_block.csv", index=False)
pd.DataFrame(excl, columns=["mesh_id", "measure_group", "reason"]).to_csv(T / "exclusions.csv", index=False)


def grid(n):
    nc = 4
    nr = int(np.ceil(n / nc))
    fig, axs = plt.subplots(nr, nc, figsize=(4.2 * nc, 3.2 * nr), squeeze=False)
    for ax in axs.ravel()[n:]:
        ax.axis("off")
    return fig, axs.ravel()


def style(ax, lab):
    ax.set_title(lab, fontsize=9)
    ax.set_xlabel("Age (years)", fontsize=8)
    ax.tick_params(labelsize=7)
    ax.grid(alpha=.25, lw=.5)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def roll(x, y, w=3):
    s = pd.Series(np.asarray(y, float), index=np.asarray(x)).groupby(level=0).mean().sort_index()
    return s.index, s.rolling(w, center=True, min_periods=1).mean().to_numpy()


# STEP 6.5 — figures
# 1 raw: all meshes, per-age mean (3-yr rolling)
fig, axs = grid(len(MEAS))
for ax, (k, lab) in zip(axs, MEAS):
    for sx in "FM":
        q = d[(d.sex == sx) & d[k].notna()]
        ax.scatter(q.age, q[k], s=4, color=COL[sx], alpha=.25, lw=0)
        x, y = roll(q.age, q[k])
        ax.plot(x, y, color=COL[sx], lw=2, label=LAB[sx])
    style(ax, lab)
axs[0].legend(fontsize=8, frameon=False)
fig.suptitle("All meshes (BMI not held fixed): points = meshes, line = 3-year rolling mean of per-age means", fontsize=10)
fig.tight_layout(rect=(0, 0, 1, 0.97)); fig.savefig(FIG / "fig_raw.png", dpi=120); plt.close(fig)

# 2 common support
fig, axs = grid(len(MEAS))
for ax, (k, lab) in zip(axs, MEAS):
    for sx in "FM":
        q = pairs[["age", f"{k}_F", f"{k}_M"]].dropna()
        g = q.groupby("age")[f"{k}_{sx}"]
        ax.scatter(g.mean().index, g.mean(), s=10, color=COL[sx], alpha=.5, lw=0)
        x, y = roll(q.age, q[f"{k}_{sx}"])
        ax.plot(x, y, color=COL[sx], lw=2, label=LAB[sx])
    style(ax, lab)
axs[0].legend(fontsize=8, frameon=False)
fig.suptitle("BMI held fixed: common support only (same age and BMI in both sexes); dots = per-age mean, line = 3-yr rolling", fontsize=10)
fig.tight_layout(rect=(0, 0, 1, 0.97)); fig.savefig(FIG / "fig_common.png", dpi=120); plt.close(fig)

# 3 sex difference vs age
fig, axs = grid(len(MEAS))
for ax, (k, lab) in zip(axs, MEAS):
    q = pairs[["age", f"{k}_F", f"{k}_M"]].dropna()
    q = q.assign(diff=q[f"{k}_M"] - q[f"{k}_F"])
    g = q.groupby("age")["diff"]
    ax.scatter(g.mean().index, g.mean(), s=10, color="0.35", lw=0, label="per-age mean of pairs")
    bl = blocks[blocks.measure == k]
    for _, r in bl.iterrows():
        lo, hi = map(int, r.age_block.split("-"))
        ax.errorbar((lo + hi) / 2, r.diff_M_minus_F, yerr=[[r.diff_M_minus_F - r.lo], [r.hi - r.diff_M_minus_F]],
                    fmt="s", ms=4, color="#009E73", capsize=2, lw=1, label="age block (paired) ± 95% CI" if lo == 46 else None)
    ax.axhline(0, color="0.5", lw=.8)
    style(ax, "M − F: " + lab)
axs[0].legend(fontsize=7, frameon=False)
fig.suptitle("Sex difference (men − women) by age", fontsize=10)
fig.tight_layout(rect=(0, 0, 1, 0.97)); fig.savefig(FIG / "fig_sexdiff.png", dpi=120); plt.close(fig)

pd.set_option("display.width", 250)
show = blocks[blocks.measure.isin(["ao_d_area", "lvot2_d_area", "lvot4cc_d_area", "lvot2_ellip_ratio", "lvot4cc_ellip_ratio"])]
print(show.round(3).to_string(index=False))
print(f"exclusions logged: {len(excl)}")
log_run(__file__, input=a.inp, n_meshes=len(d), n_pairs=len(pairs), n_exclusions=len(excl),
        outputs=["outputs/tables/sexdiff_by_age_block.csv",
                 "outputs/tables/exclusions.csv", "outputs/figures/analysis/"])
