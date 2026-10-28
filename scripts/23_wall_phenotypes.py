"""STEP 9.6-9.8 — Which cardiac features go with a more elliptical LVOT, and two phenotypes.

Every feature and the LVOT ellipticity (minor-major axis ellipse at 2 and 4 mm) is turned into a z-score within sex:
residual after a fit on age and BMI (cubic surface with >= 200 meshes per sex, linear otherwise), divided by its SD.
So z = how unusual a mesh is for its sex, age and BMI.
  9.6 feature map: Spearman rho of each feature's z with the ellipticity z (women, men, both)
  9.7 phenotypes (config wall.normal_z / abnormal_z):
      septum z = mean z of the septal measures; other walls: mean and max z of the 9 non-septal segments
      normal walls         |septum z| and |other mean z| < normal_z, other max z < abnormal_z
      reference            normal walls and all angles |z| < normal_z
      A thick elsewhere    |septum z| < normal_z, some non-septal segment z >= abnormal_z
      B angle low / high   normal walls, one angle z <= -abnormal_z / >= abnormal_z
      C thick septum       septum z >= abnormal_z (for comparison)
      B (z 1.0)            sensitivity: as B with the cut-off at wall.angle_cut_sensitivity_z instead of abnormal_z
      ellipticity z in each group vs reference: median difference, bootstrap 95% CI, Mann-Whitney p
  9.7b angles as continuous measures, normal-wall meshes only (no cut-off): Spearman rho (women, men, both) and OLS
      ellipticity z ~ angle z + female, then + aortic annulus diameter z, septum z, other-wall mean z
  9.8 figure: feature map
Caveat: representative (cohort-average) meshes; z describes residual shape variation, not individual patients.
Outputs: outputs/tables/<inp>_feature_map.csv, <inp>_phenotypes.csv, <inp>_angle_continuous.csv, outputs/figures/analysis/fig_<inp>_feature_map.png
Usage: python scripts/23_wall_phenotypes.py [--inp wall_all|wall_subset]
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src.config import load_config  # noqa: E402
from src.provenance import log_run  # noqa: E402
from src.wall import SECTORS  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--inp", default="wall_all")
a = ap.parse_args()
cfg = load_config()
W_, SEXC = cfg["wall"], cfg["analysis"]["colours"]
T = ROOT / "outputs/tables"
LEV = list(W_["levels"])
SEPTAL = ["wt_subaortic_septum"] + [f"wt_{lv}_septal" for lv in LEV]
OTHER = [f"wt_{lv}_{s}" for lv in LEV for s in SECTORS if s != "septal"]
ANGLES = {"ao_mv_angle": "Aorto-mitral plane angle", "aosa_deg": "Aorto-septal angle",
          "ao_axis_vs_lv_long_angle": "Aortic axis vs LV long axis"}
FEAT = ({c: c.replace("wt_", "").replace("_", " ") + " thickness" for c in SEPTAL + OTHER}
        | {"sep_bulge_ratio": "Septal bulge (subaortic / mid septum)", "sep_lat_ratio_basal": "Septal / lateral, basal",
           "sep_lat_ratio_mid": "Septal / lateral, mid", "wt_max_p95": "Max wall thickness (95th pct)"}
        | ANGLES
        | {"ao_mv_centre_dist": "Aortic-mitral centre distance", "mv_ellip_ratio": "Mitral annulus ellipticity",
           "ao_d_area": "Aortic annulus diameter", "lv_myo_ml": "LV myocardial volume", "lv_cavity_ml": "LV EDV",
           "sphericity": "LV sphericity"})

w = pd.read_csv(T / f"{a.inp}.csv")
d = pd.read_csv(T / "lvot_all.csv")
fl = d.qc_flags.astype("string").fillna("")
d.loc[fl.str.contains("ao_size_implausible"), ["ao_d_area", "ao_axis_vs_lv_long_angle", "ao_mv_angle", "ao_mv_centre_dist"]] = np.nan
d.loc[fl.str.contains("mv_size_implausible"), ["mv_ellip_ratio", "ao_mv_angle", "ao_mv_centre_dist"]] = np.nan
d = d.copy()
d["sphericity"] = d.lv_cavity_ml * 1000 / (np.pi / 6 * d.lv_length ** 3)
A = pd.read_csv(T / "area_methods_by_mesh.csv")
A = A[A.level.str.startswith("LVOT")].assign(e=lambda x: x.ct_long / x.ct_short)
E = A.pivot(index="mesh_id", columns="level", values="e")
E.columns = ["ellip_2mm" if "2 mm" in c else "ellip_4mm" for c in E.columns]
D = w.merge(d[["mesh_id"] + [c for c in FEAT if c in d]], on="mesh_id").merge(E, left_on="mesh_id", right_index=True)
OUT = ["ellip_2mm", "ellip_4mm"]
print(f"{len(D)} meshes ({(D.sex == 'F').sum()} F / {(D.sex == 'M').sum()} M)")


def zscore(df, c):
    x = df[["age", "bmi", c]].dropna()
    u, v = (x.age - 65) / 10, (x.bmi - 27) / 5
    cols = [np.ones(len(x)), u, v]
    if len(x) >= 200:
        cols += [u * u, v * v, u * v, u ** 3, v ** 3, u * u * v, u * v * v]
    X = np.column_stack(cols)
    r = x[c] - X @ np.linalg.lstsq(X, x[c], rcond=None)[0]
    return r / r.std()


Z = pd.concat([pd.DataFrame({c: zscore(D[D.sex == s], c) for c in list(FEAT) + OUT}) for s in "FM"])
Z = D[["mesh_id", "sex"]].join(Z)

# 9.6 feature map
rows = []
for c, lab in FEAT.items():
    r = {"feature": c, "label": lab}
    for o in OUT:
        for s, sub in (("F", Z[Z.sex == "F"]), ("M", Z[Z.sex == "M"]), ("all", Z)):
            x = sub[[c, o]].dropna()
            r[f"rho_{o}_{s}"], r[f"p_{o}_{s}"] = stats.spearmanr(x[c], x[o]) if len(x) > 5 else (np.nan, np.nan)
    rows.append(r)
MAP = pd.DataFrame(rows).sort_values("rho_ellip_4mm_all")

# 9.7 phenotypes
nz, az = W_["normal_z"], W_["abnormal_z"]
sep_z = Z[SEPTAL].mean(1)                      # septum: mean z of its measures
oth_mean, oth_max = Z[OTHER].mean(1), Z[OTHER].max(1)
sep_ok = sep_z.abs() < nz
walls_ok = sep_ok & (oth_mean.abs() < nz) & (oth_max < az)
ang_ok = (Z[list(ANGLES)].abs() < nz).all(1)
G = {"reference (all normal)": walls_ok & ang_ok,
     "A: normal septum, thicker elsewhere": sep_ok & (oth_max >= az),
     "C: thicker septum": sep_z >= az}
for c, lab in ANGLES.items():
    G[f"B: normal walls, {lab} low"] = walls_ok & (Z[c] <= -az)
    G[f"B: normal walls, {lab} high"] = walls_ok & (Z[c] >= az)
sz = W_["angle_cut_sensitivity_z"]
for c, lab in ANGLES.items():
    G[f"B (z {sz}): normal walls, {lab} low"] = walls_ok & (Z[c] <= -sz)
    G[f"B (z {sz}): normal walls, {lab} high"] = walls_ok & (Z[c] >= sz)
rng = np.random.default_rng(cfg["seed"])
ref = G["reference (all normal)"]
rows = []
for g, m in G.items():
    r = {"group": g, "n": int(m.sum()), "n_F": int((m & (Z.sex == "F")).sum()), "n_M": int((m & (Z.sex == "M")).sum())}
    for o in OUT:
        x, y = Z.loc[m, o].dropna().values, Z.loc[ref, o].dropna().values
        raw = D.loc[m.values, o].dropna()
        r[f"{o}_raw_median"] = raw.median() if len(raw) else np.nan
        r[f"{o}_z_median"] = np.median(x) if len(x) else np.nan
        if g.startswith("reference") or len(x) < 5 or len(y) < 5:
            continue
        bs = [np.median(rng.choice(x, len(x))) - np.median(rng.choice(y, len(y))) for _ in range(2000)]
        r[f"{o}_z_diff_vs_ref"] = np.median(x) - np.median(y)
        r[f"{o}_lo"], r[f"{o}_hi"] = np.percentile(bs, [2.5, 97.5])
        r[f"{o}_p"] = stats.mannwhitneyu(x, y).pvalue
    rows.append(r)
PH = pd.DataFrame(rows)


# 9.7b angles as continuous measures within normal-wall meshes
def ols(y, X):
    X = np.column_stack([np.ones(len(y))] + X)
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    res = y - X @ b
    df_ = len(y) - X.shape[1]
    se = np.sqrt(np.diag(res @ res / df_ * np.linalg.inv(X.T @ X)))
    return b[1], se[1], 2 * stats.t.sf(abs(b[1] / se[1]), df_)


Zw = Z[walls_ok].assign(female=lambda x: (x.sex == "F").astype(float), septum=sep_z[walls_ok], other=oth_mean[walls_ok])
rows = []
for c, lab in ANGLES.items():
    for o in OUT:
        r = {"angle": lab, "outcome": o}
        for s, sub in (("F", Zw[Zw.sex == "F"]), ("M", Zw[Zw.sex == "M"]), ("all", Zw)):
            x = sub[[c, o]].dropna()
            r[f"n_{s}"] = len(x)
            r[f"rho_{s}"], r[f"rho_p_{s}"] = stats.spearmanr(x[c], x[o])
        for tag, adj in (("", []), ("_adj", ["ao_d_area", "septum", "other"])):
            x = Zw[[o, c, "female"] + adj].dropna()
            b, se, p = ols(x[o].values, [x[c].values, x.female.values] + [x[k].values for k in adj])
            r[f"slope{tag}"], r[f"lo{tag}"], r[f"hi{tag}"], r[f"p{tag}"], r[f"n{tag}"] = b, b - 1.96 * se, b + 1.96 * se, p, len(x)
        rows.append(r)
CT = pd.DataFrame(rows)

pd.set_option("display.width", 250, "display.max_columns", 30)
print(MAP[["label", "rho_ellip_2mm_F", "rho_ellip_2mm_M", "rho_ellip_4mm_F", "rho_ellip_4mm_M", "rho_ellip_4mm_all"]].round(2).to_string(index=False))
print(PH.round(3).to_string(index=False))
print(f"normal-wall meshes: {int(walls_ok.sum())}")
print(CT.round(3).to_string(index=False))
CT.to_csv(T / f"{a.inp}_angle_continuous.csv", index=False)
MAP.to_csv(T / f"{a.inp}_feature_map.csv", index=False)
PH.to_csv(T / f"{a.inp}_phenotypes.csv", index=False)

# 9.8 figure: feature map (rho with ellipticity z at 4 mm, women and men)
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import font_manager  # noqa: E402
for f_ in sorted((ROOT / "assets/fonts").glob("Roboto-*.ttf")):
    font_manager.fontManager.addfont(str(f_))
plt.rcParams["font.family"] = "Roboto"
M = MAP.dropna(subset=["rho_ellip_4mm_all"])
fig, axs = plt.subplots(1, 2, figsize=(12, .32 * len(M) + 1.5), sharey=True)
for ax, o, ttl in zip(axs, OUT, ("2 mm below the annulus", "4 mm below the annulus")):
    y = np.arange(len(M))
    for s, dy in (("F", .12), ("M", -.12)):
        ax.scatter(M[f"rho_{o}_{s}"], y + dy, color=SEXC[s], s=36, label={"F": "Women", "M": "Men"}[s], zorder=3)
    ax.axvline(0, color="#5A6B77", lw=1)
    ax.set_yticks(y, M.label, fontsize=10)
    ax.set_title(ttl, loc="left", fontweight="bold", color="#102D40")
    ax.set_xlabel("Spearman ρ with LVOT ellipticity (both as z for sex, age, BMI)")
    ax.set_xlim(-.8, .8); ax.grid(axis="x", alpha=.25)
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)
axs[0].legend(frameon=False, loc="lower right")
fig.suptitle("Which cardiac features go with a more elliptical LVOT?", x=.01, ha="left", fontweight="bold", color="#102D40")
fig.tight_layout()
fp = ROOT / f"outputs/figures/analysis/fig_{a.inp}_feature_map.png"
fig.savefig(fp, dpi=200)
log_run(__file__, inp=a.inp, outputs=[f"outputs/tables/{a.inp}_feature_map.csv", f"outputs/tables/{a.inp}_phenotypes.csv",
                                      f"outputs/tables/{a.inp}_angle_continuous.csv",
                                      str(fp.relative_to(ROOT))], n=len(D))
