"""STEP 8.1 — Which LV / valve geometry goes with an elliptical LVOT, and does it explain the sex gap?

Table-only (no raw meshes): shape measures already in outputs/tables/lvot_all.csv (STEPS 1.2, 2.2, 3.5, 3.6, 4.7)
plus sphericity, against LVOT ellipticity (minor-major axis ellipse, ct_long / ct_short; STEP 6.6) at 2 mm and 4 mm.
Matched same-age, same-BMI pairs only (paired flag from scripts/06_area_methods.py).
  1. women vs men: median (IQR) of each shape measure, paired difference
  2. Spearman rho of each measure with ellipticity, within each sex
  3. sex gap in ellipticity: OLS ellip ~ female + age_c + bmi_c + age_c x female, then + one standardised
     shape measure at a time, then + all; % of the sex coefficient explained
Outputs: outputs/tables/shape_by_sex.csv, shape_vs_ellipticity.csv, shape_sexgap_attenuation.csv
Usage: python scripts/21_shape_vs_ellipticity.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src.config import load_config  # noqa: E402
from src.provenance import log_run  # noqa: E402

A_ = load_config()["analysis"]
T = ROOT / "outputs/tables"
SHAPE = {  # column: label
    "ao_axis_vs_lv_long_angle": "Aortic axis vs LV long axis (deg)",
    "ao_mv_angle": "Aortic vs mitral annulus plane (deg)",
    "ao_mv_centre_dist": "Aortic-mitral centre distance (mm)",
    "lvot_centreline_vs_ao_normal": "LVOT centreline tilt vs aortic axis (deg)",
    "mv_ellip_ratio": "Mitral annulus ellipticity",
    "mv_d_area": "Mitral annulus diameter (mm)",
    "ao_d_area": "Aortic annulus diameter (mm)",
    "lv_length": "LV length (mm)",
    "lv_cavity_ml": "LV EDV (mL)",
    "lv_myo_ml": "LV myocardial volume (mL)",
    "sphericity": "LV sphericity (EDV / sphere of LV length)",
}
LEVELS = {"LVOT 2 mm (closed)": "2 mm", "LVOT 4 mm (chord-closed)": "4 mm"}

d = pd.read_csv(T / "lvot_all.csv")
flags = d.qc_flags.astype("string").fillna("")
d.loc[flags.str.contains("ao_size_implausible"), ["ao_d_area", "ao_axis_vs_lv_long_angle", "ao_mv_angle",
                                                   "ao_mv_centre_dist", "lvot_centreline_vs_ao_normal"]] = np.nan
d.loc[flags.str.contains("mv_size_implausible"), ["mv_ellip_ratio", "mv_d_area", "ao_mv_angle", "ao_mv_centre_dist"]] = np.nan
d = d.copy()
d["sphericity"] = d.lv_cavity_ml * 1000 / (np.pi / 6 * d.lv_length ** 3)
A = pd.read_csv(T / "area_methods_by_mesh.csv")
A = A[A.paired & A.level.isin(LEVELS)].assign(ellip=lambda x: x.ct_long / x.ct_short)
E = A.pivot(index="mesh_id", columns="level", values="ellip").rename(columns=LEVELS)
pairs_ids = A.mesh_id.unique()
D = d[d.mesh_id.isin(pairs_ids)].merge(E, left_on="mesh_id", right_index=True)
D = D.copy()
D["female"] = (D.sex == "F").astype(float)
D["age_c"] = D.age - A_["age_centre"]
D["bmi_c"] = D.bmi - A_["bmi_ref"]
print(f"{len(D)} meshes in matched pairs ({(D.sex == 'F').sum()} F / {(D.sex == 'M').sum()} M)")

# 1. women vs men (matched pairs: pair by age + bmi, ranks within cell as in 06_area_methods)
D["k"] = D.groupby(["sex", "age", "bmi"]).cumcount()
P = D[D.sex == "F"].merge(D[D.sex == "M"], on=["age", "bmi", "k"], suffixes=("_F", "_M"))
rows = []
for c, lab in SHAPE.items():
    f, m = D.loc[D.sex == "F", c], D.loc[D.sex == "M", c]
    dif = (P[f"{c}_M"] - P[f"{c}_F"]).dropna()
    q = lambda s: f"{s.median():.2f} ({s.quantile(.25):.2f}-{s.quantile(.75):.2f})"
    rows.append({"measure": lab, "women": q(f), "men": q(m), "paired_M_minus_F_median": round(dif.median(), 3),
                 "n_pairs": len(dif), "wilcoxon_p": stats.wilcoxon(dif).pvalue})
S1 = pd.DataFrame(rows)

# 2. Spearman with ellipticity, within sex
rows = []
for c, lab in SHAPE.items():
    r = {"measure": lab}
    for lv in LEVELS.values():
        for sx in "FM":
            x = D[D.sex == sx][[c, lv]].dropna()
            rho, p = stats.spearmanr(x[c], x[lv])
            r[f"rho_{lv}_{sx}"], r[f"p_{lv}_{sx}"] = round(rho, 3), p
    rows.append(r)
S2 = pd.DataFrame(rows)


# 3. sex-gap attenuation
def ols(y, X):
    X = np.column_stack([np.ones(len(y))] + X)
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    res = y - X @ b
    cov = res @ res / (len(y) - X.shape[1]) * np.linalg.inv(X.T @ X)
    return b, np.sqrt(np.diag(cov))


rows = []
for lv in LEVELS.values():
    base_cols = ["female", "age_c", "bmi_c"]
    for extra in [None] + list(SHAPE) + ["all"]:
        cols = list(SHAPE) if extra == "all" else ([] if extra is None else [extra])
        x = D[[lv] + base_cols + cols].dropna()
        Z = [((x[c] - x[c].mean()) / x[c].std()).values for c in cols]
        b, se = ols(x[lv].values, [x.female.values, x.age_c.values, x.bmi_c.values,
                                   (x.age_c * x.female).values] + Z)
        rows.append({"level": lv, "added": "none" if extra is None else ("all" if extra == "all" else SHAPE[extra]),
                     "n": len(x), "female_coef": b[1], "female_lo": b[1] - 1.96 * se[1], "female_hi": b[1] + 1.96 * se[1],
                     "shape_coef_per_sd": b[5] if len(cols) == 1 else np.nan,
                     "shape_lo": b[5] - 1.96 * se[5] if len(cols) == 1 else np.nan,
                     "shape_hi": b[5] + 1.96 * se[5] if len(cols) == 1 else np.nan})
S3 = pd.DataFrame(rows)
for lv in LEVELS.values():
    m = S3.level == lv
    b0 = S3.loc[m & (S3.added == "none"), "female_coef"].iloc[0]
    S3.loc[m, "pct_gap_explained"] = 100 * (1 - S3.loc[m, "female_coef"] / b0)

pd.set_option("display.width", 220, "display.max_columns", 20, "display.max_colwidth", 45)
print(S1.to_string(index=False))
print(S2[["measure"] + [c for c in S2 if c.startswith("rho")]].to_string(index=False))
print(S3.round(3).to_string(index=False))
outs = {"shape_by_sex.csv": S1, "shape_vs_ellipticity.csv": S2, "shape_sexgap_attenuation.csv": S3}
for k, v in outs.items():
    v.to_csv(T / k, index=False)
log_run(__file__, outputs=[f"outputs/tables/{k}" for k in outs], n=len(D))
