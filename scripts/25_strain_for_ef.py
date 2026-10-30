"""STEP 10 — Strain each end-diastolic geometry needs to reach a given EF, and how it relates to LVOT ellipticity.

Model estimate, NOT measured strain: the meshes are end-diastolic only. Same idea as Stokke et al. (JACC 2017,
doi:10.1016/j.jacc.2017.06.046): EF depends on GLS, GCS, wall thickness and short-axis diameter, so at a fixed EF the
strain needed depends on geometry. Our own derivation (midwall strain, incompressible wall; checked against Stokke in STEP 10.6):
  Model functions in src/strain.py.
  10.1 LV = half prolate ellipsoid (cut at the base). Cavity: long semi-axis L = lv_length, short semi-axis r from
       EDV = 2/3 pi r^2 L. Wall thickness T (uniform) from myocardial volume, outer semi-axes r+T, L+T.
  10.2 Systole: midwall semi-axes shorten by GCS (r + T/2) and GLS (L + T/2); the wall thickens so myocardial volume is
       unchanged (incompressible); EF = 1 - ESV / EDV.
  10.3 For each mesh, solve the GLS (with GCS = ratio x GLS) that gives the target EF, for each EF x ratio scenario.
  10.4 Needed strain magnitude (|GLS|) as z within sex for age and BMI (as STEP 9.6) vs LVOT ellipticity z:
       Spearman rho (women, men, both) and OLS slope (+ female), then + EDV z and myocardial volume z.
  10.5 figure: needed |GLS| vs ellipticity at 2 mm (primary scenario), and the slope in every scenario.
Caveats: global model (says nothing about basal / septal strain near the LVOT); LVOT shape is not in the model, so any link
comes through LV size and wall thickness; representative (cohort-average) meshes.
Outputs: outputs/tables/strain_for_ef_by_mesh.csv, strain_for_ef_vs_ellipticity.csv, outputs/figures/analysis/fig_strain_for_ef.png
Usage: python scripts/25_strain_for_ef.py
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
from src.strain import ef, geometry, gls_needed  # noqa: E402

cfg = load_config()
S_, SEXC = cfg["strain_model"], cfg["analysis"]["colours"]
T_ = ROOT / "outputs/tables"
d = pd.read_csv(T_ / "lvot_all.csv")
A = pd.read_csv(T_ / "area_methods_by_mesh.csv")
A = A[A.level.str.startswith("LVOT")].assign(e=lambda x: x.ct_long / x.ct_short)
E = A.pivot(index="mesh_id", columns="level", values="e")
E.columns = ["ellip_2mm" if "2 mm" in c else "ellip_4mm" for c in E.columns]
D = d[["mesh_id", "sex", "age", "bmi", "lv_cavity_ml", "lv_myo_ml", "lv_length"]].merge(E, left_on="mesh_id", right_index=True)
print(f"{len(D)} meshes ({(D.sex == 'F').sum()} F / {(D.sex == 'M').sum()} M)")

# model check on the median mesh: EF at GLS = GCS = -20 %, and EF change per 1 % point of each strain
med = D[["lv_cavity_ml", "lv_myo_ml", "lv_length"]].median()
r0, T0 = geometry(med.lv_cavity_ml * 1000, med.lv_myo_ml * 1000, med.lv_length)
e0 = ef(r0, med.lv_length, T0, -.2, -.2)
dl = ef(r0, med.lv_length, T0, -.21, -.2) - e0
dc = ef(r0, med.lv_length, T0, -.2, -.21) - e0
print(f"median mesh: r {r0:.1f} mm, L {med.lv_length:.1f} mm, T {T0:.1f} mm; EF at GLS = GCS = -20 %: {e0:.3f}; "
      f"EF per 1 % point GLS {dl * 100:.2f}, GCS {dc * 100:.2f} (ratio {dc / dl:.2f})")

rows = []
for m in D.itertuples():
    r, T = geometry(m.lv_cavity_ml * 1000, m.lv_myo_ml * 1000, m.lv_length)
    row = {"mesh_id": m.mesh_id, "r_mm": r, "T_mm": T}
    for f in S_["ef"]:
        for k in S_["gcs_gls_ratio"]:
            row[f"gls_ef{round(f * 100)}_k{k}"] = gls_needed(r, m.lv_length, T, f, k)
    rows.append(row)
G = pd.DataFrame(rows)
D = D.merge(G, on="mesh_id")
SC = [c for c in G if c.startswith("gls_")]
PRIM = f"gls_ef{round(S_['primary_ef'] * 100)}_k{S_['primary_ratio']}"
print(D[["r_mm", "T_mm", PRIM]].describe().round(3).to_string())


def zscore(df, c):  # as scripts/23 (STEP 9.6)
    x = df[["age", "bmi", c]].dropna()
    u, v = (x.age - 65) / 10, (x.bmi - 27) / 5
    cols = [np.ones(len(x)), u, v]
    if len(x) >= 200:
        cols += [u * u, v * v, u * v, u ** 3, v ** 3, u * u * v, u * v * v]
    X = np.column_stack(cols)
    r = x[c] - X @ np.linalg.lstsq(X, x[c], rcond=None)[0]
    return r / r.std()


def ols(y, X):
    X = np.column_stack([np.ones(len(y))] + X)
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    res = y - X @ b
    s2 = res @ res / (len(y) - X.shape[1])
    se = np.sqrt(np.diag(s2 * np.linalg.inv(X.T @ X)))
    return b[1], se[1], 2 * stats.t.sf(abs(b[1] / se[1]), len(y) - X.shape[1])


# STEP 10.4 — needed strain magnitude vs ellipticity, all as z within sex
D[[f"abs_{c}" for c in SC]] = -D[SC]
ZC = [f"abs_{c}" for c in SC] + ["ellip_2mm", "ellip_4mm", "lv_cavity_ml", "lv_myo_ml"]
Z = pd.concat([pd.concat({c: zscore(g, c) for c in ZC}, axis=1) for _, g in D.groupby("sex")]).sort_index()
Z["sex"], Z["female"] = D.sex, (D.sex == "F").astype(float)
out = []
for c in SC:
    for o in ("ellip_2mm", "ellip_4mm"):
        r = {"scenario": c, "ef": float(c.split("_ef")[1].split("_")[0]) / 100, "gcs_gls_ratio": float(c.split("_k")[1]),
             "ellipticity": o, "median_gls_F": D.loc[D.sex == "F", c].median(), "median_gls_M": D.loc[D.sex == "M", c].median()}
        for s, sub in (("F", Z[Z.sex == "F"]), ("M", Z[Z.sex == "M"]), ("all", Z)):
            x = sub[[f"abs_{c}", o]].dropna()
            r[f"rho_{s}"], r[f"rho_p_{s}"] = stats.spearmanr(x[f"abs_{c}"], x[o])
        for tag, adj in (("", []), ("_adj", ["lv_cavity_ml", "lv_myo_ml"])):
            x = Z[[f"abs_{c}", o, "female"] + adj].dropna()
            b, se, p = ols(x[f"abs_{c}"].values, [x[o].values, x.female.values] + [x[k].values for k in adj])
            r[f"slope{tag}"], r[f"lo{tag}"], r[f"hi{tag}"], r[f"p{tag}"] = b, b - 1.96 * se, b + 1.96 * se, p
        out.append(r)
R = pd.DataFrame(out)
pd.set_option("display.width", 250, "display.max_columns", 30)
print(R.drop(columns=["scenario"]).round(3).to_string(index=False))
# raw scale for the primary scenario: needed |GLS| (% points) per 0.1 ellipticity, within sex
for s in ("F", "M"):
    x = D[D.sex == s].dropna(subset=["ellip_2mm"])
    print(s, f"needed |GLS| per +0.1 ellipticity at 2 mm: {-np.polyfit(x.ellip_2mm, x[PRIM], 1)[0] * 10:.2f} % points")
D[["mesh_id", "sex", "age", "bmi", "ellip_2mm", "ellip_4mm", "r_mm", "T_mm"] + SC].to_csv(T_ / "strain_for_ef_by_mesh.csv", index=False)
R.to_csv(T_ / "strain_for_ef_vs_ellipticity.csv", index=False)

# STEP 10.5 — figure
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import font_manager  # noqa: E402
for f_ in sorted((ROOT / "assets/fonts").glob("Roboto-*.ttf")):
    font_manager.fontManager.addfont(str(f_))
plt.rcParams["font.family"] = "Roboto"
INK, MUTED = "#102D40", "#5A6B77"
fig, axs = plt.subplots(1, 2, figsize=(12, 5), gridspec_kw={"width_ratios": [1.1, 1]})
ax = axs[0]
for s, lab in (("F", "Women"), ("M", "Men")):
    x = D[D.sex == s].dropna(subset=["ellip_2mm"])
    ax.scatter(x.ellip_2mm, -x[PRIM] * 100, s=8, alpha=.35, color=SEXC[s], lw=0, label=lab)
    b = np.polyfit(x.ellip_2mm, -x[PRIM] * 100, 1)
    xx = np.linspace(*x.ellip_2mm.quantile([.02, .98]), 2)
    ax.plot(xx, np.polyval(b, xx), color=SEXC[s], lw=2)
ax.set_xlabel("LVOT ellipticity 2 mm below the annulus (major / minor)")
ax.set_ylabel("|GLS| needed (%)")
ax.set_title(f"Strain needed for EF {S_['primary_ef']:.0%} (GCS = {S_['primary_ratio']} × GLS)", loc="left", fontweight="bold", color=INK)
ax.legend(frameon=False)
ax = axs[1]
P = R[R.ellipticity == "ellip_2mm"].reset_index(drop=True)
y = np.arange(len(P))[::-1]
for tag, dy, col, lab in (("", .15, INK, "sex only"), ("_adj", -.15, MUTED, "+ EDV, myocardial volume")):
    ax.errorbar(P[f"slope{tag}"], y + dy, xerr=[P[f"slope{tag}"] - P[f"lo{tag}"], P[f"hi{tag}"] - P[f"slope{tag}"]],
                fmt="o", color=col, ms=5, capsize=0, label=lab)
ax.axvline(0, color=MUTED, lw=1)
ax.set_yticks(y, [f"EF {e:.0%}, GCS = {k} × GLS" for e, k in zip(P.ef, P.gcs_gls_ratio)], fontsize=9)
ax.set_xlabel("SD change in needed |GLS| per SD of ellipticity (z for sex, age, BMI)")
ax.set_title("Same link in every scenario?", loc="left", fontweight="bold", color=INK)
ax.legend(frameon=False, loc="lower right", fontsize=9)
for a_ in axs:
    a_.grid(alpha=.25)
    for s_ in ("top", "right"):
        a_.spines[s_].set_visible(False)
fig.suptitle("Model estimate from end-diastolic shape, not measured strain", x=.01, ha="left", fontsize=10, color=MUTED)
fig.tight_layout()
fp = ROOT / "outputs/figures/analysis/fig_strain_for_ef.png"
fig.savefig(fp, dpi=200)
log_run(__file__, outputs=["outputs/tables/strain_for_ef_by_mesh.csv", "outputs/tables/strain_for_ef_vs_ellipticity.csv",
                           str(fp.relative_to(ROOT))], n=len(D))
