"""STEP 12 — Synthetic aortic stenosis: remodel each UKB mesh as in Chan 2024, then strain needed for EF by LVOT ellipticity.

Model estimate, NOT measured strain (as STEP 10). The strain models (src/strain.py) use only EDV, myocardial volume and LV
length, so a "synthetic AS mesh" here is each UKB mesh with those three scaled; LVOT shape is kept as measured (the Chan
meshes are LV only, so they say nothing about how the LVOT changes in AS). Equivalent to thickening the wall uniformly.
  12.1 remodelling ratios from the Chan meshes (STEP 11, outputs/tables/chan2024_by_mesh.csv): group statistic
       (as_remodelling.statistic) of EDV, myocardial volume and LV length in each of as_remodelling.groups, divided by
       the reference group (controls)
  12.2 scenarios per UKB mesh: normal (as STEP 10), each Chan group's ratios, and hypertrophy only (myocardial volume x
       as_remodelling.hypertrophy_myo_ratio, EDV and length unchanged); r and T as STEP 10.1
  12.3 |GLS| needed for each EF in strain_model.ef (GCS = strain_model.primary_ratio x GLS), our model and Stokke Equation 1
       (published f, c); change from normal per mesh
  12.4 ellipticity analysis as STEP 10.4 in each scenario: z within sex for age and BMI, OLS slope (+ female),
       then + EDV z and myocardial volume z (of the remodelled geometry); raw % points per +0.1 ellipticity
  12.5 heterogeneity: per-mesh myocardial ratio = AS ratio x exp(N(0, s)), s = excess log-SD of AS over control myocardial
       volume in Chan; slope over as_remodelling.heterogeneity_reps draws (Stokke model, EF primary)
  12.6 figure
Caveats: global model (no basal-septal bulge, which is common in AS and sits next to the LVOT); same strain everywhere;
ratios from 36 AS and 35 control hearts of a different cohort; afterload is not in the model (only geometry).
Outputs: outputs/tables/synthetic_as_by_mesh.csv, synthetic_as_summary.csv, outputs/figures/analysis/fig_synthetic_as_strain.png
Usage: python scripts/28_synthetic_as_strain.py
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
from src.strain import geometry, gls_needed, stokke_gls_needed  # noqa: E402

cfg = load_config()
S_, A_, SEXC = cfg["strain_model"], cfg["as_remodelling"], cfg["analysis"]["colours"]
F0, C0, K0, EF0 = S_["stokke_f"], S_["stokke_c"], S_["primary_ratio"], S_["primary_ef"]
T_ = ROOT / "outputs/tables"
GEO = ["edv_ml", "myo_ml", "lv_length"]

# STEP 12.1 — remodelling ratios from the Chan meshes
C = pd.read_csv(T_ / "chan2024_by_mesh.csv")
ref = C[C.group == A_["reference_group"]][GEO].agg(A_["statistic"])
SCEN = {"normal": dict(edv_ml=1.0, myo_ml=1.0, lv_length=1.0)}
for g in A_["groups"]:
    SCEN[g] = (C[C.group == g][GEO].agg(A_["statistic"]) / ref).to_dict()
for k in A_["hypertrophy_myo_ratio"]:
    SCEN[f"hypertrophy x{k}"] = dict(edv_ml=1.0, myo_ml=k, lv_length=1.0)
for n, v in SCEN.items():
    print(f"{n:18s}" + "  ".join(f"{c} x{v[c]:.3f}" for c in GEO))
lm = lambda g: np.log(C.loc[C.group == g, "myo_ml"])  # noqa: E731
S_HET = np.sqrt(max(lm(A_["groups"][0]).var() - lm(A_["reference_group"]).var(), 0))
print(f"excess log-SD of AS myocardial volume over controls: {S_HET:.3f}")

# UKB meshes, ellipticity as STEP 10
d = pd.read_csv(T_ / "lvot_all.csv")
A = pd.read_csv(T_ / "area_methods_by_mesh.csv")
A = A[A.level.str.startswith("LVOT")].assign(e=lambda x: x.ct_long / x.ct_short)
E = A.pivot(index="mesh_id", columns="level", values="e")
E.columns = ["ellip_2mm" if "2 mm" in c else "ellip_4mm" for c in E.columns]
D = d[["mesh_id", "sex", "age", "bmi", "lv_cavity_ml", "lv_myo_ml", "lv_length"]].merge(E, left_on="mesh_id", right_index=True)
D = D.rename(columns={"lv_cavity_ml": "edv_ml", "lv_myo_ml": "myo_ml"}).reset_index(drop=True)
print(f"{len(D)} UKB meshes")


# STEP 12.2-12.3 — remodelled geometry and needed strain
def needed(edv, myo, L, target):
    r, T = geometry(edv * 1000, myo * 1000, L)
    return r, T, gls_needed(r, L, T, target, K0), stokke_gls_needed(2 * r, T, target, K0, F0, C0)


rows = []
for n, v in SCEN.items():
    for m in D.itertuples():
        edv, myo, L = m.edv_ml * v["edv_ml"], m.myo_ml * v["myo_ml"], m.lv_length * v["lv_length"]
        row = {"mesh_id": m.mesh_id, "scenario": n, "edv_ml": edv, "myo_ml": myo, "lv_length": L}
        for f in S_["ef"]:
            row["r_mm"], row["T_mm"], row[f"ours_ef{round(f * 100)}"], row[f"stokke_ef{round(f * 100)}"] = needed(edv, myo, L, f)
        rows.append(row)
B = pd.DataFrame(rows).merge(D[["mesh_id", "sex", "age", "bmi", "ellip_2mm", "ellip_4mm"]], on="mesh_id")
NC = [c for c in B if c.startswith(("ours_ef", "stokke_ef"))]
B = B.merge(B[B.scenario == "normal"][["mesh_id"] + NC].rename(columns={c: f"{c}_normal" for c in NC}), on="mesh_id")
for c in NC:
    B[f"{c}_change"] = (B[c] - B[f"{c}_normal"]) * -100  # change in |GLS|, % points (+ = more shortening needed)


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


def slopes(x, y, o="ellip_2mm"):
    """STEP 10.4 analysis on one scenario: needed |GLS| z vs ellipticity z (unadjusted, + EDV and myocardial volume)."""
    x = x.assign(abs_=-x[y])
    ZC = ["abs_", o, "edv_ml", "myo_ml"]
    Z = pd.concat([pd.concat({c: zscore(g, c) for c in ZC}, axis=1) for _, g in x.groupby("sex")]).sort_index()
    Z["female"] = (x.sex == "F").astype(float)
    out = {}
    for tag, adj in (("", []), ("_adj", ["edv_ml", "myo_ml"])):
        z = Z[["abs_", o, "female"] + adj].dropna()
        b, se, p = ols(z.abs_.values, [z[o].values, z.female.values] + [z[k].values for k in adj])
        out.update({f"slope{tag}": b, f"lo{tag}": b - 1.96 * se, f"hi{tag}": b + 1.96 * se, f"p{tag}": p})
    for s in ("F", "M"):
        g = x[x.sex == s].dropna(subset=[o])
        out[f"per0.1_{s}"] = np.polyfit(g[o], g.abs_ * 100, 1)[0] / 10
    return out


# STEP 12.4 — summary and ellipticity analysis per scenario
out = []
for n in SCEN:
    x = B[B.scenario == n].reset_index(drop=True)
    for c in NC:
        for o in ("ellip_2mm", "ellip_4mm"):
            r = {"scenario": n, "model": c.split("_ef")[0], "ef": float(c.split("_ef")[1]) / 100, "ellipticity": o,
                 "T_mm_median": x.T_mm.median(), "myo_ml_median": x.myo_ml.median(), "edv_ml_median": x.edv_ml.median()}
            for s in ("F", "M"):
                r[f"gls_median_{s}"] = 100 * x.loc[x.sex == s, c].median()
            r["change_median"], r["change_p5"], r["change_p95"] = x[f"{c}_change"].quantile([.5, .05, .95])
            r["rho_vs_normal"] = stats.spearmanr(x[c], x[f"{c}_normal"])[0]
            r.update(slopes(x, c, o))
            out.append(r)
R = pd.DataFrame(out)

# STEP 12.5 — heterogeneous hypertrophy (Stokke model, primary EF, 2 mm)
rng = np.random.default_rng(A_["seed"])
v = SCEN[A_["groups"][0]]
PC = f"stokke_ef{round(EF0 * 100)}"
het = []
for _ in range(A_["heterogeneity_reps"]):
    k = v["myo_ml"] * np.exp(rng.normal(0, S_HET, len(D)))
    x = D.assign(edv_ml=D.edv_ml * v["edv_ml"], myo_ml=D.myo_ml * k, lv_length=D.lv_length * v["lv_length"])
    x[PC] = [needed(m.edv_ml, m.myo_ml, m.lv_length, EF0)[3] for m in x.itertuples()]
    het.append(slopes(x, PC))
H = pd.DataFrame(het)
hrow = {"scenario": f"{A_['groups'][0]} heterogeneous (log-SD {S_HET:.2f})", "model": "stokke", "ef": EF0,
        "ellipticity": "ellip_2mm", "slope": H.slope.median(), "lo": H.slope.quantile(.025), "hi": H.slope.quantile(.975),
        "slope_adj": H.slope_adj.median(), "lo_adj": H.slope_adj.quantile(.025), "hi_adj": H.slope_adj.quantile(.975),
        "per0.1_F": H["per0.1_F"].median(), "per0.1_M": H["per0.1_M"].median(), "n_reps": len(H)}
R = pd.concat([R, pd.DataFrame([hrow])], ignore_index=True)
pd.set_option("display.width", 250, "display.max_columns", 30)
P = R[(R.ef == EF0) & (R.ellipticity == "ellip_2mm")]
print(P[["scenario", "model", "T_mm_median", "gls_median_F", "gls_median_M", "change_median", "rho_vs_normal", "slope",
         "lo", "hi", "slope_adj", "per0.1_F", "per0.1_M"]].round(3).to_string(index=False))
B[["mesh_id", "scenario", "sex", "age", "bmi", "ellip_2mm", "ellip_4mm", "edv_ml", "myo_ml", "lv_length", "r_mm", "T_mm"]
  + NC + [f"{c}_change" for c in NC]].to_csv(T_ / "synthetic_as_by_mesh.csv", index=False)
R.to_csv(T_ / "synthetic_as_summary.csv", index=False)

# STEP 12.6 — figure
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import font_manager  # noqa: E402
for f_ in sorted((ROOT / "assets/fonts").glob("Roboto-*.ttf")):
    font_manager.fontManager.addfont(str(f_))
plt.rcParams["font.family"] = "Roboto"
INK, MUTED = "#102D40", "#5A6B77"
AS0 = A_["groups"][0]
fig, axs = plt.subplots(1, 2, figsize=(12.5, 5.6), gridspec_kw={"width_ratios": [1.1, 1]})
ax = axs[0]
for n, ls, a in (("normal", "--", .15), (AS0, "-", .35)):
    for s, lab in (("F", "Women"), ("M", "Men")):
        x = B[(B.scenario == n) & (B.sex == s)].dropna(subset=["ellip_2mm"])
        ax.scatter(x.ellip_2mm, -x[PC] * 100, s=7, alpha=a, color=SEXC[s], lw=0)
        b = np.polyfit(x.ellip_2mm, -x[PC] * 100, 1)
        xx = np.linspace(*x.ellip_2mm.quantile([.02, .98]), 2)
        ax.plot(xx, np.polyval(b, xx), color=SEXC[s], lw=2, ls=ls, label=f"{lab}, {'synthetic AS' if n == AS0 else 'as measured'}")
ax.set_xlabel("LVOT ellipticity 2 mm below the annulus (major / minor)")
ax.set_ylabel("|GLS| needed (%)")
ax.set_title(f"Strain needed for EF {EF0:.0%}, Stokke Equation 1 (GCS = {K0} × GLS)", loc="left", fontweight="bold", color=INK, fontsize=11)
ax.legend(frameon=False, fontsize=8.5, ncol=2, loc="upper left")
ax = axs[1]
P = R[(R.ef == EF0) & (R.ellipticity == "ellip_2mm")].reset_index(drop=True)
labs = list(dict.fromkeys(P.scenario))
y = {n: i for i, n in enumerate(labs[::-1])}
for mdl, dy, col in (("ours", .18, MUTED), ("stokke", -.18, INK)):
    for tag, mk, lab in (("", "o", "sex only"), ("_adj", "s", "+ EDV, myocardial volume")):
        q = P[P.model == mdl]
        yy = np.array([y[n] for n in q.scenario]) + dy + (.08 if tag else -.08)
        ax.errorbar(q[f"slope{tag}"], yy, xerr=[q[f"slope{tag}"] - q[f"lo{tag}"], q[f"hi{tag}"] - q[f"slope{tag}"]], fmt=mk,
                    color=col, mfc="white" if tag else col, ms=4.5, capsize=0,
                    label=f"{'our model' if mdl == 'ours' else 'Stokke'}, {lab}")
ax.axvline(0, color=MUTED, lw=1)
ax.set_yticks(list(y.values()), [n.replace("normal", "as measured (UKB)") for n in y], fontsize=9)
ax.set_xlabel("SD change in needed |GLS| per SD of ellipticity at 2 mm")
ax.set_title("Does the ellipticity link change with AS remodelling?", loc="left", fontweight="bold", color=INK, fontsize=11)
ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(.4, -.13), ncol=2, fontsize=8)
for a_ in axs:
    a_.grid(alpha=.25)
    for s_ in ("top", "right"):
        a_.spines[s_].set_visible(False)
fig.suptitle("Model estimate: UKB geometry remodelled with Chan 2024 AS / control ratios, not measured strain",
             x=.01, ha="left", fontsize=10, color=MUTED)
fig.tight_layout()
fp = ROOT / "outputs/figures/analysis/fig_synthetic_as_strain.png"
fig.savefig(fp, dpi=200)
log_run(__file__, outputs=["outputs/tables/synthetic_as_by_mesh.csv", "outputs/tables/synthetic_as_summary.csv",
                           str(fp.relative_to(ROOT))], n=len(D), ratios=SCEN, het_log_sd=S_HET)
