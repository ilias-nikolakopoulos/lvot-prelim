"""STEP 10.6 — Check our strain-for-EF model (STEP 10) against Stokke et al.'s published Equation 1.

Stokke et al., JACC 2017;70:942-954 (doi:10.1016/j.jacc.2017.06.046), validated in 100 patients (EF 16-72 %, r = 0.95).
Their equation is in src/strain.py (stokke_ef). Three checks:
  a. transcription: our code of Equation 1 reproduces the numbers in the paper (Figure 4A, Figure 5A as stated in the text)
  b. like-for-like: Equation 1 with our assumptions (strain at midwall f = 0.5, incompressible c = 1) vs our model,
     same mesh geometry (d = 2 r, w = T from STEP 10.1): EF at the same strain, and strain needed for EF 60 %
  c. published settings (f = strain_model.stokke_f, c = strain_model.stokke_c): strain needed for EF 60 % per mesh, and the
     ellipticity analysis of STEP 10.4 repeated with it
Outputs: outputs/tables/strain_model_check.csv, outputs/figures/analysis/fig_strain_model_check.png
Usage: python scripts/26_strain_model_check.py
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
from src.strain import ef, geometry, gls_needed, stokke_ef, stokke_gls_needed  # noqa: E402

cfg = load_config()
S_, SEXC = cfg["strain_model"], cfg["analysis"]["colours"]
F0, C0 = S_["stokke_f"], S_["stokke_c"]
EF0, K0 = S_["primary_ef"], S_["primary_ratio"]
T_ = ROOT / "outputs/tables"
rows = []

# a. transcription: paper values (text of Results); d = 51 mm, w = 9 mm = control means (Table 2), d = 54 mm = all-subject mean
for name, paper, ours in [
        ("Fig 4A: EF at GLS = GCS = -20 % (d 51, w 9)", 59, 100 * stokke_ef(51, 9, -.2, -.2, F0, C0)),
        ("Fig 4A: EF at GLS 0, GCS -20 % (d 51, w 9)", 43, 100 * stokke_ef(51, 9, 0, -.2, F0, C0)),
        ("Fig 4A: EF at GCS 0, GLS -20 % (d 51, w 9)", 23, 100 * stokke_ef(51, 9, -.2, 0, F0, C0)),
        ("Fig 5A: strain for EF 50 %, w 10 (d 54)", -16.0, 100 * stokke_gls_needed(54, 10, .5, 1, F0, C0)),
        ("Fig 5A: strain for EF 50 %, w 20 (d 54)", -13.0, 100 * stokke_gls_needed(54, 20, .5, 1, F0, C0))]:
    rows.append({"check": "a transcription", "item": name, "paper": paper, "value": ours})

# b, c: per mesh
d = pd.read_csv(T_ / "strain_for_ef_by_mesh.csv")
g = pd.read_csv(T_ / "lvot_all.csv")[["mesh_id", "lv_cavity_ml", "lv_myo_ml", "lv_length"]]
D = d.merge(g, on="mesh_id")
PRIM = f"gls_ef{round(EF0 * 100)}_k{K0}"
D["ours_ef_at_prim"] = [ef(m.r_mm, m.lv_length, m.T_mm, -.2, -.2) for m in D.itertuples()]
D["stk_mid_ef_at_prim"] = [stokke_ef(2 * m.r_mm, m.T_mm, -.2, -.2, .5, 1) for m in D.itertuples()]
D["stk_mid_gls"] = [stokke_gls_needed(2 * m.r_mm, m.T_mm, EF0, K0, .5, 1) for m in D.itertuples()]
D["stk_pub_gls"] = [stokke_gls_needed(2 * m.r_mm, m.T_mm, EF0, K0, F0, C0) for m in D.itertuples()]
D["ours_gls"] = D[PRIM]


def agree(a, b):
    diff = (b - a) * 100
    return {"mean_diff": diff.mean(), "sd_diff": diff.std(), "min_diff": diff.min(), "max_diff": diff.max(),
            "pearson_r": np.corrcoef(a, b)[0, 1], "spearman_rho": stats.spearmanr(a, b)[0]}


rows.append({"check": "b like-for-like", "item": "EF at GLS = GCS = -20 % (% points): Stokke f 0.5, c 1 minus ours",
             **agree(D.ours_ef_at_prim, D.stk_mid_ef_at_prim)})
rows.append({"check": "b like-for-like", "item": f"GLS for EF {EF0:.0%} (% points): Stokke f 0.5, c 1 minus ours",
             **agree(D.ours_gls, D.stk_mid_gls)})
rows.append({"check": "c published settings", "item": f"GLS for EF {EF0:.0%} (% points): Stokke f {F0}, c {C0} minus ours",
             **agree(D.ours_gls, D.stk_pub_gls)})
for s in ("F", "M"):
    for c in ("ours_gls", "stk_pub_gls"):
        rows.append({"check": "c published settings", "item": f"median GLS for EF {EF0:.0%}, {s}: {c}",
                     "value": 100 * D.loc[D.sex == s, c].median()})


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


# c. ellipticity analysis (as STEP 10.4) with both models
D["abs_ours"], D["abs_stk"] = -D.ours_gls, -D.stk_pub_gls
ZC = ["abs_ours", "abs_stk", "ellip_2mm", "ellip_4mm", "lv_cavity_ml", "lv_myo_ml"]
Z = pd.concat([pd.concat({c: zscore(x, c) for c in ZC}, axis=1) for _, x in D.groupby("sex")]).sort_index()
Z["female"] = (D.sex == "F").astype(float)
for y, lab in (("abs_ours", "ours"), ("abs_stk", f"Stokke f {F0}, c {C0}")):
    for o in ("ellip_2mm", "ellip_4mm"):
        for tag, adj in (("", []), (" + EDV, myo volume", ["lv_cavity_ml", "lv_myo_ml"])):
            x = Z[[y, o, "female"] + adj].dropna()
            b, se, p = ols(x[y].values, [x[o].values, x.female.values] + [x[k].values for k in adj])
            rows.append({"check": "c ellipticity slope", "item": f"{lab}: |GLS| z per ellipticity z, {o}{tag}",
                         "value": b, "lo": b - 1.96 * se, "hi": b + 1.96 * se, "p": p, "n": len(x)})
R = pd.DataFrame(rows)
pd.set_option("display.width", 250, "display.max_columns", 20, "display.max_colwidth", 80)
print(R.round(3).to_string(index=False))
R.to_csv(T_ / "strain_model_check.csv", index=False)

# figure: needed GLS per mesh, ours vs Stokke (like-for-like and published settings)
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import font_manager  # noqa: E402
for f_ in sorted((ROOT / "assets/fonts").glob("Roboto-*.ttf")):
    font_manager.fontManager.addfont(str(f_))
plt.rcParams["font.family"] = "Roboto"
INK, MUTED = "#102D40", "#5A6B77"
fig, axs = plt.subplots(1, 2, figsize=(11, 4.8))
for ax, c, ttl in zip(axs, ("stk_mid_gls", "stk_pub_gls"),
                      ("Same assumptions (midwall, incompressible)", f"Stokke's published settings (f {F0}, c {C0})")):
    for s, lab in (("F", "Women"), ("M", "Men")):
        x = D[D.sex == s]
        ax.scatter(-x.ours_gls * 100, -x[c] * 100, s=8, alpha=.4, color=SEXC[s], lw=0, label=lab)
    lim = np.array([min(-D.ours_gls.max(), -D[c].max()), max(-D.ours_gls.min(), -D[c].min())]) * 100 + [-.3, .3]
    ax.plot(lim, lim, color=MUTED, lw=1, ls="--", label="identity")
    a_ = agree(D.ours_gls, D[c])
    ax.text(.03, .97, f"r = {a_['pearson_r']:.3f}\n|GLS| Stokke − ours: {-a_['mean_diff']:+.2f} % points (SD {a_['sd_diff']:.2f})",
            transform=ax.transAxes, va="top", fontsize=9, color=INK)
    ax.set_xlabel(f"|GLS| for EF {EF0:.0%}, our model (%)")
    ax.set_ylabel(f"|GLS| for EF {EF0:.0%}, Stokke Equation 1 (%)")
    ax.set_title(ttl, loc="left", fontweight="bold", color=INK, fontsize=11)
    ax.grid(alpha=.25)
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)
axs[0].legend(frameon=False, loc="lower right", fontsize=9)
fig.suptitle(f"Strain needed for EF {EF0:.0%} per mesh (GCS = GLS): our model vs Stokke et al. 2017",
             x=.01, ha="left", fontweight="bold", color=INK)
fig.tight_layout()
fp = ROOT / "outputs/figures/analysis/fig_strain_model_check.png"
fig.savefig(fp, dpi=200)
log_run(__file__, outputs=["outputs/tables/strain_model_check.csv", str(fp.relative_to(ROOT))], n=len(D))
