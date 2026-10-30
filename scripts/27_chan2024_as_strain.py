"""STEP 11 — Strain models (STEP 10) in real aortic stenosis hearts: Chan et al. 2024 public data.

Data: Chan et al., Circ Cardiovasc Imaging 2024 (doi:10.1161/CIRCIMAGING.123.016489), FigShare 10.6084/m9.figshare.25428604,
CC BY 4.0. 109 end-diastolic LV meshes from CMR (severe AS, AS 8 months after surgical AVR, hypertension, healthy controls),
plus a per-case sheet with CMR feature-tracking strain (AS, post-AVR, part of HTN; none for controls) and echo EF after AVR.
LV only (short-axis contours): no LVOT or valves, so nothing here is about LVOT shape.
  11.1 download the sheet and VTK meshes if missing (config external.chan2024_as)
  11.2 per mesh: EDV (endocardium capped at its basal rim), myocardial volume (closed endo + epi + basal ring surface),
       LV length (basal rim centre to farthest endocardial point), median wall thickness (endocardial points to epicardium)
  11.3 forward test: EF predicted from measured GLS, GCS and the mesh geometry (our model; Stokke Equation 1 with strain
       layer f in external.chan2024_as.strain_layer_f, c = strain_model.stokke_c) vs measured EF:
       per case for post-AVR (echo EF after AVR), group means for AS and HTN (CMR EF, paper Table 1),
       controls from the paper's group mean strain (Table 2) on each control geometry
  11.4 backward use (as STEP 10): strain needed for EF 60 % (GCS = GLS) by group, vs measured strain;
       per case, needed vs measured |GLS| and |GCS| (Spearman)
Outputs: outputs/tables/chan2024_by_mesh.csv, chan2024_summary.csv, outputs/figures/analysis/fig_chan2024_strain_check.png
Usage: python scripts/27_chan2024_as_strain.py
"""
import sys
import urllib.request
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import pyvista as pv
import vtk
from scipy import stats
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src.config import load_config  # noqa: E402
from src.provenance import log_run  # noqa: E402
from src.strain import ef, geometry, gls_needed, stokke_ef, stokke_gls_needed  # noqa: E402

vtk.vtkObject.GlobalWarningDisplayOff()
cfg = load_config()
X_, S_ = cfg["external"]["chan2024_as"], cfg["strain_model"]
DIR = ROOT / X_["dir"]
T_ = ROOT / "outputs/tables"
GROUPS = {1: "AS", 2: "AS", 3: "AS post-AVR", 4: "Control", 5: "HTN"}  # 1/2 = asymptomatic / symptomatic severe AS
# paper Tables 1-2 (group means; CMR EF and feature-tracking strain); post-AVR not reported there
PAPER = {"AS": {"ef": 73, "edv": 141, "gls": -14.4, "gcs": -17.9}, "Control": {"ef": 69, "edv": 137, "gls": -16.0, "gcs": -20.0},
         "HTN": {"ef": 69, "edv": 146, "gls": -11.0, "gcs": -17.7}}

# STEP 11.1 — download
DIR.mkdir(parents=True, exist_ok=True)
for name, url in X_["files"].items():
    if not (DIR / name).exists():
        urllib.request.urlretrieve(url, DIR / name)
MD = DIR / "VTKmeshesAlignedDofs"
if not MD.exists():
    zipfile.ZipFile(DIR / "VTKmeshesAlignedDofs.zip").extractall(DIR)
LAB = np.loadtxt(MD / "NodeLabels.txt", dtype=int)  # 1 endocardium, 2 epicardium, 0 basal ring


# STEP 11.2 — geometry per mesh
def tri_surface(p, f):
    return pv.PolyData(p, np.hstack([np.full((len(f), 1), 3), f]).ravel()).clean(tolerance=1e-4).triangulate()


def signed_volume(s):
    s = s.compute_normals(auto_orient_normals=True, consistent_normals=True, cell_normals=False)
    t, q = s.faces.reshape(-1, 4)[:, 1:], s.points
    return abs(np.einsum("ij,ij->i", q[t[:, 0]], np.cross(q[t[:, 1]], q[t[:, 2]])).sum() / 6)


def boundary(s):
    return s.extract_feature_edges(boundary_edges=True, feature_edges=False, manifold_edges=False, non_manifold_edges=False)


def capped(s):  # close the open base with a fan to the rim centre
    e = boundary(s)
    c = e.points.mean(0)
    ln = e.lines.reshape(-1, 3)[:, 1:]
    idx = cKDTree(s.points).query(e.points)[1]
    tri = np.vstack([s.faces.reshape(-1, 4)[:, 1:], np.column_stack([idx[ln[:, 0]], idx[ln[:, 1]], np.full(len(ln), s.n_points)])])
    return tri_surface(np.vstack([s.points, c]), tri), c


def measure(fp):
    m = pv.read(fp)
    f = m.faces.reshape(-1, 4)[:, 1:]
    L = LAB[f]
    endo = tri_surface(m.points, f[(L == 1).any(1) & ~(L == 2).any(1)])
    epi = tri_surface(m.points, f[(L == 2).any(1) & ~(L == 1).any(1)])
    shell = tri_surface(m.points, f)
    cav, c = capped(endo)
    assert boundary(cav).n_cells == 0 and boundary(shell).n_cells == 0, fp.name
    wt = cKDTree(epi.points).query(endo.points)[0]  # nearest epicardial node (node spacing ~2 mm, so slightly high)
    return {"edv_ml": signed_volume(cav) / 1000, "myo_ml": signed_volume(shell) / 1000,
            "lv_length": np.linalg.norm(endo.points - c, axis=1).max(), "wt_median": np.median(wt)}


x = pd.read_excel(DIR / "infocases.xls", engine="xlrd")
x = x.rename(columns={"Mesh ID": "mesh", "Subject ID": "subject", "circum Stain": "gcs", "longit Strain": "gls",
                      "Post AVR EF (ECHO) ": "ef_echo_post", "CORRECT 1=AS asym; 2=AS sym; 3=AS postAVR; 4=control; 5=HYT": "code"})
x["group"] = x.code.map(GROUPS)
x = x[x.group.notna()][["mesh", "subject", "group", "code", "gls", "gcs", "ef_echo_post"]]
x[["gls", "gcs"]] /= 100
rows = []
for m in x.itertuples():
    fp = MD / f"Case{int(m.mesh)}.vtk"
    if fp.exists():
        rows.append({"mesh": m.mesh, **measure(fp)})
D = x.merge(pd.DataFrame(rows), on="mesh")
print(f"{len(D)} meshes with a group; {len(x) - len(D)} sheet rows without a mesh")
D["r_mm"], D["T_mm"] = zip(*[geometry(a * 1000, b * 1000, c) for a, b, c in zip(D.edv_ml, D.myo_ml, D.lv_length)])

# STEP 11.3 — forward: EF from strain + geometry (controls: paper group-mean strain)
gl = D.gls.where(D.gls.notna() | (D.group != "Control"), PAPER["Control"]["gls"] / 100)
gc = D.gcs.where(D.gcs.notna() | (D.group != "Control"), PAPER["Control"]["gcs"] / 100)
D["strain_source"] = np.where(D.gls.notna(), "per case", np.where(D.group == "Control", "paper group mean", "none"))
ok = gl.notna() & gc.notna()
D["ef_ours"] = np.where(ok, [ef(r, L, T, a, b) if k else np.nan for r, L, T, a, b, k in zip(D.r_mm, D.lv_length, D.T_mm, gl, gc, ok)], np.nan)
MODELS = ["ef_ours"]
for f in X_["strain_layer_f"]:
    c = f"ef_stokke_f{f}"
    D[c] = np.where(ok, stokke_ef(2 * D.r_mm, D.T_mm, gl, gc, f, S_["stokke_c"]), np.nan)
    MODELS.append(c)
D[MODELS] *= 100

# STEP 11.4 — backward: strain needed for EF 60 %
EF0, K0 = S_["primary_ef"], S_["primary_ratio"]
D["gls_needed_ours"] = [gls_needed(r, L, T, EF0, K0) for r, L, T in zip(D.r_mm, D.lv_length, D.T_mm)]
D["gls_needed_stokke"] = [stokke_gls_needed(2 * r, T, EF0, K0, S_["stokke_f"], S_["stokke_c"]) for r, T in zip(D.r_mm, D.T_mm)]

out = []
for g, s in D.groupby("group"):
    r = {"group": g, "n": len(s), "n_strain": int(s.gls.notna().sum())}
    for c in ["edv_ml", "myo_ml", "lv_length", "wt_median", "T_mm", "r_mm", "gls", "gcs", "gls_needed_ours", "gls_needed_stokke"] + MODELS:
        r[c] = s[c].mean() * (100 if c.startswith("gls") or c == "gcs" else 1)
        r[c + "_sd"] = s[c].std() * (100 if c.startswith("gls") or c == "gcs" else 1)
    r.update({f"paper_{k}": v for k, v in PAPER.get(g, {}).items()})
    out.append(r)
G = pd.DataFrame(out)
pd.set_option("display.width", 250, "display.max_columns", 40)
show = ["group", "n", "n_strain", "edv_ml", "paper_edv", "myo_ml", "lv_length", "wt_median", "T_mm", "gls", "paper_gls", "gcs",
        "paper_gcs", "gls_needed_ours", "gls_needed_stokke"] + MODELS + ["paper_ef"]
print(G[show].round(1).to_string(index=False))

# per case: post-AVR predicted vs echo EF; needed vs measured strain
P = D[(D.group == "AS post-AVR") & D.ef_echo_post.notna() & D.gls.notna()]
acc = []
for c in MODELS:
    dd = P[c] - P.ef_echo_post
    acc.append({"group": "AS post-AVR, per case vs echo EF", "model": c, "n": len(P), "mean_diff": dd.mean(), "sd_diff": dd.std(),
                "r": np.corrcoef(P[c], P.ef_echo_post)[0, 1] if len(P) > 2 else np.nan})
for g in ("AS", "HTN", "Control"):
    s = D[D.group == g]
    for c in MODELS:
        acc.append({"group": f"{g}, group mean vs paper CMR EF", "model": c, "n": int(s[c].notna().sum()),
                    "mean_diff": s[c].mean() - PAPER[g]["ef"]})
A = pd.DataFrame(acc)
print(A.round(3).to_string(index=False))
H = D[D.gls.notna()]
for y in ("gls", "gcs"):
    rho, p = stats.spearmanr(-H.gls_needed_ours, -H[y])
    print(f"needed |GLS| (ours) vs measured |{y.upper()}|, all cases with strain (n={len(H)}): rho {rho:.2f}, p {p:.3f}")
    for g, s in H.groupby("group"):
        rho, p = stats.spearmanr(-s.gls_needed_ours, -s[y])
        print(f"   {g} (n={len(s)}): rho {rho:.2f}, p {p:.3f}")
# paired pre / post AVR
pre = D[D.group == "AS"].set_index("subject")
post = D[D.group == "AS post-AVR"].set_index("subject")
pp = pre.join(post, lsuffix="_pre", rsuffix="_post", how="inner")
for c in ("edv_ml", "myo_ml", "T_mm", "gls_needed_ours", "gls", "gcs"):
    d_ = (pp[f"{c}_post"] - pp[f"{c}_pre"]).dropna() * (100 if c.startswith("gls") or c == "gcs" else 1)
    print(f"post - pre AVR {c}: {d_.mean():+.2f} (SD {d_.std():.2f}, n {len(d_)}, Wilcoxon p {stats.wilcoxon(d_).pvalue:.3f})")

D.to_csv(T_ / "chan2024_by_mesh.csv", index=False)
pd.concat([G, A]).to_csv(T_ / "chan2024_summary.csv", index=False)

# figure
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import font_manager  # noqa: E402
for f_ in sorted((ROOT / "assets/fonts").glob("Roboto-*.ttf")):
    font_manager.fontManager.addfont(str(f_))
plt.rcParams["font.family"] = "Roboto"
INK, MUTED = "#102D40", "#5A6B77"
GC = {"Control": "#009E73", "AS": "#D55E00", "AS post-AVR": "#E69F00", "HTN": "#0072B2"}
fig, axs = plt.subplots(1, 2, figsize=(12, 5))
ax = axs[0]
SF = f"ef_stokke_f{S_['stokke_f']}"
ax.scatter(P.ef_echo_post, P[SF], color=GC["AS post-AVR"], s=30, label="post-AVR, per case (echo EF)", zorder=3)
for g in ("AS", "HTN", "Control"):
    s = D[D.group == g]
    ax.errorbar(PAPER[g]["ef"], s[SF].mean(), yerr=s[SF].std(), fmt="s", ms=8, color=GC[g], capsize=3,
                label=f"{g}, group mean (paper CMR EF)" + (" *" if g == "Control" else ""))
ax.plot([35, 95], [35, 95], ls="--", color=MUTED, lw=1)
ax.set_xlabel("Measured EF (%)")
ax.set_ylabel("EF predicted from measured strain + mesh geometry (%)")
ax.set_title(f"Forward test: Stokke Equation 1 (f {S_['stokke_f']}, c {S_['stokke_c']})", loc="left", fontweight="bold", color=INK, fontsize=11)
ax.legend(frameon=False, fontsize=8, loc="upper left")
ax.text(.99, .02, "* controls: paper group-mean strain", transform=ax.transAxes, ha="right", fontsize=8, color=MUTED)
ax = axs[1]
for g, s in H.groupby("group"):
    ax.scatter(-s.gls_needed_ours * 100, -s.gls * 100, color=GC[g], s=26, label=g)
ax.set_xlabel(f"|GLS| needed for EF {EF0:.0%} from geometry (our model, %)")
ax.set_ylabel("|GLS| measured, CMR feature tracking (%)")
ax.set_title("Backward use: does geometry track measured strain?", loc="left", fontweight="bold", color=INK, fontsize=11)
ax.legend(frameon=False, fontsize=9)
for a_ in axs:
    a_.grid(alpha=.25)
    for s_ in ("top", "right"):
        a_.spines[s_].set_visible(False)
fig.suptitle("Strain models in aortic stenosis hearts (Chan et al. 2024 public data, CC BY 4.0)", x=.01, ha="left", fontweight="bold", color=INK)
fig.tight_layout()
fp = ROOT / "outputs/figures/analysis/fig_chan2024_strain_check.png"
fig.savefig(fp, dpi=200)
log_run(__file__, outputs=["outputs/tables/chan2024_by_mesh.csv", "outputs/tables/chan2024_summary.csv", str(fp.relative_to(ROOT))], n=len(D))
