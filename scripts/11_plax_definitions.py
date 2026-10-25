"""STEP 6.9 (deferred analysis) — How the echo PLAX plane is defined, and how much alternatives matter.

Anchors (all from the mesh):
  aortic axis    line through the aortic rim centre, perpendicular to the rim plane (STEP 3.2)
  mitral centre  centre of the ellipse fitted to the mitral rim nodes (LV nodes touching the mitral plug),
                 projected onto the rim's best-fit plane (STEP 3.1-3.4)
  LV apex        LV endocardial point farthest from the mitral centre (STEP 3.6)
PLAX plane definitions compared:
  A  aortic axis + mitral centre      (used in STEP 3.7)
  B  aortic axis + LV apex
  C  aortic centre + mitral centre + LV apex (3-point plane)
For each, the PLAX line in the LVOT section (2 mm) = where the plane cuts the section; D_PLAX = edge-to-edge
length through the outline's area centroid; echo area = pi/4 * D_PLAX^2 vs the outline area.
Outputs: outputs/tables/plax_definitions.csv (per mesh), outputs/figures/analysis/fig_plax_definitions.png
Usage: python scripts/11_plax_definitions.py [--demo female_age_65_bmi_26] [--ids ...]
       (default meshes: the 8-mesh subset + the 30 pilot pairs, seed from config)
"""
import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import font_manager
from matplotlib.patches import Ellipse, Polygon
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src.config import load_config  # noqa: E402
from src.lvot import chord_through, fit_ellipse, measure_mesh, moment_ellipse, polygon_area_perim  # noqa: E402
from src.provenance import log_run  # noqa: E402

for f_ in sorted((ROOT / "assets/fonts").glob("Roboto-*.ttf")):
    font_manager.fontManager.addfont(str(f_))
plt.rcParams["font.family"] = "Roboto"
NAME_COL = "#102D40"
DEF = {"A": ("A  aortic axis + mitral centre (used)", "#B94F00"),
       "B": ("B  aortic axis + LV apex", "#0072B2"),
       "C": ("C  aortic centre + mitral centre + apex", "#009E73")}

ap = argparse.ArgumentParser()
ap.add_argument("--demo", default="female_age_65_bmi_26")
ap.add_argument("--ids", nargs="+")
a = ap.parse_args()
cfg = load_config()
labels = pd.read_csv(ROOT / cfg["labels"]).set_index("mesh_id")
if a.ids:
    ids = a.ids
else:
    lab = labels.reset_index()
    pr = lab[lab.sex == "F"].merge(lab[lab.sex == "M"], on=["age", "bmi"], suffixes=("_F", "_M")).sample(30, random_state=cfg["seed"])
    ids = list(dict.fromkeys(cfg["subset"] + list(pr.mesh_id_F) + list(pr.mesh_id_M)))


def unit(v):
    return v / np.linalg.norm(v)


def line_angle(a_, b_):
    """Angle between two undirected 2D lines, 0-90 deg."""
    return abs(np.degrees(np.angle(np.exp(2j * (np.arctan2(a_[1], a_[0]) - np.arctan2(b_[1], b_[0])))))) / 2


def analyse(mid):
    r = measure_mesh(labels.loc[mid, "folder"], keep_geometry=True)
    G = r["_geom"]
    ao, mv = G["rings"]["ao"], G["rings"]["mv"]
    n, e1, e2, o = ao["normal"], ao["e1"], ao["e2"], G["origin"]
    proj = lambda v: unit(np.array([v @ e1, v @ e2]))  # noqa: E731
    pn_c = unit(np.cross(mv["centre"] - ao["centre"], G["apex"] - ao["centre"]))
    U = {"A": proj(mv["centre"] - o), "B": proj(G["apex"] - o), "C": proj(np.cross(pn_c, n))}
    s = next(x for x in G["sections"] if x["depth"] == 2.0)
    out = {"mesh_id": mid, "sex": r["sex"], "age": r["age"], "bmi": r["bmi"]}
    geo = {"r": r, "G": G, "U": U, "pn_c": pn_c}
    if not s.get("kept"):
        return out, geo
    xy = s["xy"]
    el = moment_ellipse(xy)
    c = np.array([el["cx"], el["cy"]])
    A_ = polygon_area_perim(xy)[0]
    minor = np.array([np.cos(el["theta"] + np.pi / 2), np.sin(el["theta"] + np.pi / 2)])
    for k, u in U.items():
        D = chord_through(xy, c, u)
        out[f"D_{k}"] = D
        out[f"ang_{k}_vs_A"] = line_angle(u, U["A"])
        out[f"ang_{k}_vs_minor"] = line_angle(u, minor)
        out[f"err_{k}"] = 100 * (np.pi / 4 * D * D - A_) / A_
    geo.update(xy=xy, c=c, minor=minor, el=el, A=A_)
    return out, geo


rows, demo = [], None
for i, mid in enumerate(ids, 1):
    row, geo = analyse(mid)
    rows.append(row)
    if mid == a.demo:
        demo = geo
    print(f"[{i}/{len(ids)}] {mid}", flush=True)
if demo is None:
    demo = analyse(a.demo)[1]
R = pd.DataFrame(rows)
R.round(3).to_csv(ROOT / "outputs/tables/plax_definitions.csv", index=False)

# ---------------- figure ----------------
G, U, r = demo["G"], demo["U"], demo["r"]
ao, mv = G["rings"]["ao"], G["rings"]["mv"]
n, e1, e2, o = ao["normal"], ao["e1"], ao["e2"], G["origin"]
apex = G["apex"]
fig = plt.figure(figsize=(17, 14))

# (a) 3D anchors and planes, in a heart frame (long axis up)
ax = fig.add_subplot(2, 2, 1, projection="3d")
zc = unit(mv["centre"] - apex)
xc = ao["centre"] - mv["centre"]; xc = unit(xc - (xc @ zc) * zc)
yc = np.cross(zc, xc)
H = lambda X: np.column_stack([(np.atleast_2d(X) - mv["centre"]) @ v for v in (xc, yc, zc)])  # noqa: E731
E = np.asarray(G["endo"].points)
pick = np.random.default_rng(0).choice(len(E), 6000, replace=False)
ax.scatter(*H(E[pick]).T, s=.5, color="0.75", alpha=.4, depthshade=False)
ax.scatter(*H(G["ao_ring_nodes"]).T, s=4, color="k", depthshade=False, label="aortic rim")
fe = fit_ellipse(mv["xy"])
mc0 = mv["centre"] - fe["cx"] * mv["e1"] - fe["cy"] * mv["e2"]  # node mean (origin of mv['xy'])
MV3 = mc0 + mv["xy"][:, :1] * mv["e1"] + mv["xy"][:, 1:] * mv["e2"]
ax.scatter(*H(MV3).T, s=4, color="#CC79A7", depthshade=False, label="mitral rim")
t = np.linspace(0, 2 * np.pi, 200)
ex = fe["dmax"] / 2 * np.cos(t); ey = fe["dmin"] / 2 * np.sin(t)
rx = ex * np.cos(fe["theta"]) - ey * np.sin(fe["theta"]); ry = ex * np.sin(fe["theta"]) + ey * np.cos(fe["theta"])
ax.plot(*H(mc0 + np.outer(rx + fe["cx"], mv["e1"]) + np.outer(ry + fe["cy"], mv["e2"])).T, color="#CC79A7", lw=1.5)
ax.scatter(*H(mv["centre"]).T, s=160, marker="*", color="#CC79A7", edgecolor="k", depthshade=False, label="mitral centre")
ax.scatter(*H(apex).T, s=80, marker="^", color="0.2", depthshade=False, label="LV apex")
axl = H(np.vstack([ao["centre"] - 15 * n, ao["centre"] + 30 * n]))
ax.plot(*axl.T, color="k", lw=2, label="aortic axis")
for k, (lab_, col) in DEF.items():
    if k == "C":
        pn = demo["pn_c"]
    else:
        tgt = mv["centre"] if k == "A" else apex
        pn = unit(np.cross(n, tgt - o))
    b1 = zc - (zc @ pn) * pn; b1 = unit(b1); b2 = np.cross(pn, b1)
    cen = ao["centre"] - 30 * zc
    Q = np.array([cen + s1 * 50 * b1 + s2 * 32 * b2 for s1, s2 in ((-1, -1), (1, -1), (1, 1), (-1, 1))])
    ax.add_collection3d(Poly3DCollection([H(Q)], facecolor=col, alpha=.05, edgecolor=col, lw=1.8))
    ax.plot([], [], color=col, lw=6, alpha=.4, label=f"plane {lab_}")
lim = H(E)
ax.set_xlim(lim[:, 0].min(), lim[:, 0].max()); ax.set_ylim(lim[:, 1].min(), lim[:, 1].max()); ax.set_zlim(lim[:, 2].min(), lim[:, 2].max() + 10)
ax.set_box_aspect(np.ptp(lim, axis=0) + [0, 0, 10]); ax.view_init(elev=22, azim=-25); ax.set_axis_off()
ax.legend(fontsize=8.5, loc="lower left", bbox_to_anchor=(-0.05, 0.0), frameon=False)
ax.set_title("(a) Anchors and the three candidate PLAX planes", fontsize=13, fontweight="bold", color=NAME_COL, loc="left")

# (b) mitral centre definition
ax = fig.add_subplot(2, 2, 2)
ax.scatter(*mv["xy"].T, s=6, color="#CC79A7", label="mitral rim nodes (LV nodes touching the mitral plug), projected on their plane")
ax.add_patch(Ellipse((fe["cx"], fe["cy"]), fe["dmax"], fe["dmin"], angle=np.degrees(fe["theta"]), fill=False,
                     ec="#B94F00", lw=2, ls="--", label=f"fitted ellipse ({fe['dmax']:.1f} × {fe['dmin']:.1f} mm)"))
ax.scatter(fe["cx"], fe["cy"], s=260, marker="*", color="#CC79A7", edgecolor="k", zorder=4, label="mitral centre = ellipse centre")
ac = ao["centre"] - mc0
ax.scatter(ac @ mv["e1"], ac @ mv["e2"], s=80, marker="o", color="k", zorder=4, label="aortic centre (projected)")
ax.set_aspect("equal"); ax.set_xlabel("mm, in the mitral rim plane"); ax.grid(alpha=.2)
ax.legend(fontsize=8.5, loc="lower left", framealpha=.9)
ax.set_title("(b) Mitral centre: centre of the ellipse fitted to the mitral rim", fontsize=13, fontweight="bold",
             color=NAME_COL, loc="left")

# (c) PLAX lines on the LVOT 2 mm outline
ax = fig.add_subplot(2, 2, 3)
xy, c, minor, A_ = demo["xy"], demo["c"], demo["minor"], demo["A"]
ax.add_patch(Polygon(xy, closed=True, fc="#0072B2", alpha=.15, ec="#0072B2", lw=2, label=f"LVOT 2 mm outline ({A_:.0f} mm²)"))
ax.plot(*np.vstack([c - minor * 9, c + minor * 9]).T, ls=":", color="0.3", lw=1.5, label="true minor axis")
from matplotlib.path import Path as MPath  # noqa: E402
tt = np.linspace(-30, 30, 12001)
for k, (lab_, col) in DEF.items():
    u = U[k]
    ins = MPath(xy).contains_points(c + tt[:, None] * u)
    p, q = c + tt[ins].min() * u, c + tt[ins].max() * u
    D = np.linalg.norm(q - p)
    err = 100 * (np.pi / 4 * D * D - A_) / A_
    ax.plot(*np.vstack([p, q]).T, color=col, lw=2.5, marker="|", ms=12, mew=2,
            label=f"{lab_.split('  ')[0]}: D {D:.1f} mm, circle {err:+.0f}%")
mdir = unit(np.array([(mv["centre"] - o) @ e1, (mv["centre"] - o) @ e2]))
ax.annotate("", xy=c + mdir * 15, xytext=c + mdir * 12, arrowprops=dict(arrowstyle="->", color="tab:red"))
ax.text(*(c + mdir * 15 + np.array([0, 1.2])), "mitral", color="tab:red", ha="center", fontsize=9)
ax.set_xlim(c[0] - 17, c[0] + 17); ax.set_ylim(c[1] - 17, c[1] + 17); ax.set_aspect("equal")
ax.set_xticks([]); ax.set_yticks([])
ax.legend(fontsize=8.5, loc="lower right", framealpha=.9)
ax.set_title("(c) Where each plane cuts the LVOT (2 mm, looking down the aortic axis)", fontsize=13,
             fontweight="bold", color=NAME_COL, loc="left")

# (d) summary across meshes
ax = fig.add_subplot(2, 2, 4)
Rk = R.dropna(subset=["err_A"])
x = np.arange(3)
for j, sx in enumerate("FM"):
    h = Rk[Rk.sex == sx]
    med = [h[f"err_{k}"].median() for k in DEF]
    q1 = [h[f"err_{k}"].quantile(.25) for k in DEF]; q3 = [h[f"err_{k}"].quantile(.75) for k in DEF]
    xs = x + (j - .5) * .3
    ax.bar(xs, med, width=.28, color=[DEF[k][1] for k in DEF], alpha={"F": .35, "M": .75}[sx],
           edgecolor=[DEF[k][1] for k in DEF], label={"F": "Women", "M": "Men"}[sx])
    ax.errorbar(xs, med, yerr=[np.subtract(med, q1), np.subtract(q3, med)], fmt="none", ecolor="0.3", capsize=3, lw=1)
    for xi, m_ in zip(xs, med):
        ax.text(xi, m_ - 1.5, f"{m_:.0f}%", ha="center", va="top", fontsize=9, color=NAME_COL, fontweight="bold")
angA = [f"{Rk[f'ang_{k}_vs_A'].median():.0f}° (max {Rk[f'ang_{k}_vs_A'].max():.0f}°)" for k in DEF]
angM = [f"{Rk[f'ang_{k}_vs_minor'].median():.0f}°" for k in DEF]
ax.set_xticks(x)
ax.set_xticklabels([f"{k}\nfrom A: {'—' if k == 'A' else angA[i]}\nfrom minor axis: {angM[i]}" for i, k in enumerate(DEF)],
                   fontsize=9, color=NAME_COL)
ax.axhline(0, color="0.3", lw=1)
ax.set_ylim(-45, 3); ax.set_ylabel("Echo circle vs true area, % (median, IQR)")
for s_ in ("top", "right"):
    ax.spines[s_].set_visible(False)
from matplotlib.patches import Patch  # noqa: E402
ax.legend(handles=[Patch(fc="0.5", alpha=.35, label="Women"), Patch(fc="0.5", alpha=.75, label="Men")],
          fontsize=9, frameon=False, loc="lower right")
ax.set_title(f"(d) All three definitions give ~−30% ({len(Rk)} meshes, LVOT 2 mm)", fontsize=13, fontweight="bold",
             color=NAME_COL, loc="left")

fig.suptitle("Defining the echo PLAX view in the mesh", fontsize=17, fontweight="bold", y=0.995)
fig.tight_layout(rect=(0, 0, 1, 0.97), h_pad=3)
out = ROOT / "outputs/figures/analysis/fig_plax_definitions.png"
fig.savefig(out, dpi=115)
print("wrote", out)
log_run(__file__, n_meshes=len(R), demo=a.demo, outputs=["outputs/tables/plax_definitions.csv", str(out.relative_to(ROOT))])
