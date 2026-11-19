"""STEP 7d — 'What is a mesh?' explainer figure (one mesh, four panels, house style).

(a) the whole heart: outer surface coloured by region (LV wall, RV wall, four valve caps)
(b) zoom on the surface: it is made of small triangles (faces of the tetrahedra)
(c) a thin slice through the wall: the wall is filled with tetrahedra, not just a shell
(d) one tetrahedron: 4 corner points (nodes) with x, y, z coordinates, plus a region tag
Output: outputs/figures/mesh_explainer.png
Usage: python scripts/16_mesh_explainer.py [mesh_id]
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyvista as pv
from matplotlib import font_manager
from matplotlib.collections import PolyCollection
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src.config import load_config, mesh_folder  # noqa: E402
from src.provenance import log_run  # noqa: E402

for f_ in sorted((ROOT / "assets/fonts").glob("Roboto-*.ttf")):
    font_manager.fontManager.addfont(str(f_))
plt.rcParams["font.family"] = "Roboto"
NAVY, INK = "#102D40", "#1E2A33"
TAGC = {0: "#56B4E9", 1: "#D55E00", 4: "#CC79A7", 5: "#E6C229", 6: "#009E73", 7: "#7B4FA6"}
TAGN = {0: "RV wall", 1: "LV wall", 4: "mitral valve cap", 5: "aortic valve cap", 6: "tricuspid valve cap",
        7: "pulmonary valve cap"}

mid = sys.argv[1] if len(sys.argv) > 1 else "female_age_65_bmi_26"
lab = pd.read_csv(ROOT / load_config()["labels"]).set_index("mesh_id")
g = pv.read(mesh_folder(lab.loc[mid, "folder"]) / f"{mid}.fiber.vtu")
tags = np.asarray(g.cell_data["elemTag"])
conn = g.cells.reshape(-1, 5)[:, 1:]
g.points = np.asarray(g.points) * 1e-3
P = np.asarray(g.points)
n_tet, n_pt = len(tags), g.n_points

# heart frame: z = long axis (apex -> mitral centre), x towards the aortic cap
mvc, aoc = P[np.unique(conn[tags == 4])].mean(0), P[np.unique(conn[tags == 5])].mean(0)
lvp = P[np.unique(conn[tags == 1])]
apex = lvp[np.argmax(np.linalg.norm(lvp - mvc, axis=1))]
z = (mvc - apex) / np.linalg.norm(mvc - apex)
x = aoc - mvc; x -= (x @ z) * z; x /= np.linalg.norm(x)
y = np.cross(z, x)
H = lambda X: (np.atleast_2d(X) - apex) @ np.vstack([x, y, z]).T  # noqa: E731

surf = g.extract_surface(pass_cellid=True)
surf.cell_data["tag"] = tags[np.asarray(surf.cell_data["vtkOriginalCellIds"])]
reg = surf.connectivity("largest").extract_surface()
rtag = np.asarray(reg.cell_data["tag"])
T = reg.faces.reshape(-1, 4)[:, 1:]
V = H(np.asarray(reg.points))

fig = plt.figure(figsize=(16, 9))

# (a) whole heart
ax = fig.add_axes([0.0, 0.08, 0.36, 0.84], projection="3d")
ax.add_collection3d(Poly3DCollection(V[T], facecolors=[TAGC[t] for t in rtag], edgecolor="none", alpha=1))
lo, hi = V.min(0), V.max(0)
ax.set_xlim(lo[0], hi[0]); ax.set_ylim(lo[1], hi[1]); ax.set_zlim(lo[2], hi[2])
ax.set_box_aspect(hi - lo); ax.view_init(elev=18, azim=-55); ax.set_axis_off()
for t in (1, 0, 4, 5, 6, 7):
    ax.scatter([], [], [], s=60, color=TAGC[t], marker="s", label=TAGN[t])
ax.legend(fontsize=10, frameon=False, loc="lower left", bbox_to_anchor=(0.02, -0.02))
fig.text(0.02, 0.915, "1  The whole heart", fontsize=15, fontweight="bold", color=NAVY)
fig.text(0.02, 0.88, f"outer surface, coloured by region tag  ·  ~{n_pt / 1000:.0f}k points in total",
         fontsize=10.5, color=INK)

# (b) zoom on surface triangles (LV free wall patch)
ax = fig.add_axes([0.38, 0.5, 0.28, 0.36])
cen = V[T].mean(1)
pick = (rtag == 1) & (np.abs(cen[:, 1]) < 9) & (np.abs(cen[:, 2] - 45) < 9) & (cen[:, 0] < 0)
tri2 = V[T[pick]][:, :, [1, 2]]
ax.add_collection(PolyCollection(tri2, facecolors=TAGC[1], edgecolors="white", linewidths=.8, alpha=.85))
ax.set_xlim(-8, 8); ax.set_ylim(37, 53); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
for s_ in ax.spines.values():
    s_.set_color("0.7")
ax.plot([-7, -2], [38.3, 38.3], color=INK, lw=2.5)
ax.text(-4.5, 38.8, "5 mm", ha="center", fontsize=10, color=INK, fontweight="bold", bbox=dict(fc="white", ec="none", pad=1.5))
fig.text(0.38, 0.915, "2  Zoom on the surface", fontsize=15, fontweight="bold", color=NAVY)
fig.text(0.38, 0.88, "made of small triangles, about 1 mm across", fontsize=10.5, color=INK)

# (c) thin slice through the wall: tetrahedra cut by a short-axis plane, mid-ventricle
ax = fig.add_axes([0.70, 0.5, 0.28, 0.36])
Hp = H(P)
zc = 0.55 * (H(mvc)[0, 2])
tz = Hp[conn][:, :, 2]
cut = (tz.min(1) < zc) & (tz.max(1) > zc) & np.isin(tags, (0, 1))
quads = []
for tet, t in zip(conn[cut], tags[cut]):
    pts = Hp[tet]; d = pts[:, 2] - zc; xs = []
    for i in range(4):
        for j in range(i + 1, 4):
            if d[i] * d[j] < 0:
                xs.append(pts[i] + d[i] / (d[i] - d[j]) * (pts[j] - pts[i]))
    xs = np.array(xs)[:, :2]
    c_ = xs.mean(0); xs = xs[np.argsort(np.arctan2(*(xs - c_).T[::-1]))]
    quads.append((xs, t))
wx = np.concatenate([q[:, 0] for q, _ in quads]); wy = np.concatenate([q[:, 1] for q, _ in quads])
ax.add_collection(PolyCollection([q for q, _ in quads], facecolors=[TAGC[t] for _, t in quads], edgecolors="white",
                                 linewidths=.15, alpha=.9))
ax.set_xlim(wx.min() - 3, wx.max() + 3); ax.set_ylim(wy.min() - 3, wy.max() + 3); ax.set_aspect("equal")
ax.set_xticks([]); ax.set_yticks([])
for s_ in ax.spines.values():
    s_.set_color("0.7")
ax.text(0.02, 0.03, "LV wall (orange) and RV wall (blue);\nthe empty middle is the blood pool", transform=ax.transAxes,
        fontsize=9, color=INK)
fig.text(0.70, 0.915, "3  A slice through the wall", fontsize=15, fontweight="bold", color=NAVY)
fig.text(0.70, 0.88, f"the wall is filled with ~{n_tet / 1e6:.1f} million tetrahedra", fontsize=10.5, color=INK)

# (d) one tetrahedron
ax = fig.add_axes([0.37, 0.02, 0.29, 0.38], projection="3d")
cand = np.flatnonzero(tags == 1)[::200]  # pick a well-shaped tetrahedron: volume / (mean edge)^3 closest to regular
def quality(tet):
    q = P[tet]; v = abs(np.linalg.det(q[1:] - q[0])) / 6
    e = np.mean([np.linalg.norm(q[i] - q[j]) for i in range(4) for j in range(i + 1, 4)])
    return v / e ** 3
k = cand[np.argmax([quality(conn[c]) for c in cand])]
Q = Hp[conn[k]]; Q = Q - Q.mean(0)
faces = [[Q[a], Q[b], Q[c]] for a, b, c in ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3))]
ax.add_collection3d(Poly3DCollection(faces, facecolor=TAGC[1], edgecolor=NAVY, lw=1.5, alpha=.35))
ax.scatter(*Q.T, s=70, color=NAVY, depthshade=False)
for i, q in enumerate(Q):
    ax.text(*(q * 1.85), f"node {i + 1}", fontsize=10.5, color=NAVY, ha="center", va="center", fontweight="bold")
r = np.abs(Q).max() * 1.9
ax.set_xlim(-r, r); ax.set_ylim(-r, r); ax.set_zlim(-r, r); ax.set_box_aspect((1, 1, 1))
# choose the view where the four corners are most spread out on screen (no node hidden behind a face)
from mpl_toolkits.mplot3d import proj3d  # noqa: E402
best = (-1, 20, 35)
for el in (10, 20, 30):
    for az in range(0, 360, 15):
        ax.view_init(elev=el, azim=az)
        xy2 = np.array([proj3d.proj_transform(*q, ax.get_proj())[:2] for q in Q])
        sep = min(np.linalg.norm(xy2[i] - xy2[j]) for i in range(4) for j in range(i + 1, 4))
        hull_inside = any(np.all(np.cross(np.roll(np.delete(xy2, i, 0), -1, 0) - np.delete(xy2, i, 0),
                                          xy2[i] - np.delete(xy2, i, 0)) > 0) for i in range(4))
        if not hull_inside and sep > best[0]:
            best = (sep, el, az)
ax.view_init(elev=best[1], azim=best[2])
ax.set_axis_off()
fig.text(0.38, 0.44, "4  One tetrahedron", fontsize=15, fontweight="bold", color=NAVY)
fig.text(0.38, 0.405, "the basic building block", fontsize=10.5, color=INK)
nodes = P[conn[k]]
txt = ["Each tetrahedron stores:",
       "• 4 corner points (nodes), each with x, y, z",
       f"   e.g. node 1 = ({nodes[0, 0]:.1f}, {nodes[0, 1]:.1f}, {nodes[0, 2]:.1f}) mm",
       "• a region tag: LV wall, RV wall or a valve cap",
       "• fibre direction (not used here)",
       "",
       "Neighbouring tetrahedra share nodes, so",
       "together they form one continuous solid heart.",
       "",
       "In this project, every measurement is made",
       "directly on these points and triangles."]
for i, t_ in enumerate(txt):
    fig.text(0.665, 0.40 - i * 0.032, t_, fontsize=11.5 if i else 12.5, color=NAVY if i in (0, 6, 9) else INK,
             fontweight="bold" if i == 0 else "normal")

fig.suptitle("What is a mesh?", fontsize=22, fontweight="bold", color=NAVY, y=0.998)
out = ROOT / "outputs/figures/mesh_explainer.png"
fig.savefig(out, dpi=150)
print("wrote", out, n_tet, n_pt)
log_run(__file__, mesh=mid, output=str(out.relative_to(ROOT)), n_tet=int(n_tet), n_points=int(n_pt))
