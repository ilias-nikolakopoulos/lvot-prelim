"""STEP 5.5 — Every step from raw mesh to measured LVOT outline, for one mesh (3 x 3 panels).

1 mesh and tags            2 three closed surfaces       3 LV endocardium + aortic rim plane
4 cutting planes           5 plane x triangles -> segs   6 segments -> closed chain (2 mm) + checks
7 open chain (4 mm)        8 chord closure (4 mm)        9 polygon -> area (shoelace)
Usage: python scripts/08_outline_steps.py [mesh_id] -> outputs/figures/qc/outline_steps_<id>.png
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import LineCollection, PolyCollection
from matplotlib.patches import Circle, Polygon, Wedge
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src.config import load_config  # noqa: E402
from src.lvot import chain_loops, measure_mesh, moment_ellipse, polygon_area_perim  # noqa: E402
from src.provenance import log_run  # noqa: E402

mid = sys.argv[1] if len(sys.argv) > 1 else "female_age_65_bmi_26"
cfg = load_config()
labels = pd.read_csv(ROOT / cfg["labels"]).set_index("mesh_id")
r = measure_mesh(labels.loc[mid, "folder"], keep_geometry=True)
G = r["_geom"]
ao, mv = G["rings"]["ao"], G["rings"]["mv"]
n, e1, e2, org = ao["normal"], ao["e1"], ao["e2"], G["origin"]
endo, cav, reg = G["endo"], G["cavity"], G["regions"]
S = {s["depth"]: s for s in G["sections"]}

TAGC = {0: "#56B4E9", 1: "#D55E00", 4: "#CC79A7", 5: "#E69F00", 6: "#009E73", 7: "#0072B2"}
TAGN = {0: "RV wall", 1: "LV wall", 4: "mitral plug", 5: "aortic plug", 6: "tricuspid plug", 7: "pulmonary plug"}

# heart frame: z = long axis (apex -> mitral centre, base up), x towards the aortic side
apex = G["apex"]
zc = mv["centre"] - apex; zc /= np.linalg.norm(zc)
xc = (ao["centre"] - mv["centre"]); xc -= (xc @ zc) * zc; xc /= np.linalg.norm(xc)
yc = np.cross(zc, xc)
H0 = mv["centre"]


def hf(X):
    D = X - H0
    return np.column_stack([D @ xc, D @ yc, D @ zc])


def tris(poly):
    return poly.faces.reshape(-1, 4)[:, 1:]


fig = plt.figure(figsize=(19, 19))
k = 0


def panel(title, proj=None):
    global k
    k += 1
    ax = fig.add_subplot(3, 3, k, projection=proj)
    ax.set_title(f"{k}  {title}", fontsize=11, loc="left")
    return ax


# ---------- 1 mesh and tags ----------
ax = panel("Mesh: ~2 M tetrahedra, each tagged\n(outer surface shown, coloured by tag)", "3d")
rid = np.asarray(reg.cell_data["RegionId"]); rtag = np.asarray(reg.cell_data["tag"])
outer = np.bincount(rid).argmax()
T = tris(reg)
sel = np.flatnonzero(rid == outer)
VV = hf(np.asarray(reg.points))
C = VV[T[sel]].mean(1)
fc = [TAGC[t] for t in rtag[sel]]
ax.add_collection3d(Poly3DCollection(VV[T[sel]], facecolors=fc, edgecolor="none", alpha=1))
for t in (0, 1, 4, 5, 6, 7):
    ax.scatter([], [], [], s=30, color=TAGC[t], label=TAGN[t])
ax.set_xlim(C[:, 0].min(), C[:, 0].max()); ax.set_ylim(C[:, 1].min(), C[:, 1].max()); ax.set_zlim(C[:, 2].min(), C[:, 2].max())
ax.view_init(elev=28, azim=-60)
ax.set_box_aspect(np.ptp(C, axis=0)); ax.set_axis_off()
ax.legend(fontsize=8, loc="upper left", frameon=False)

# ---------- 2 three closed surfaces (long-axis cut through both valves) ----------
ax = panel("Outer skin of all tets = 3 closed surfaces\n(long-axis cut through aortic and mitral openings)")
pn = np.cross(zc, xc)  # plane containing the long axis and the aortic direction
sl = reg.slice(normal=pn, origin=(ao["centre"] + mv["centre"]) / 2)
P2 = hf(np.asarray(sl.points))[:, [0, 2]]
segs = sl.lines.reshape(-1, 3)[:, 1:]
srid, stag = np.asarray(sl.cell_data["RegionId"]), np.asarray(sl.cell_data["tag"])
regs = [x for x in np.unique(srid) if x != outer]
names = {outer: ("outer surface", "0.2")}
for x in regs:
    names[x] = ("LV cavity (LV wall + mitral/aortic plugs)", "#D55E00") if x == G["lv_region"] else \
        ("RV cavity (RV wall + tricuspid/pulmonary plugs)", "#56B4E9")
for x, (nm, col) in names.items():
    m = (srid == x) & (stag < 4)
    ax.add_collection(LineCollection(P2[segs[m]], colors=col, lw=1.6, label=nm))
    m = (srid == x) & (stag >= 4)
    ax.add_collection(LineCollection(P2[segs[m]], colors=col, lw=3.5, linestyles=":"))
ax.plot([], [], ":", color="0.4", lw=3, label="plug faces (dotted) close each opening")
ax.autoscale(); ax.set_aspect("equal"); ax.set_xlabel("mm"); ax.set_ylabel("mm (long axis, base up)")
ax.legend(fontsize=8, loc="lower left", frameon=False)
ax.text(0.98, 0.02, f"LV cavity volume = {r['lv_cavity_ml']:.0f} mL", transform=ax.transAxes, ha="right", fontsize=9)

# ---------- 3 LV endocardium + aortic rim plane ----------
ax = panel("LV cavity seen from inside, looking up the aortic axis:\ndrop plug faces → LV endocardium; fit a plane to the aortic rim", "3d")
ct = np.asarray(cav.cell_data["tag"])
CT = tris(cav)
V = hf(np.asarray(cav.points))
cen = V[CT].mean(1)
dep3 = (np.asarray(cav.points)[CT].mean(1) - org) @ n  # depth below the aortic rim along its axis
near = (np.linalg.norm(cen - hf(ao["centre"][None])[0], axis=1) < 38) & (dep3 < 6)  # basal part only
for t, a in ((1, .25), (4, .7), (5, .9)):
    m = near & (ct == t)
    ax.add_collection3d(Poly3DCollection(V[CT[m]], facecolor=TAGC[t] if t != 1 else "0.7",
                                         edgecolor="none", alpha=a))
ax.scatter([], [], color="0.7", label="LV endocardium (kept)")
ax.scatter([], [], color=TAGC[5], label="aortic plug faces (dropped)")
ax.scatter([], [], color=TAGC[4], label="mitral plug faces (dropped)")
RN = hf(G["ao_ring_nodes"])
ax.scatter(*RN.T, s=3, color="k", label="aortic rim nodes (LV ∩ aortic plug)")
sq = np.array([[-12, -12], [12, -12], [12, 12], [-12, 12]])
Q = hf(ao["centre"] + sq[:, :1] * e1 + sq[:, 1:] * e2)
ax.add_collection3d(Poly3DCollection([Q], facecolor="k", alpha=.12, edgecolor="k"))
a0 = hf(ao["centre"][None])[0]; nn = (hf((ao["centre"] + 12 * n)[None])[0] - a0)
ax.quiver(*a0, *nn, color="k", arrow_length_ratio=.2)
ax.text(*(a0 + nn * 1.15), "aortic axis\n(into LV)", fontsize=8)
lim = np.vstack([V[CT[near]].reshape(-1, 3)])
ax.set_xlim(lim[:, 0].min(), lim[:, 0].max()); ax.set_ylim(lim[:, 1].min(), lim[:, 1].max()); ax.set_zlim(lim[:, 2].min(), lim[:, 2].max())
nh = hf((ao["centre"] + n)[None])[0] - a0  # camera inside the LV, looking up the aortic axis at the base
ax.set_box_aspect(np.ptp(lim, axis=0)); ax.view_init(elev=np.degrees(np.arcsin(nh[2])), azim=np.degrees(np.arctan2(nh[1], nh[0])))
ax.set_axis_off()
ax.legend(fontsize=8, loc="upper left", frameon=False, markerscale=3)

# ---------- 4 cutting planes ----------
ax = panel("Cutting planes perpendicular to the aortic axis, every 1 mm below\nthe ventricular edge of the aortic rim (side view)")
md = mv["centre"] - org
u = np.array([md @ e1, md @ e2]); u /= np.linalg.norm(u)
u3 = u[0] * e1 + u[1] * e2
sl = cav.slice(normal=np.cross(n, u3), origin=org)
D = np.asarray(sl.points) - org
P4 = np.column_stack([D @ u3, -(D @ n)])
segs = sl.lines.reshape(-1, 3)[:, 1:]
stg = np.asarray(sl.cell_data["tag"])
for t, col, lw, lab in ((1, "0.3", 1.6, "LV endocardium"), (5, TAGC[5], 3, "aortic plug face"), (4, TAGC[4], 3, "mitral plug face")):
    ax.add_collection(LineCollection(P4[segs[stg == t]], colors=col, lw=lw, label=lab))
for d in range(0, 9):
    col = "#0072B2" if d == 2 else ("#009E73" if d == 4 else "0.75")
    ax.axhline(-d, color=col, lw=1.4 if d in (2, 4) else .6, ls="-" if d in (2, 4) else "--")
    ax.text(-19.5, -d + .15, f"{d} mm", fontsize=7, color=col)
ax.text(10, -2 + .3, "2 mm (primary)", color="#0072B2", fontsize=8)
ax.text(10, -4 + .3, "4 mm (secondary)", color="#009E73", fontsize=8)
ax.set_xlim(-20, 24); ax.set_ylim(-14, 6); ax.set_aspect("equal")
ax.set_xlabel("mm, towards mitral →"); ax.set_ylabel("mm along aortic axis (aorta up)")
ax.legend(fontsize=8, loc="lower left", frameon=False)

# ---------- 5 plane x triangles -> segments ----------
ax = panel("Each triangle the plane crosses gives one short\nsegment (zoom on the wall at 2 mm, en face)")
o = org + 2.0 * n
ET = tris(endo); EP = np.asarray(endo.points)
sd = (EP - o) @ n
crosses = (np.sign(sd[ET]).min(1) < 0) & (np.sign(sd[ET]).max(1) > 0)
loop2 = S[2.0]["xy"]
focus = loop2[np.argmax(loop2 @ np.array([1.0, 0.3]))]
inplane = lambda X: np.column_stack([(X - o) @ e1, (X - o) @ e2])
Tq = ET[crosses]
Pq = inplane(EP[Tq.reshape(-1)]).reshape(-1, 3, 2)
keep = np.linalg.norm(Pq.mean(1) - focus, axis=1) < 7
ax.add_collection(PolyCollection(Pq[keep], facecolors="0.9", edgecolors="0.55", lw=.6))
segl = []
for tri in Tq[keep]:
    dv = sd[tri]; pts = []
    for a_, b_ in ((0, 1), (1, 2), (2, 0)):
        if np.sign(dv[a_]) != np.sign(dv[b_]):
            t_ = dv[a_] / (dv[a_] - dv[b_])
            pts.append(EP[tri[a_]] + t_ * (EP[tri[b_]] - EP[tri[a_]]))
    if len(pts) == 2:
        segl.append(inplane(np.array(pts)))
cols = ["#0072B2", "#D55E00"]
ax.add_collection(LineCollection(segl, colors=[cols[i % 2] for i in range(len(segl))], lw=2.6))
allp = np.vstack(segl)
ax.scatter(*allp.T, s=10, color="k", zorder=3)
ax.plot([], [], color="0.55", lw=1, label="triangles crossed by the plane (projected)")
ax.plot([], [], color="#0072B2", lw=2.6, label="one segment per triangle (alternating colours)")
ax.scatter([], [], s=10, color="k", label="crossing points = outline points")
ax.set_xlim(focus[0] - 6, focus[0] + 6); ax.set_ylim(focus[1] - 6, focus[1] + 6); ax.set_aspect("equal")
ax.set_xlabel("mm"); ax.legend(fontsize=8, loc="lower left", framealpha=.9)

# ---------- 6 chain -> closed loop + checks ----------
ax = panel("Join segments end to end → chain. At 2 mm it closes\non itself → closed tube; check it is the LVOT")
xy = loop2
sc = ax.scatter(*xy.T, c=np.arange(len(xy)), cmap="viridis", s=12, zorder=3)
ax.plot(*np.vstack([xy, xy[:1]]).T, color="0.5", lw=.8)
for i in range(0, len(xy), max(1, len(xy) // 12)):
    j = (i + 3) % len(xy)
    ax.annotate("", xy=xy[j], xytext=xy[i], arrowprops=dict(arrowstyle="->", color="0.2", lw=1.2))
ax.scatter(*xy[0], s=120, marker="*", color="k", zorder=4, label="start = end (chain closes)")
ax.plot(0, 0, "k+", ms=12, mew=2, label="aortic axis (must be inside)")
ra = ao["d_area"] / 2
ax.add_patch(Circle((0, 0), ra, fill=False, ls="--", color="0.4", label=f"annulus-size circle (Ø {ao['d_area']:.1f} mm)"))
ax.add_patch(Circle((0, 0), 1.6 * ra, fill=False, ls=":", color="tab:red",
                    label=f"limit 1.6× annulus (Ø {1.6 * ao['d_area']:.1f}) — larger = cavity"))
cb = plt.colorbar(sc, ax=ax, fraction=.04, pad=.01); cb.set_label("order along the chain", fontsize=8)
ax.set_xlim(-17, 17); ax.set_ylim(-17, 17); ax.set_aspect("equal"); ax.set_xlabel("mm")
ax.legend(fontsize=7.5, loc="lower left", framealpha=.9)
ax.text(0.98, 0.98, f"{len(xy)} points, closed ✓\naxis inside ✓\nØ {S[2.0]['d_area']:.1f} < {1.6 * ao['d_area']:.1f} mm ✓",
        transform=ax.transAxes, ha="right", va="top", fontsize=9)

# ---------- 7 open chain at 4 mm ----------
ax = panel("At 4 mm the chain does not close: two loose ends\nwhere the plane runs into the mitral opening")
o4 = org + 4.0 * n
sl4 = endo.slice(normal=n, origin=o4)
loops = chain_loops(np.asarray(sl4.points), sl4.lines, o4, e1, e2)
for closed, q in loops:
    ax.plot(*q.T, "-", color="#009E73", lw=1.8)
    ax.scatter(*q.T, s=6, color="#009E73")
    if not closed:
        ax.scatter(*q[[0, -1]].T, s=110, marker="X", color="k", zorder=4)
ax.scatter([], [], s=80, marker="X", color="k", label="loose ends of the chain")
ax.plot([], [], color="#009E73", lw=1.8, label="cut of the LV endocardium at 4 mm")
ax.annotate("", xy=u * 21, xytext=u * 14, arrowprops=dict(arrowstyle="->", color="tab:red", lw=1.5))
ax.text(*(u * 23), "mitral\nopening", color="tab:red", ha="center", va="center", fontsize=8)
ax.plot(0, 0, "k+", ms=12, mew=2)
ax.set_xlim(-26, 26); ax.set_ylim(-26, 26); ax.set_aspect("equal"); ax.set_xlabel("mm")
ax.legend(fontsize=8, loc="lower right", framealpha=.9)

# ---------- 8 chord closure ----------
ax = panel("Keep wall points near the axis, find the largest empty\nangle (gap), bridge it with a straight chord")
s4 = S[4.0]
rmax = cfg["measure"]["cc_radius_factor"] * ao["d_area"]
axis_pt = np.zeros(2)
kept2 = [s for s in G["sections"] if s.get("kept")]
from src.lvot import resample_loop  # noqa: E402
axis_pt = resample_loop(kept2[-1]["xy"]).mean(0)
X = np.vstack([q for _, q in loops])
far = np.linalg.norm(X - axis_pt, axis=1) >= rmax
ax.scatter(*X[far].T, s=5, color="0.75", label="wall points too far from the axis (ignored)")
arc = s4["cc_xy"]
ax.add_patch(Polygon(arc, closed=True, fc="#009E73", alpha=.18, ec="none"))
ax.scatter(*arc.T, s=6, color="#009E73", label="wall points used")
ang = np.degrees(np.arctan2(*(arc - axis_pt).T[::-1]))
occ = np.zeros(72, bool); occ[np.floor((ang + 180) / 5).astype(int).clip(0, 71)] = True
e = np.concatenate([occ, occ]); best = cur = end = 0
for i, v_ in enumerate(e):
    cur = 0 if v_ else cur + 1
    if cur > best:
        best, end = cur, i
g0 = ((end - best + 1) % 72) * 5 - 180
ax.add_patch(Wedge(axis_pt, rmax, g0, g0 + best * 5, color="tab:orange", alpha=.18, lw=0, label=f"gap {s4['gap_deg']:.0f}° (≤ 120° allowed)"))
ax.add_patch(Circle(axis_pt, rmax, fill=False, ls=":", color="0.4"))
ax.plot(*np.vstack([arc[-1], arc[0]]).T, "-", color="#E69F00", lw=3, label="chord")
ax.plot(*axis_pt, "k+", ms=12, mew=2)
ax.set_xlim(-26, 26); ax.set_ylim(-26, 26); ax.set_aspect("equal"); ax.set_xlabel("mm")
ax.legend(fontsize=8, loc="lower right", framealpha=.9)

# ---------- 9 polygon -> area ----------
ax = panel("Ordered points = polygon → area by the shoelace formula\n(sum of triangles from the centre), then diameters")
xy = loop2
A, Pm = polygon_area_perim(xy)
el = moment_ellipse(xy)
c = np.array([el["cx"], el["cy"]])
ax.add_patch(Polygon(xy, closed=True, fc="#0072B2", alpha=.15, ec="#0072B2", lw=2))
step = max(1, len(xy) // 24)
for i in range(0, len(xy), step):
    j = (i + step) % len(xy)
    ax.add_patch(Polygon([c, xy[i], xy[j]], closed=True, fc="#0072B2" if (i // step) % 2 else "#56B4E9",
                         alpha=.25, ec="0.4", lw=.4))
ax.plot(*c, "k+", ms=10)
ax.set_xlim(-14, 14); ax.set_ylim(-14, 14); ax.set_aspect("equal"); ax.set_xlabel("mm")
ax.text(0.02, 0.02,
        f"A = ½ |Σ (xᵢ·yᵢ₊₁ − xᵢ₊₁·yᵢ)| = {A:.0f} mm²\n"
        f"area-derived Ø = 2√(A/π) = {2 * np.sqrt(A / np.pi):.1f} mm\n"
        f"moment ellipse {el['dmax']:.1f} × {el['dmin']:.1f} mm → ratio {el['dmax'] / el['dmin']:.2f}\n"
        f"echo PLAX Ø {S[2.0]['d_plax']:.1f} mm",
        transform=ax.transAxes, fontsize=9.5, va="bottom", bbox=dict(fc="white", ec="0.7"))

fig.suptitle(f"From mesh to LVOT outline, step by step — {mid} (end-diastole)", fontsize=14, y=0.995)
fig.tight_layout(rect=(0, 0, 1, 0.985))
out = ROOT / "outputs/figures/qc" / f"outline_steps_{mid}.png"
fig.savefig(out, dpi=110)
print("wrote", out)
log_run(__file__, mesh=mid, output=str(out.relative_to(ROOT)))
