"""STEP 5.4 — Illustrate how an LVOT outline is obtained (one mesh).

(a) 3D: LV endocardium around the LVOT (end-diastole), aortic ring, and the cutting planes at
    2 mm and 4 mm below the ring band, with the curves where each plane cuts the wall.
(b) 2 mm, en face: the cut is a closed chain of points -> polygon -> area.
(c) 4 mm, en face: the cut is open on the mitral side -> straight chord closes it -> polygon -> area.
Usage: python scripts/05_outline_demo.py <mesh_id>  -> outputs/figures/qc/outline_demo_<id>.png
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Polygon
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src.config import load_config  # noqa: E402
from src.lvot import measure_mesh, polygon_area_perim  # noqa: E402
from src.provenance import log_run  # noqa: E402

mid = sys.argv[1] if len(sys.argv) > 1 else "female_age_65_bmi_26"
labels = pd.read_csv(ROOT / load_config()["labels"]).set_index("mesh_id")
r = measure_mesh(labels.loc[mid, "folder"], keep_geometry=True)
G = r["_geom"]
ao, mv = G["rings"]["ao"], G["rings"]["mv"]
n, e1, e2, org = ao["normal"], ao["e1"], ao["e2"], G["origin"]
S = {s["depth"]: s for s in G["sections"]}
s2, s4 = S[2.0], S[4.0]
C2, C4 = "#0072B2", "#009E73"  # 2 mm blue, 4 mm green
CHORD = "#E69F00"


def frame(X):
    """mm coordinates in the ring frame: x, y in the ring plane, z up towards the aorta."""
    D = X - org
    return np.column_stack([D @ e1, D @ e2, -(D @ n)])


md = mv["centre"] - org
mdir = np.array([md @ e1, md @ e2]); mdir /= np.linalg.norm(mdir)

fig = plt.figure(figsize=(17, 7))
# ---------- (a) 3D ----------
ax = fig.add_subplot(1, 3, 1, projection="3d")
endo = G["endo"]
tri = endo.faces.reshape(-1, 4)[:, 1:]
V = frame(np.asarray(endo.points))
cen = V[tri].mean(1)
keep = (np.hypot(cen[:, 0], cen[:, 1]) < 17) & (cen[:, 2] > -9) & (cen[:, 2] < 3)
# cut away the half facing the viewer so the inside of the tube is visible
perp = np.array([-mdir[1], mdir[0]])
keep &= (cen[:, :2] @ perp) > -2
ax.add_collection3d(Poly3DCollection(V[tri[keep]], facecolor="0.75", edgecolor="0.55", linewidths=0.1, alpha=0.5))
cxy = np.array([ao["ell"]["cx"], ao["ell"]["cy"]])
rp = ao["poly"] - cxy  # ring outline relative to the ring centre (= section origin in-plane)
rp = np.vstack([rp, rp[:1]])
ax.plot(rp[:, 0], rp[:, 1], np.full(len(rp), ao["band_max_depth"]), "k-", lw=2.5, label="aortic ring (plug–LV contact)")
for dep, col, sec in ((2, C2, s2), (4, C4, s4)):
    L = 15
    sq = np.array([[-L, -L], [L, -L], [L, L], [-L, L]])
    ax.add_collection3d(Poly3DCollection([np.column_stack([sq, np.full(4, -dep)])], facecolor=col, alpha=0.08, edgecolor=col))
    xy = sec["xy"] if sec.get("kept") else sec["cc_xy"]
    if sec.get("kept"):
        xy = np.vstack([xy, xy[:1]])
    ax.plot(xy[:, 0], xy[:, 1], np.full(len(xy), -dep), "-", color=col, lw=2.5,
            label=f"plane {dep} mm below ring ∩ LV wall")
    if not sec.get("kept"):
        ch = np.vstack([xy[-1], xy[0]])
        ax.plot(ch[:, 0], ch[:, 1], [-dep, -dep], "-", color=CHORD, lw=3, label="chord across the mitral-side gap")
ax.quiver(0, 0, -6, *(mdir * 14), 0, color="tab:red", arrow_length_ratio=0.15)
ax.text(*(mdir * 16), -6, "mitral", color="tab:red", fontsize=9)
ax.quiver(0, 0, 3, 0, 0, 6, color="k", arrow_length_ratio=0.2)
ax.text(0, 0, 9.5, "aorta", fontsize=9, ha="center")
ax.set_xlim(-16, 16); ax.set_ylim(-16, 16); ax.set_zlim(-9, 9)
ax.set_box_aspect((1, 1, 0.6))
az = np.degrees(np.arctan2(perp[1], perp[0]))
ax.view_init(elev=22, azim=az - 200)
ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])
ax.set_title("(a) LV endocardium around the LVOT (cut open)\nflat planes perpendicular to the aortic axis", fontsize=10)
ax.legend(fontsize=7.5, loc="lower center", bbox_to_anchor=(0.5, -0.02), ncol=2, frameon=False)


def enface(ax, sec, dep, col, title):
    ax.plot(*(rp.T), "-", color="0.6", lw=1, label="aortic ring (projected)")
    if sec.get("kept"):
        xy = sec["xy"]
        area = polygon_area_perim(xy)[0]
        ax.add_patch(Polygon(xy, closed=True, fc=col, alpha=0.18, ec="none"))
        ax.plot(*np.vstack([xy, xy[:1]]).T, "-", color=col, lw=1.2)
        ax.plot(*xy.T, "o", color=col, ms=2.5, label=f"cut points ({len(xy)}): closed chain")
        sub = f"closed outline → polygon\narea {area:.0f} mm² → area-derived diameter {sec['d_area']:.1f} mm"
    else:
        xy = sec["cc_xy"]
        area = polygon_area_perim(xy)[0]
        ax.add_patch(Polygon(xy, closed=True, fc=col, alpha=0.18, ec="none"))
        ax.plot(*xy.T, "o", color=col, ms=2.5, label=f"cut points ({len(xy)}): open chain")
        ch = np.vstack([xy[-1], xy[0]])
        ax.plot(*ch.T, "-", color=CHORD, lw=3, label=f"chord bridging the {sec['gap_deg']:.0f}° gap")
        ax.plot(*ch.T, "o", color=CHORD, ms=7, mfc="white", mew=2)
        sub = (f"open on the mitral side → closed by a chord\n"
               f"area {area:.0f} mm² → area-derived diameter {sec['cc_d_area']:.1f} mm")
    ax.annotate("", xy=mdir * 19, xytext=mdir * 12, arrowprops=dict(arrowstyle="->", color="tab:red", lw=1.5))
    ax.text(*(mdir * 21), "mitral", color="tab:red", ha="center", va="center", fontsize=9)
    ax.set_title(f"{title}\n{sub}", fontsize=10)
    ax.set_aspect("equal"); ax.set_xlim(-22, 22); ax.set_ylim(-22, 22)
    ax.set_xlabel("mm"); ax.grid(alpha=.2)
    ax.legend(fontsize=7.5, loc="lower right")


enface(fig.add_subplot(1, 3, 2), s2, 2, C2, "(b) 2 mm below the ring, looking down the aortic axis")
enface(fig.add_subplot(1, 3, 3), s4, 4, C4, "(c) 4 mm below the ring, looking down the aortic axis")
fig.suptitle(f"How an LVOT outline is obtained — {mid} (end-diastole). Shaded = region whose area is measured.",
             fontsize=11, y=0.99)
fig.tight_layout(rect=(0, 0, 1, 0.95))
out = ROOT / "outputs/figures/qc" / f"outline_demo_{mid}.png"
fig.savefig(out, dpi=130)
print("wrote", out)
log_run(__file__, mesh=mid, output=str(out.relative_to(ROOT)))
