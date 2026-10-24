"""STEP 5.3 — Illustrate how open LVOT sections are chord-closed (one mesh, several depths).

Usage: python scripts/04_gap_demo.py <mesh_id>  -> outputs/figures/qc/gap_demo_<id>.png
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Circle, Polygon, Wedge

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src.config import load_config  # noqa: E402
from src.lvot import chain_loops, measure_mesh, resample_loop  # noqa: E402
from src.provenance import log_run  # noqa: E402

mid = sys.argv[1]
labels = pd.read_csv(ROOT / load_config()["labels"]).set_index("mesh_id")
r = measure_mesh(labels.loc[mid, "folder"], keep_geometry=True)
G = r["_geom"]
ao, mv = G["rings"]["ao"], G["rings"]["mv"]
n, e1, e2, org = ao["normal"], ao["e1"], ao["e2"], G["origin"]
rmax = load_config()["measure"]["cc_radius_factor"] * ao["d_area"]
kept = [s for s in G["sections"] if s.get("kept")]
axis_pt = resample_loop(kept[-1]["xy"]).mean(0)
md = mv["centre"] - org
mang = np.arctan2(md @ e2, md @ e1)


def ell_xy(el, m=200):
    t = np.linspace(0, 2 * np.pi, m)
    a, b, th = el["dmax"] / 2, el["dmin"] / 2, el["theta"]
    x, y = a * np.cos(t), b * np.sin(t)
    return np.column_stack([el["cx"] + x * np.cos(th) - y * np.sin(th), el["cy"] + x * np.sin(th) + y * np.cos(th)])


depths = [2, 4, 6, 8]
fig, axs = plt.subplots(1, 4, figsize=(18, 5.2))
for ax, dep in zip(axs, depths):
    s = next(q for q in G["sections"] if q["depth"] == dep)
    o = org + dep * n
    sl = G["endo"].slice(normal=n, origin=o)
    loops = chain_loops(np.asarray(sl.points), sl.lines, o, e1, e2)
    for _, xy in loops:  # full intersection of the plane with the LV endocardium
        ax.plot(*xy.T, "-", color="0.75", lw=1, zorder=1)
    ax.add_patch(Circle(axis_pt, rmax, fill=False, ls=":", color="0.5", lw=.8))
    ax.plot(*axis_pt, "k+", ms=8)
    ax.annotate("", xy=axis_pt + 1.15 * rmax * np.array([np.cos(mang), np.sin(mang)]),
                xytext=axis_pt + 0.75 * rmax * np.array([np.cos(mang), np.sin(mang)]),
                arrowprops=dict(arrowstyle="->", color="tab:red"))
    ax.text(*(axis_pt + 1.2 * rmax * np.array([np.cos(mang), np.sin(mang)])), "mitral", color="tab:red",
            ha="center", va="center", fontsize=8)
    if s.get("kept"):
        xy = s["xy"]
        ax.add_patch(Polygon(xy, closed=True, fc="#0072B2", alpha=.2, ec="#0072B2", lw=2))
        title = f"{dep} mm: closed loop (strict)\nd_area {s['d_area']:.1f} mm, {s['dmax']:.1f} × {s['dmin']:.1f}"
        el = s["ell"]
    else:
        X = np.vstack([xy for _, xy in loops])
        X = X[np.linalg.norm(X - axis_pt, axis=1) < rmax]
        ax.scatter(*X.T, s=5, color="k", zorder=3, label="wall points used (inside dotted circle)")
        # gap sector: complement of the occupied 5-degree bins
        ang = np.degrees(np.arctan2(*(X - axis_pt).T[::-1]))
        occ = np.zeros(72, bool)
        occ[np.floor((ang + 180) / 5).astype(int).clip(0, 71)] = True
        e = np.concatenate([occ, occ]); best = cur = end = 0
        for i, v in enumerate(e):
            cur = 0 if v else cur + 1
            if cur > best:
                best, end = cur, i
        g0 = ((end - best + 1) % 72) * 5 - 180
        ax.add_patch(Wedge(axis_pt, rmax, g0, g0 + best * 5, color="tab:orange", alpha=.15, lw=0,
                           label=f"gap = largest empty sector ({best*5}°)"))
        if s.get("cc"):
            arc = s["cc_xy"]
            ax.add_patch(Polygon(arc, closed=True, fc="#009E73", alpha=.2, ec="none"))
            ax.plot(*np.vstack([arc[-1], arc[0]]).T, "-", color="tab:orange", lw=2.5, label="chord (bridges the gap)")
            title = (f"{dep} mm: open, gap {s['gap_deg']:.0f}° → chord-closed\n"
                     f"d_area {s['cc_d_area']:.1f} mm, {s['cc_dmax']:.1f} × {s['cc_dmin']:.1f}")
            el = s["cc_ell"]
        else:
            title = f"{dep} mm: open, gap {s['gap_deg']:.0f}° > 120° → rejected"
            el = None
        ax.legend(fontsize=7, loc="lower left")
    if el is not None:
        ax.plot(*ell_xy(el).T, "--", color="k", lw=.8)
    ax.set_title(title, fontsize=9)
    ax.set_aspect("equal")
    ax.set_xlim(axis_pt[0] - 26, axis_pt[0] + 26); ax.set_ylim(axis_pt[1] - 26, axis_pt[1] + 26)
    ax.set_xticks([]); ax.set_yticks([])
fig.suptitle(f"{mid}: LVOT cross-sections viewed down the aortic axis. Grey = full plane ∩ LV endocardium; "
             "dashed = moment-equivalent ellipse", fontsize=10)
fig.tight_layout()
out = ROOT / "outputs/figures/qc" / f"gap_demo_{mid}.png"
fig.savefig(out, dpi=120)
print("wrote", out)
log_run(__file__, mesh=mid, output=str(out.relative_to(ROOT)))
