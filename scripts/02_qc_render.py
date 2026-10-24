"""STEP 5.2 — QC render: aortic ring, mitral direction, strict LVOT sections and chord-closed sections, per mesh.

Left: en-face view along the aortic-ring normal (looking from the aorta into the LV).
Right: side view (plane containing the ring normal and the mitral direction) of LV endo nodes near the LVOT.
Usage: python scripts/02_qc_render.py <mesh_id> [...]  (or --subset)  -> outputs/figures/qc/<id>.png
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src.config import load_config  # noqa: E402
from src.lvot import measure_mesh  # noqa: E402
from src.provenance import log_run  # noqa: E402

cfg = load_config()
labels = pd.read_csv(ROOT / cfg["labels"]).set_index("mesh_id")
ids = cfg["subset"] if sys.argv[1] == "--subset" else sys.argv[1:]
odir = ROOT / "outputs/figures/qc"
odir.mkdir(parents=True, exist_ok=True)


def ell_xy(el, n=200):
    t = np.linspace(0, 2 * np.pi, n)
    a, b, th = el["dmax"] / 2, el["dmin"] / 2, el["theta"]
    x, y = a * np.cos(t), b * np.sin(t)
    return np.column_stack([el["cx"] + x * np.cos(th) - y * np.sin(th), el["cy"] + x * np.sin(th) + y * np.cos(th)])


for mid in ids:
    r = measure_mesh(labels.loc[mid, "folder"], keep_geometry=True)
    G = r["_geom"]
    ao, mv = G["rings"]["ao"], G["rings"]["mv"]
    n, e1, e2, org = ao["normal"], ao["e1"], ao["e2"], G["origin"]
    fig, axs = plt.subplots(1, 2, figsize=(12, 5.6))
    ax = axs[0]
    # ring (projected relative to section origin, which is on the ring axis)
    rc = ao["centre"] - (org - ao["band_max_depth"] * n)  # zero
    ax.scatter(ao["xy"][:, 0], ao["xy"][:, 1], s=4, c="k", alpha=.4, label=f"Ao ring nodes (d_area {ao['d_area']:.1f})")
    ax.plot(*ell_xy(ao["ell"]).T, "k-", lw=1.5)
    cmap = plt.get_cmap("viridis")
    for s in G["sections"]:
        col = cmap(s["depth"] / 10)
        if s["depth"] > 10:
            continue
        if s.get("kept"):
            xy = s["xy"]
            if s["depth"] % 2 == 0:
                ax.plot(*ell_xy(s["ell"]).T, "-.", color=col, lw=.6)
            ax.plot(*np.vstack([xy, xy[:1]]).T, "-", color=col, lw=1.2,
                    label=f"{s['depth']:g} mm closed: {s['d_area']:.1f}, ratio {s['dmax']/s['dmin']:.2f}" if s["depth"] % 2 == 0 else None)
        elif s.get("cc"):
            xy = s["cc_xy"]
            ax.plot(*xy.T, ".", color=col, ms=2)
            ax.plot(*np.vstack([xy[-1], xy[0]]).T, ":", color=col, lw=1.5)
            if s["depth"] % 2 == 0:
                ax.plot(*ell_xy(s["cc_ell"]).T, "--", color=col, lw=.8,
                        label=f"{s['depth']:g} mm chord-closed: {s['cc_d_area']:.1f}, gap {s['gap_deg']:.0f}°")
    md = mv["centre"] - org
    mxy = np.array([md @ e1, md @ e2])
    mxy = mxy / np.linalg.norm(mxy) * 16
    ax.annotate("to mitral", xy=mxy * 1.25, xytext=mxy * 0.75, arrowprops=dict(arrowstyle="->", color="tab:red"),
                color="tab:red", ha="center")
    ax.set_aspect("equal"); ax.set_xlim(-22, 22); ax.set_ylim(-22, 22)
    ax.set_title("En face along aortic-ring normal (mm)")
    ax.legend(fontsize=6.5, loc="lower left")
    # side view
    ax = axs[1]
    u = mxy / np.linalg.norm(mxy); u3 = u[0] * e1 + u[1] * e2
    E = np.asarray(G["endo"].points) - org
    dep, side, perp = E @ n, E @ u3, E @ np.cross(n, u3)
    sel = (np.abs(perp) < 1.5) & (dep > -8) & (dep < 30) & (np.abs(side) < 40)
    ax.scatter(side[sel], -dep[sel], s=1, c="tab:blue", alpha=.5, label="LV endo (|off-plane| < 1.5 mm)")
    for name, R, c in (("Ao ring", ao, "k"), ("MV ring", mv, "tab:red")):
        Xr = R["xy"][:, 0][:, None] * R["e1"] + R["xy"][:, 1][:, None] * R["e2"] + (R["centre"] - R["xy"].mean(0)[0] * R["e1"] - R["xy"].mean(0)[1] * R["e2"]) - org
        ax.scatter(Xr @ u3, -(Xr @ n), s=2, c=c, alpha=.4, label=name)
    for d in (0, 2, 4, 6, 8, 10):
        ax.axhline(-d, color=cmap(d / 10), lw=.6, ls="--")
    ap = G["apex"] - org
    ax.annotate("to apex", xy=(ap @ u3 / np.linalg.norm(ap) * 25, -(ap @ n) / np.linalg.norm(ap) * 25), xytext=(0, 0),
                arrowprops=dict(arrowstyle="->"))
    ax.set_aspect("equal"); ax.set_title("Side view: depth below ring (down) vs towards mitral (right)")
    ax.legend(fontsize=7, loc="lower right")
    fig.suptitle(f"{mid}  Ao {r['ao_d_area']:.1f} mm (ratio {r['ao_ellip_ratio']:.2f})  MV {r['mv_d_area']:.1f} mm  "
                 f"Ao–MV {r['ao_mv_angle']:.0f}°  Ao axis–LV long axis {r['ao_axis_vs_lv_long_angle']:.0f}°  "
                 f"closed tube {r['lvot_tube_len']} mm  flags: {r['qc_flags'] or '-'}", fontsize=9)
    fig.tight_layout()
    fig.savefig(odir / f"{mid}.png", dpi=110)
    plt.close(fig)
    print("wrote", odir / f"{mid}.png")
log_run(__file__, meshes=list(ids), output_dir=str(odir.relative_to(ROOT)))
