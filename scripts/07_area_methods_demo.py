"""STEP 6.7 — Illustrate the five area estimates on one LVOT outline, plus population averages.

Top row (one mesh, strict LVOT section at 2 mm, looking down the aortic axis):
  1 moment ellipse        same area spread as the outline (best-possible ellipse)
  2 Feret ellipse         widest x narrowest calliper width (CT / 3D max-min approach)
  3 CT ellipse            longest diameter through the centre x the one at 90 deg to it (CT reading)
  4 two-view ellipse      PLAX diameter x perpendicular diameter (echo biplane approach)
  5 echo circle           PLAX diameter only, circular LVOT assumed (standard echo)
Bottom: mean % difference from the outline area by method, level and sex (all meshes).
Usage: python scripts/07_area_methods_demo.py [mesh_id] -> outputs/figures/analysis/fig_area_methods_demo.png
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Ellipse, Polygon
from matplotlib.path import Path as MPath

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src.config import load_config  # noqa: E402
from src.lvot import ct_diams, measure_mesh, moment_ellipse, polygon_area_perim  # noqa: E402
from src.provenance import log_run  # noqa: E402

mid = sys.argv[1] if len(sys.argv) > 1 else "female_age_65_bmi_26"
cfg = load_config()
COL = cfg["analysis"]["colours"]
labels = pd.read_csv(ROOT / cfg["labels"]).set_index("mesh_id")
r = measure_mesh(labels.loc[mid, "folder"], keep_geometry=True)
G = r["_geom"]
ao, mv = G["rings"]["ao"], G["rings"]["mv"]
n, e1, e2, org = ao["normal"], ao["e1"], ao["e2"], G["origin"]
sec = next(s for s in G["sections"] if s["depth"] == 2.0)
xy = sec["xy"]
A = polygon_area_perim(xy)[0]
el = moment_ellipse(xy)
c = np.array([el["cx"], el["cy"]])
md = mv["centre"] - org
u = np.array([md @ e1, md @ e2]); u /= np.linalg.norm(u)  # PLAX direction (towards mitral)
v = np.array([-u[1], u[0]])


def line_through(p, d):
    """End points of the line through p along d, clipped to the outline."""
    t = np.linspace(-30, 30, 12001)
    ins = MPath(xy).contains_points(p + t[:, None] * d)
    return p + t[ins].min() * d, p + t[ins].max() * d


def feret(theta):
    d = np.array([np.cos(theta), np.sin(theta)])
    proj = xy @ d
    return proj.max() - proj.min(), proj.min(), proj.max(), d


angs = np.linspace(0, np.pi, 720, endpoint=False)
widths = np.array([feret(t)[0] for t in angs])
tmax, tmin = angs[widths.argmax()], angs[widths.argmin()]

GREY, EST, CAL = "#0072B2", "#D55E00", "0.15"
fig = plt.figure(figsize=(22, 10))
gs = fig.add_gridspec(2, 5, height_ratios=[1.25, 1])


def base(ax, title):
    ax.add_patch(Polygon(xy, closed=True, fc=GREY, alpha=.18, ec=GREY, lw=2, label=f"true outline area {A:.0f} mm²"))
    ax.annotate("", xy=c + u * 15, xytext=c + u * 11.5, arrowprops=dict(arrowstyle="->", color="tab:red"))
    ax.text(*(c + u * 16.8), "mitral", color="tab:red", ha="center", va="center", fontsize=8)
    ax.set_aspect("equal"); ax.set_xlim(c[0] - 17, c[0] + 17); ax.set_ylim(c[1] - 17, c[1] + 17)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(title, fontsize=10.5, loc="left")


def finish(ax, est, formula):
    pct = 100 * (est - A) / A
    ax.text(0.02, 0.02, f"{formula}\n= {est:.0f} mm²  ({pct:+.0f}% vs true)", transform=ax.transAxes,
            fontsize=9.5, va="bottom", bbox=dict(fc="white", ec="0.7", alpha=.9))
    ax.legend(fontsize=7.5, loc="upper right", framealpha=.9)


def cal(ax, p, q, lab, ls="-"):
    ax.plot(*np.vstack([p, q]).T, ls, color=CAL, lw=1.8, marker="|", ms=10, mew=2, label=lab)


# 1 moment ellipse
ax = fig.add_subplot(gs[0, 0])
base(ax, "1  Moment ellipse — best-possible ellipse\n(uses the whole outline; no scanner reports it)\naxes at 90°, fitted to the whole shape")
th = el["theta"]; dM = np.array([np.cos(th), np.sin(th)]); dm = np.array([-dM[1], dM[0]])
cal(ax, c - dM * el["dmax"] / 2, c + dM * el["dmax"] / 2, f"major axis {el['dmax']:.1f} mm")
cal(ax, c - dm * el["dmin"] / 2, c + dm * el["dmin"] / 2, f"minor axis {el['dmin']:.1f} mm", "--")
ax.add_patch(Ellipse(c, el["dmax"], el["dmin"], angle=np.degrees(th), fill=False, ec=EST, lw=2.2, ls="--",
                     label="assumed ellipse"))
finish(ax, np.pi / 4 * el["dmax"] * el["dmin"], "π/4 × major × minor")

# 2 Feret ellipse
ax = fig.add_subplot(gs[0, 1])
base(ax, "2  Feret ellipse\n(widest and narrowest calliper widths)")
for t_, lab, ls in ((tmax, "widest width", "-"), (tmin, "narrowest width", "--")):
    w, lo, hi, d = feret(t_)
    perp = np.array([-d[1], d[0]])
    off = (c @ perp)
    p, q = lo * d + off * perp, hi * d + off * perp
    cal(ax, p, q, f"{lab} {w:.1f} mm", ls)
    for e_ in (lo, hi):  # calliper jaws
        j = e_ * d + off * perp
        ax.plot(*np.vstack([j - perp * 4, j + perp * 4]).T, color=CAL, lw=1)
Fmax, Fmin = widths.max(), widths.min()
ax.add_patch(Ellipse(c, Fmax, Fmin, angle=np.degrees(tmax), fill=False, ec=EST, lw=2.2, ls="--", label="assumed ellipse"))
ang_between = abs(np.degrees(np.angle(np.exp(2j * (tmax - tmin))))) / 2
ax.set_title(ax.get_title(loc="left") + f"\nthe two widths are {ang_between:.0f}° apart, not 90°", fontsize=10.5, loc="left")
finish(ax, np.pi / 4 * Fmax * Fmin, "π/4 × widest × narrowest")

# 3 CT ellipse (longest through the centre x perpendicular)
ax = fig.add_subplot(gs[0, 2])
base(ax, "3  CT ellipse — CT reading\n(longest diameter through the centre +\nthe diameter at 90° to it)")
cd = ct_diams(xy)
dl = np.array([np.cos(cd["ct_theta"]), np.sin(cd["ct_theta"])]); ds = np.array([-dl[1], dl[0]])
pl, ql = line_through(c, dl); ps, qs = line_through(c, ds)
cal(ax, pl, ql, f"longest {cd['ct_long']:.1f} mm")
cal(ax, ps, qs, f"perpendicular {cd['ct_short']:.1f} mm", "--")
ax.add_patch(Ellipse(c, cd["ct_long"], cd["ct_short"], angle=np.degrees(cd["ct_theta"]), fill=False, ec=EST, lw=2.2,
                     ls="--", label="assumed ellipse"))
finish(ax, np.pi / 4 * cd["ct_long"] * cd["ct_short"], "π/4 × longest × perpendicular")

# 3 two-view (PLAX x perpendicular)
ax = fig.add_subplot(gs[0, 3])
base(ax, "4  Two-view ellipse — echo biplane style\n(PLAX diameter + the one at 90° to it)")
p1, q1 = line_through(c, u); p2, q2 = line_through(c, v)
Dp, Dq = np.linalg.norm(q1 - p1), np.linalg.norm(q2 - p2)
cal(ax, p1, q1, f"PLAX diameter {Dp:.1f} mm")
cal(ax, p2, q2, f"perpendicular {Dq:.1f} mm", "--")
ax.add_patch(Ellipse(c, Dp, Dq, angle=np.degrees(np.arctan2(u[1], u[0])), fill=False, ec=EST, lw=2.2, ls="--",
                     label="assumed ellipse"))
off_deg = abs(np.degrees(np.angle(np.exp(2j * (np.arctan2(u[1], u[0]) - (th + np.pi / 2)))))) / 2
ax.set_title(ax.get_title(loc="left") + f"\nthe PLAX line is {off_deg:.0f}° off the true minor axis", fontsize=10.5, loc="left")
finish(ax, np.pi / 4 * Dp * Dq, "π/4 × D_PLAX × D_perp")

# 4 echo circle
ax = fig.add_subplot(gs[0, 4])
base(ax, "5  Echo circle — standard echo\n(PLAX diameter only, LVOT assumed round)\nthe long side-to-side axis is ignored")
cal(ax, p1, q1, f"PLAX diameter {Dp:.1f} mm")
ax.add_patch(Ellipse((p1 + q1) / 2, Dp, Dp, fill=True, fc=EST, alpha=.12, ec=EST, lw=2.2, ls="--",
                     label="assumed circle"))
finish(ax, np.pi / 4 * Dp ** 2, "π/4 × D_PLAX²")

# bottom: population means
S = pd.read_csv(ROOT / "outputs/tables/area_methods_summary.csv")
M = [("ellipse_moment", "1 Moment\nellipse"), ("ellipse_feret", "2 Feret\nellipse"), ("ellipse_ct", "3 CT\nellipse"),
     ("ellipse_2view", "4 Two-view\nellipse"), ("echo_circle", "5 Echo\ncircle")]
sub = gs[1, :].subgridspec(1, 3)
LEV = ["Aortic annulus", "LVOT 2 mm (closed)", "LVOT 4 mm (chord-closed)"]
for j, lev in enumerate(LEV):
    ax = fig.add_subplot(sub[0, j])
    x = np.arange(len(M))
    for k, sx in enumerate("FM"):
        row = S[(S.level == lev) & (S.sex == sx)].iloc[0]
        vals = [row[f"{m}_pct_mean"] for m, _ in M]
        sds = [row[f"{m}_pct_sd"] for m, _ in M]
        b = ax.bar(x + (k - .5) * .36, vals, width=.34, color=COL[sx], yerr=sds, capsize=2,
                   error_kw=dict(lw=.8, ecolor="0.4"), zorder=3, label={"F": "Women", "M": "Men"}[sx])
        for xi, val in zip(x + (k - .5) * .36, vals):
            ax.text(xi, val + (1.2 if val >= 0 else -1.2), f"{val:+.0f}", ha="center",
                    va="bottom" if val >= 0 else "top", fontsize=8, color="0.2")
    ax.axhline(0, color="0.3", lw=1)
    ax.axhspan(-7, 7, color="0.9", lw=0, zorder=0)
    ax.set_xticks(x); ax.set_xticklabels([l for _, l in M], fontsize=8.5)
    ax.set_ylim(-45, 15); ax.set_ylabel("% vs true area (mean ± SD)")
    ax.set_title(f"All meshes — {lev}", fontsize=10, loc="left")
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)
    if j == 0:
        ax.legend(fontsize=8, frameon=False, loc="lower left")
    nct = S[(S.level == lev)].n_ct.sum(); nall = S[(S.level == lev)].n.sum()
    if nct < nall:
        ax.text(0.99, 0.03, f"CT ellipse: {nct} of {nall} meshes so far", transform=ax.transAxes, ha="right", fontsize=8, color="0.3")
    if j == 2:
        ax.text(0.01, 0.97, "grey band = ±7%: all ellipse methods fall inside it; the circle does not",
                transform=ax.transAxes, fontsize=8.5, va="top")
fig.suptitle(f"Five ways to estimate LVOT area from diameters — top: {mid}, LVOT 2 mm below the annulus "
             "(end-diastole), looking down the aortic axis", fontsize=11.5, y=0.995)
fig.tight_layout(rect=(0, 0, 1, 0.97))
out = ROOT / "outputs/figures/analysis/fig_area_methods_demo.png"
fig.savefig(out, dpi=125)
print("wrote", out)
log_run(__file__, mesh=mid, output=str(out.relative_to(ROOT)))
