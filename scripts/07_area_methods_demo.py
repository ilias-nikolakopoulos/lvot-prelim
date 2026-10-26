"""STEP 6.7 — Illustrate the five area estimates on one LVOT outline, plus population averages.

Top row (one mesh, strict LVOT section at 2 mm, looking down the aortic axis):
  1 moment ellipse        same area spread as the outline (best-possible ellipse)
  2 Feret ellipse         widest x narrowest calliper width (CT / 3D max-min approach)
  3 minor-major axis ellipse  longest diameter through the centre x the one at 90 deg to it (as in a CT reading)
  4 two-view ellipse      PLAX diameter x perpendicular diameter (echo biplane approach)
  5 echo circle           PLAX diameter only, circular LVOT assumed (standard echo)
Bottom: mean % difference from the outline area by method, level and sex (matched same-age, same-BMI pairs).
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
from matplotlib import font_manager  # noqa: E402

# Roboto (Apache 2.0, googlefonts/roboto v2.138) shipped in assets/fonts so the figure renders the same anywhere
for f_ in sorted((ROOT / "assets/fonts").glob("Roboto-*.ttf")):
    font_manager.fontManager.addfont(str(f_))
plt.rcParams["font.family"] = "Roboto"

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

GREY, EST, CAL = "#0072B2", "#B94F00", "0.15"
NAME_COL = "#102D40"  # dark blue-grey for method names
SHOW_TWO_VIEW = False  # PLAX x perpendicular ellipse: kept in the code and tables, hidden from this figure
LEVEL_COL = {"Aortic annulus": "#E6C229", "LVOT 2 mm (closed)": "#009E73", "LVOT 4 mm (chord-closed)": "#7B4FA6"}
LEVEL_TITLE = {"Aortic annulus": "Aortic annulus", "LVOT 2 mm (closed)": "2 mm from annulus",
               "LVOT 4 mm (chord-closed)": "4 mm from annulus"}
ncol = 5 if SHOW_TWO_VIEW else 4
fig = plt.figure(figsize=(4.4 * ncol, 10))
gs = fig.add_gridspec(2, ncol, height_ratios=[1.25, 1])
panels = iter(range(ncol))


def base(ax, name, desc):
    """Bold dark-blue method name with a plain description underneath, above the panel."""
    ax.add_patch(Polygon(xy, closed=True, fc=GREY, alpha=.18, ec=GREY, lw=2, label=f"true outline area {A:.0f} mm²"))
    ax.annotate("", xy=c + u * 13.5, xytext=c + u * 10.5, arrowprops=dict(arrowstyle="->", color="tab:red"))
    ax.text(*(c + u * 13.5 + np.array([0, 1.4])), "mitral", color="tab:red", ha="center", va="bottom", fontsize=9)
    ax.set_aspect("equal"); ax.set_xlim(c[0] - 17, c[0] + 17); ax.set_ylim(c[1] - 17, c[1] + 17)
    ax.set_xticks([]); ax.set_yticks([])
    nl = 2  # fixed slot so all method names sit at the same height
    ax.annotate(desc, xy=(0, 1), xycoords="axes fraction", xytext=(0, 4), textcoords="offset points",
                fontsize=10, va="bottom", ha="left", color="0.15")
    ax.annotate(name, xy=(0, 1), xycoords="axes fraction", xytext=(0, 8 + 12.5 * nl), textcoords="offset points",
                fontsize=12.5, fontweight="bold", color=NAME_COL, va="bottom", ha="left")


def shape_label(est, kind="ellipse"):
    return f"assumed {kind}: {est:.0f} mm² ({100 * (est - A) / A:+.0f}%)"


def finish(ax):
    ax.legend(fontsize=7.5, loc="upper right", framealpha=.9)


def cal(ax, p, q, lab, ls="-"):
    ax.plot(*np.vstack([p, q]).T, ls, color=CAL, lw=1.8, marker="|", ms=10, mew=2, label=lab)


num = iter(range(1, ncol + 1))

# moment ellipse
ax = fig.add_subplot(gs[0, next(panels)])
base(ax, f"{next(num)}  Moment ellipse", "best-possible ellipse fitted to the whole outline")
th = el["theta"]; dM = np.array([np.cos(th), np.sin(th)]); dm = np.array([-dM[1], dM[0]])
cal(ax, c - dM * el["dmax"] / 2, c + dM * el["dmax"] / 2, f"major axis {el['dmax']:.1f} mm")
cal(ax, c - dm * el["dmin"] / 2, c + dm * el["dmin"] / 2, f"minor axis {el['dmin']:.1f} mm", "--")
ax.add_patch(Ellipse(c, el["dmax"], el["dmin"], angle=np.degrees(th), fill=False, ec=EST, lw=2.2, ls="--",
                     label=shape_label(np.pi / 4 * el["dmax"] * el["dmin"])))
finish(ax)

# Feret ellipse
ax = fig.add_subplot(gs[0, next(panels)])
ang_between = abs(np.degrees(np.angle(np.exp(2j * (tmax - tmin))))) / 2
base(ax, f"{next(num)}  Feret ellipse", f"widest and narrowest calliper widths;\nthe two widths are {ang_between:.0f}° apart, not 90°")
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
ax.add_patch(Ellipse(c, Fmax, Fmin, angle=np.degrees(tmax), fill=False, ec=EST, lw=2.2, ls="--",
                     label=shape_label(np.pi / 4 * Fmax * Fmin)))
finish(ax)

# minor-major axis ellipse (longest through the centre x perpendicular)
ax = fig.add_subplot(gs[0, next(panels)])
base(ax, f"{next(num)}  Minor–major axis ellipse", "longest diameter through the centre +\nthe diameter at 90° to it")
cd = ct_diams(xy)
dl = np.array([np.cos(cd["ct_theta"]), np.sin(cd["ct_theta"])]); ds = np.array([-dl[1], dl[0]])
pl, ql = line_through(c, dl); ps, qs = line_through(c, ds)
cal(ax, pl, ql, f"longest {cd['ct_long']:.1f} mm")
cal(ax, ps, qs, f"perpendicular {cd['ct_short']:.1f} mm", "--")
ax.add_patch(Ellipse(c, cd["ct_long"], cd["ct_short"], angle=np.degrees(cd["ct_theta"]), fill=False, ec=EST, lw=2.2,
                     ls="--", label=shape_label(np.pi / 4 * cd["ct_long"] * cd["ct_short"])))
finish(ax)

p1, q1 = line_through(c, u); p2, q2 = line_through(c, v)
Dp, Dq = np.linalg.norm(q1 - p1), np.linalg.norm(q2 - p2)
if SHOW_TWO_VIEW:  # two-view (PLAX x perpendicular)
    ax = fig.add_subplot(gs[0, next(panels)])
    off_deg = abs(np.degrees(np.angle(np.exp(2j * (np.arctan2(u[1], u[0]) - (th + np.pi / 2)))))) / 2
    base(ax, f"{next(num)}  Two-view ellipse", f"echo biplane style: PLAX diameter + the one\nat 90° to it; the PLAX line is {off_deg:.0f}° off the minor axis")
    cal(ax, p1, q1, f"PLAX diameter {Dp:.1f} mm")
    cal(ax, p2, q2, f"perpendicular {Dq:.1f} mm", "--")
    ax.add_patch(Ellipse(c, Dp, Dq, angle=np.degrees(np.arctan2(u[1], u[0])), fill=False, ec=EST, lw=2.2, ls="--",
                         label=shape_label(np.pi / 4 * Dp * Dq)))
    finish(ax)

# echo circle
ax = fig.add_subplot(gs[0, next(panels)])
base(ax, f"{next(num)}  Echocardiographic area", "approximation of the echo PLAX view")
cal(ax, p1, q1, f"PLAX diameter {Dp:.1f} mm")
ax.add_patch(Ellipse((p1 + q1) / 2, Dp, Dp, fill=True, fc=EST, alpha=.12, ec=EST, lw=2.2, ls="--",
                     label=shape_label(np.pi / 4 * Dp ** 2, "circle")))
finish(ax)

# bottom: population means, one panel per level
S = pd.read_csv(ROOT / "outputs/tables/area_methods_summary.csv")
M = [("ellipse_moment", "Moment\nellipse"), ("ellipse_feret", "Feret\nellipse"), ("ellipse_ct", "Minor–major\naxis ellipse")]
if SHOW_TWO_VIEW:
    M.append(("ellipse_2view", "Two-view\nellipse"))
M.append(("echo_circle", "Echocardiographic\narea"))
M = [(m, f"{i}\n{l}") for i, (m, l) in enumerate(M, 1)]
sub = gs[1, :].subgridspec(1, 3)
LEV = ["Aortic annulus", "LVOT 2 mm (closed)", "LVOT 4 mm (chord-closed)"]
for j, lev in enumerate(LEV):
    ax = fig.add_subplot(sub[0, j])
    x = np.arange(len(M))
    col = LEVEL_COL[lev]
    for k, sx in enumerate("FM"):
        row = S[(S.level == lev) & (S.sex == sx)].iloc[0]
        vals = [row[f"{m}_pct_mean"] for m, _ in M]
        sds = [row[f"{m}_pct_sd"] for m, _ in M]
        ax.bar(x + (k - .5) * .36, vals, width=.34, color=col, alpha={"F": .35, "M": .75}[sx], edgecolor=col,
               linewidth=1.2, yerr=sds, capsize=2, error_kw=dict(lw=.8, ecolor="0.4"),
               zorder=3, label={"F": "Women", "M": "Men"}[sx])
        for xi, val in zip(x + (k - .5) * .36, vals):
            ax.text(xi, val + (1.2 if val >= 0 else -1.2), f"{val:+.0f}", ha="center",
                    va="bottom" if val >= 0 else "top", fontsize=8, color="0.2")
    ax.axhline(0, color="0.3", lw=1)
    ax.set_xticks(x); ax.set_xticklabels([l for _, l in M], fontsize=9, fontweight="bold", color=NAME_COL)
    ax.set_ylim(-45, 15)
    if j == 0:
        ax.set_ylabel("% vs true area (mean ± SD)")
    ax.set_title(LEVEL_TITLE[lev], fontsize=13, fontweight="bold", loc="center", pad=12)
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)
    ax.legend(fontsize=8, frameon=False, loc="lower left")
    nct = S[(S.level == lev)].n_ct.sum(); nall = S[(S.level == lev)].n.sum()
    if nct < nall:
        ax.text(0.99, 0.03, f"CT ellipse: {nct} of {nall} meshes so far", transform=ax.transAxes, ha="right", fontsize=8, color="0.3")
fig.suptitle("Estimating the elliptical LVOT area", fontsize=16, fontweight="bold", y=0.998)
fig.tight_layout(rect=(0, 0, 1, 0.965), h_pad=4)
out = ROOT / "outputs/figures/analysis/fig_area_methods_demo.png"
fig.savefig(out, dpi=125)
print("wrote", out)
log_run(__file__, mesh=mid, output=str(out.relative_to(ROOT)), show_two_view=SHOW_TWO_VIEW)
