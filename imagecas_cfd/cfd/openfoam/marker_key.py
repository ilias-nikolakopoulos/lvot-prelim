"""Explainer figure for the symbols used in the LVOT cross-section figures (plane_figures.py).

Run from the imagecas folder:
  uv run --python 3.12 --with numpy --with scipy --with pandas --with matplotlib python -I cfd/openfoam/marker_key.py --tag lv_lam_h1 --case 388_real --depth 7.5
"""
import argparse, json, os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from plane_figures import section, draw, CFD  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument('--tag', default='lv_lam_h1'); ap.add_argument('--case', default='388_real'); ap.add_argument('--depth', type=float, default=7.5)
args = ap.parse_args()
case = os.path.join(CFD, 'runs', args.tag, args.case)
cj = json.load(open(os.path.join(case, 'case.json')))
info = json.load(open(os.path.join(CFD, 'geometry', str(cj['id']), 'slices.json')))
S = section(case, args.depth, info, cj)
r = S['row']

fig = plt.figure(figsize=(15, 8.2))
ax = fig.add_axes([0.08, 0.17, 0.36, 0.7])
im = draw(ax, S, 1.6)
cb = fig.colorbar(im, cax=fig.add_axes([0.1, 0.07, 0.32, 0.025]), orientation='horizontal'); cb.set_label('w / v$_{mean}$  (local axial velocity ÷ true mean velocity)')
ax.set_title(f'Case {cj["id"]} {cj["twin"]} LVOT, {args.depth:g} mm below the valve ({args.tag})', fontsize=10, pad=22)

# callouts
box = dict(boxstyle='round,pad=0.3', fc='white', ec='0.4', alpha=0.92)
arr = dict(arrowstyle='->', color='white', lw=1.2)
umin, umax = np.nanmin(S['u']), np.nanmax(S['u']); vmin, vmax = np.nanmin(S['v']), np.nanmax(S['v'])
ax.annotate('1  + centroid\n    = PW gate centre', (S['cv'], S['cu']), xytext=(vmax + 1.5, S['cu'] + 5), fontsize=8.5, bbox=box, arrowprops=arr, annotation_clip=False)
ax.annotate('2  dashed circle\n    = PW sample volume', (S['cv'] + 1.1, S['cu'] - 1.1), xytext=(vmax + 1.5, S['cu'] - 5), fontsize=8.5, bbox=box, arrowprops=arr, annotation_clip=False)
ax.annotate('3  x = velocity maximum', (S['pv'], S['pu']), xytext=(vmin - 8, umin + 1), fontsize=8.5, bbox=box, arrowprops=dict(arr, color='#ff5a5f'), annotation_clip=False, ha='left')
ax.annotate('4  o = flow centre', (S['fv'], S['fu']), xytext=(vmin - 8, umax - 1), fontsize=8.5, bbox=box, arrowprops=dict(arr, color='#ffb000'), annotation_clip=False, ha='left')
# 5: the mean-velocity line (w = v_mean): where a gate would read the true mean velocity (VTI-correct)
from plane_figures import MEANC  # noqa: E402
cs = ax.contour(S['v'], S['u'], np.nan_to_num(S['W'], nan=-1), levels=[1.0], colors=MEANC, linewidths=2.6)
seg = max((p for lev in cs.allsegs for p in lev), key=len)
k = np.argmax(seg[:, 0] + 0.3 * seg[:, 1])          # a point on the right-hand (mitral/LA-ish) part of the line
ax.annotate('5  magenta line = mean velocity\n    (w = v$_{mean}$, the VTI-correct speed)', tuple(seg[k]), xytext=(vmax - 9, umin - 3.5),
            fontsize=8.5, bbox=box, arrowprops=dict(arr, color=MEANC), annotation_clip=False)
ax.text(0.5, 1.0, '↑ toward LA / anterior mitral leaflet', transform=ax.transAxes, ha='center', va='bottom', fontsize=8.5, color='0.25')
ax.text(0.02, -0.02, '↓ septal side', transform=ax.transAxes, ha='left', va='top', fontsize=8.5, color='0.25')

# key panel
kx = fig.add_axes([0.56, 0.04, 0.44, 0.9]); kx.axis('off'); kx.set_xlim(0, 1); kx.set_ylim(0, 1)
entries = [
    (dict(marker='+', color='k', ms=13, mew=2, ls='none'), '1  +  Slice centroid (white cross)',
     'The geometric centre of the LVOT cross-section. This is where we place the PW-Doppler gate,\n'
     'as a sonographer would aim for the middle of the LVOT.'),
    (None, '2  - -  PW sample volume (dashed white circle)',
     'The 1.5 mm-radius gate around the centroid. The echo velocity in our results is the average\n'
     'over this region (centre + 6 points on the circle).'),
    (dict(marker='x', color='#ff5a5f', ms=12, mew=2.2, ls='none'), '3  x  Velocity maximum (red cross)',
     f'The single fastest point in the slice (after 0.3 mm smoothing). Here {r["peak_over_mean"]:.2f}× the mean,\n'
     f'{r["disp_mm"]:.1f} mm from the centroid. In a flat profile this can jump around because many points are\n'
     'almost equally fast, so it is a weak measure of skew on its own.'),
    (dict(marker='o', mfc='none', mec='#ffb000', ms=12, mew=2.2, ls='none'), '4  o  Flow centre (yellow/orange circle)',
     'The "centre of mass" of the flow: every point weighted by how fast blood goes through it.\n'
     f'If the profile is symmetric it sits on the geometric centre. Here it is {r["flow_centre_off_mm"]:.1f} mm off\n'
     f'(normalised flow displacement {r["flow_displacement"]:.3f}, the 4D-flow MRI measure). This is the robust\n'
     'measure of where the flow is shifted.'),
    ('MEAN', '5  Mean-velocity line (magenta)',
     'Every point on this line moves at exactly the true cross-sectional mean velocity (w = v_mean).\n'
     'This is the speed SV = CSA x VTI needs. Inside it (the whole core, including + and o) blood is\n'
     f'faster (here the gate reads {S["W"][np.argmin(np.abs(S["u"] - S["cu"])), np.argmin(np.abs(S["v"] - S["cv"]))]:.2f}x the mean); outside it, in the thin band near the wall,\n'
     'slower. Mostly it hugs the wall, where velocity changes steeply; on the side where flow is slower\n'
     '(here the right-hand, mediolateral side) it sits a few mm in. Not a practical, reproducible gate site.'),
    (None, 'Colour and contours',
     'Colour = w / v_mean: 1 means exactly the true mean velocity (a gate there reads the right velocity),\n'
     '>1 faster than the mean (gate overestimates), <1 slower (dark near the walls).\n'
     'Thin white lines = contours at 0.25, 0.5, 0.75 and 1.25 × mean.'),
]
y = 0.97
for mk, head, body in entries:
    if mk == 'MEAN':
        kx.add_line(Line2D([0.0, 0.04], [y - 0.012] * 2, transform=kx.transAxes, color=MEANC, lw=3))
    elif mk:
        kx.add_line(Line2D([0.02], [y - 0.012], transform=kx.transAxes, **mk))
    elif 'dashed' in head:
        kx.add_patch(plt.Circle((0.02, y - 0.012), 0.012, transform=kx.transAxes, fill=False, ls='--', color='k', lw=1.4))
    kx.text(0.06, y, head, fontsize=10.5, weight='bold', va='top')
    kx.text(0.06, y - 0.04, body, fontsize=8.6, va='top', linespacing=1.35)
    y -= 0.06 + 0.036 * (body.count('\n') + 1)
kx.text(0.0, y + 0.01, 'Echo reads the velocity at the + / dashed circle. The difference between that and the true mean\n'
        '(magenta line) is the "profile" part of the echo SV error; the AP-diameter area is the other part.',
        fontsize=9, style='italic', va='top', color='0.25')
out = os.path.join(CFD, 'qc', f'marker_key_{args.tag}_{args.case}.png')
fig.savefig(out, dpi=150, bbox_inches='tight'); plt.close(fig)
print('wrote', os.path.relpath(out, CFD))
