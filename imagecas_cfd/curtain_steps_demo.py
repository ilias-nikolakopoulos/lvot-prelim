#!/usr/bin/env python3
"""
Step-by-step picture of the curtain rule in README_lvot_curtain.md, on one real case.
It uses the functions and the locked parameters of lvot_curtain.py, so the picture shows what the pipeline does.

Usage (from this folder, without -I because it imports lvot_curtain):
  uv run --python 3.12 --with nibabel --with scipy --with numpy --with matplotlib \
    python curtain_steps_demo.py --id 2
Output: curtain_steps_case{id}.png
"""
import argparse, glob, os
import numpy as np
import nibabel as nib
from scipy import ndimage as ndi
from scipy.signal import peak_prominences
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import lvot_curtain as lc

LAB_DIR = 'labels_DTU/ImageCAS-STACOM2025-02-10-2025/segmentations'
C_LA, C_LV, C_AO, C_WALL, C_END, C_REJ = '#ff4fd8', '#39d353', '#ff5a4f', '#35d0e6', '#ffd400', '#ff3b30'


def sample(vol, A, origin, a, b, ga, gb, order):
    """Resample vol on the plane origin + ga*a + gb*b."""
    U, V = np.meshgrid(ga, gb, indexing='ij')
    pts = origin[None, None, :] + U[..., None] * a + V[..., None] * b
    Ai = np.linalg.inv(A)
    ijk = (Ai[:3, :3] @ pts.reshape(-1, 3).T).T + Ai[:3, 3]
    return ndi.map_coordinates(vol, ijk.T, order=order, cval=-1000 if order else 0).reshape(U.shape)


def smooth_gap(r):
    """The smoothed gap that curtain_end_plane() uses (README 2.5)."""
    v, gap = r[:, 0], r[:, 3]
    good = np.isfinite(gap)
    return v, gap, ndi.uniform_filter1d(np.interp(v, v[good], gap[good]), lc.GAP_SMOOTH)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--id', type=int, default=2); a = ap.parse_args()
    cid = a.id
    img = [p for p in glob.glob(f'ImageCAS_*/*/{cid}.img.nii.gz')][0]
    ln = nib.load(os.path.join(LAB_DIR, f'{cid}.img.nii.gz')); L = np.asarray(ln.dataobj).astype(np.uint8); A = ln.affine
    I = np.asarray(nib.load(img).dataobj).astype(np.float32)
    blood, myo = float(np.median(I[L == 3])), float(np.median(I[L == 1]))
    vmin, vmax = myo - 150, blood + 150

    # README 2.1-2.7 with the pipeline's own functions
    c, n, e1, e2 = lc.frame(L, A)
    planes = {}
    for off in lc.OFFS:
        Ls = lc.longaxis(L, A, c, n, e1, e2, off, 0); s = lc.longaxis(I, A, c, n, e1, e2, off, 1)
        r = lc.trace(s, Ls); planes[off] = (s, Ls, r, lc.curtain_end_plane(r))
    ends = [planes[o][3] for o in lc.OFFS]
    ce, rule = lc.combine(ends)
    if ce is None: raise SystemExit(f'case {cid}: {rule}. Choose a case with a curtain end.')
    off = min((o for o in lc.OFFS if planes[o][3] is not None), key=lambda o: abs(planes[o][3] - ce))   # plane to show
    s, Ls, r, e = planes[off]
    ok = r[:, 4] == 1
    ext = [lc.GU[0], lc.GU[-1], lc.GV[-1], lc.GV[0]]
    xc = float(np.median(r[ok, 1])); X0, X1 = max(lc.GU[0], xc - 12), min(lc.GU[-1], xc + 7)   # view centred on the wall

    fig = plt.figure(figsize=(21, 11.5))
    gs = fig.add_gridspec(2, 4, hspace=0.3, wspace=0.2, left=0.035, right=0.985, top=0.9, bottom=0.06)
    T = dict(fontsize=12.5, loc='left', fontweight='bold'); BB = dict(facecolor='black', alpha=0.6, pad=3, edgecolor='none')
    note = dict(fontsize=9.5, color='white', va='top', bbox=BB)

    # ---- 2.1 reference: valve plane and LVOT axis (wide long-axis view, centre plane)
    ax = fig.add_subplot(gs[0, 0])
    gu, gv = np.arange(-32, 30.01, 0.25), np.arange(-22, 34.01, 0.25)
    w = sample(I, A, c, e1, -n, gu, gv, 1); wl = sample(L, A, c, e1, -n, gu, gv, 0)
    ax.imshow(w.T, cmap='gray', vmin=vmin, vmax=vmax, extent=[gu[0], gu[-1], gv[-1], gv[0]])
    for lab, col in ((3, C_LV), (6, C_AO), (2, C_LA)):
        ax.contour(gu, gv, (wl == lab).T.astype(float), [0.5], colors=col, linewidths=1.1)
    ax.axhline(0, color='white', lw=1.6); ax.annotate('', xy=(0, -11), xytext=(0, 30), arrowprops=dict(arrowstyle='->', color=C_END, lw=1.8))
    ax.text(-31, -1, 'valve plane = 0 mm', color='white', fontsize=10, va='bottom', bbox=BB)
    ax.text(1.2, 28, 'LVOT axis', color=C_END, fontsize=10, bbox=BB)
    for txt, x, y, col in (('aorta', -14, -6, C_AO), ('LV', -16, 24, C_LV), ('LA', 19, 12, C_LA)):
        ax.text(x, y, txt, color=col, fontsize=12, fontweight='bold', bbox=BB)
    ax.add_patch(Rectangle((lc.GU[0], lc.GV[0]), lc.GU[-1] - lc.GU[0], lc.GV[-1] - lc.GV[0], fill=False, ec='white', ls=':', lw=1))
    ax.set_title('2.1  Set the reference', **T); ax.set_ylabel('mm below the valve plane'); ax.set_xlabel('mm toward the LA')
    ax.text(0.02, 0.985, 'The LV-aorta label boundary is the valve plane.\nThe axis is perpendicular to it.\nDotted box = the area of the next panels.',
            transform=ax.transAxes, **note)

    # ---- 2.2 five long-axis planes, seen from above (short-axis view 4 mm below the valve plane)
    ax = fig.add_subplot(gs[0, 1])
    g = np.arange(-24, 24.01, 0.25); d = 4.0
    sa = sample(I, A, c - d * n, e1, e2, g, g, 1); sl = sample(L, A, c - d * n, e1, e2, g, g, 0)
    ax.imshow(sa.T, cmap='gray', vmin=vmin, vmax=vmax, extent=[g[0], g[-1], g[-1], g[0]])
    for lab, col in ((3, C_LV), (2, C_LA)):
        ax.contour(g, g, (sl == lab).T.astype(float), [0.5], colors=col, linewidths=1.1)
    for o in lc.OFFS:
        ax.plot([lc.GU[0], lc.GU[-1]], [o, o], color=C_END if o == off else 'white', lw=2.2 if o == off else 1.1)
        ax.text(lc.GU[0] - 0.8, o, f'{o:+d}', color=C_END if o == off else 'white', fontsize=8, ha='right', va='center', bbox=dict(facecolor='black', alpha=0.6, pad=0.6, edgecolor='none'))
    ax.plot(0, 0, '+', color=C_END, ms=12, mew=2)
    ax.annotate('to the LA', xy=(23, 21), xytext=(8, 21), color=C_LA, fontsize=10, va='center', bbox=BB, arrowprops=dict(arrowstyle='->', color=C_LA, lw=1.5))
    ax.set_title('2.2  Make five long-axis planes', **T); ax.set_xlabel('mm toward the LA'); ax.set_ylabel('mm across the LVOT')
    ax.text(0.02, 0.985, f'View down the LVOT axis (+), {d:g} mm below the valve.\nEach line is one plane. The planes are 2 mm apart.\nYellow = plane {off:+d}, shown in the next panels.',
            transform=ax.transAxes, **note)

    # ---- 2.3 trace the posterior wall in one plane
    ax = fig.add_subplot(gs[0, 2])
    ax.imshow(s.T, cmap='gray', vmin=vmin, vmax=vmax, extent=ext)
    ax.contour(lc.GU, lc.GV, (Ls == 2).T.astype(float), [0.5], colors=C_LA, linewidths=1.0)
    la = r[:, 2]; fin = np.isfinite(la)
    ax.fill_betweenx(r[fin, 0], la[fin] - lc.WIN_FRONT * lc.H, la[fin] + lc.WIN_BACK * lc.H, color='white', alpha=0.16, lw=0)
    ax.plot(la[fin], r[fin, 0], color=C_LA, lw=2)
    ax.plot(r[ok, 1], r[ok, 0], '.', color=C_WALL, ms=3)
    ax.axhline(0, color='white', lw=1.2); ax.axhline(lc.START, color='white', lw=0.9, ls='--')
    ax.text(X0 + 0.3, lc.START - 0.2, 'start: 2 mm', color='white', fontsize=9, va='bottom', bbox=BB)
    ax.set_xlim(X0, X1); ax.set_ylim(15, -5)
    ax.set_title(f'2.3  Trace the posterior wall (plane {off:+d})', **T); ax.set_xlabel('mm toward the LA')
    ax.text(0.02, 0.985, 'Pink = front edge of the LA label.\nLight band = search area (6 mm in front, 1.5 mm behind).\nBlue = darkest point in the band = the wall.',
            transform=ax.transAxes, **note)

    # ---- 2.4 safeguards, on one row of the image
    ax = fig.add_subplot(gs[0, 3])
    sm = ndi.gaussian_filter(s, lc.SMOOTH_SIGMA)
    k = int(np.argmin(np.abs(r[:, 0] - 4.0))); vrow = r[k, 0]; j = int(np.argmin(np.abs(lc.GV - vrow)))
    x = int(np.argmin(np.abs(lc.GU - r[k, 1]))); xla = r[k, 2]
    prom, lb, rb = peak_prominences(-sm[:, j], [x], wlen=lc.PROM_WLEN)
    base = max(sm[lb[0], j], sm[rb[0], j]) if sm[lb[0], j] > sm[x, j] and sm[rb[0], j] > sm[x, j] else sm[x, j] + prom[0]
    ax.plot(lc.GU, sm[:, j], color='#333', lw=1.6)
    ax.axvspan(xla - lc.WIN_FRONT * lc.H, xla + lc.WIN_BACK * lc.H, color='#bbb', alpha=0.35, lw=0, label='search area')
    ax.axvline(xla, color=C_LA, lw=2, label='LA label edge')
    ax.plot(lc.GU[x], sm[x, j], 'o', color='#0a8fa3', ms=9, label='darkest point = wall')
    ax.annotate('', xy=(lc.GU[x], sm[x, j]), xytext=(lc.GU[x], sm[x, j] + prom[0]), arrowprops=dict(arrowstyle='<->', color='#b45309', lw=1.6))
    ax.text(lc.GU[x] - 0.4, sm[x, j] + prom[0] / 2, f'contrast {prom[0]:.0f} HU\n(must be 30 or more)', color='#b45309', fontsize=9.5, va='center', ha='right')
    ax.axvspan(lc.GU[x] - lc.STEP * lc.H, lc.GU[x] + lc.STEP * lc.H, color='#0a8fa3', alpha=0.25, lw=0, label='0.3 mm limit for the next step')
    ax.set_xlim(X0, X1); ax.set_xlabel('mm toward the LA'); ax.set_ylabel('CT number (HU), smoothed')
    seg = sm[(lc.GU >= X0) & (lc.GU <= X1), j]; ax.set_ylim(seg.min() - 0.42 * np.ptp(seg), seg.max() + 0.42 * np.ptp(seg))
    ax.legend(fontsize=8.5, frameon=False, loc='lower left', ncol=2)
    ax.set_title(f'2.4  Apply the safeguards (row at {vrow:.1f} mm)', **T)
    ax.text(0.02, 0.985, 'A point is rejected if it is not dark enough,\nor if it is too far to the side.\nThe trace stops after 1 mm of rejected points.',
            transform=ax.transAxes, fontsize=9.5, va='top')
    for sp in ('top', 'right'): ax.spines[sp].set_visible(False)

    # ---- 2.5 the gap, drawn on the image
    ax = fig.add_subplot(gs[1, 0])
    ax.imshow(s.T, cmap='gray', vmin=vmin, vmax=vmax, extent=ext)
    ax.plot(la[fin], r[fin, 0], color=C_LA, lw=2); ax.plot(r[ok, 1], r[ok, 0], '.', color=C_WALL, ms=3)
    for i in np.where(ok)[0][::5]:
        ax.plot([r[i, 1], r[i, 2]], [r[i, 0], r[i, 0]], color=C_END, lw=1.3)
    if (~ok).any(): ax.plot(r[~ok, 1], r[~ok, 0], 'x', color=C_REJ, ms=5, label='rejected point')
    ax.axhline(e, color=C_END, ls='--', lw=1.4)
    ax.text(X0 + 0.4, e - 0.25, f'gap starts to open: {e:.1f} mm', color=C_END, fontsize=10, va='bottom', bbox=BB)
    ax.text(X0 + 0.4, (lc.START + e) / 2, 'curtain', color='white', fontsize=11, va='center', fontweight='bold', bbox=BB)
    ax.text(X0 + 0.4, min(e + 3.5, 13.5), 'anterior\nmitral leaflet', color='white', fontsize=11, va='center', fontweight='bold', bbox=BB)
    ax.set_xlim(X0, X1); ax.set_ylim(15, -5)
    ax.set_title('2.5  Measure the gap', **T); ax.set_xlabel('mm toward the LA'); ax.set_ylabel('mm below the valve plane')
    ax.text(0.02, 0.985, 'Yellow lines = gap from the wall (blue) to the LA edge (pink).\nAbove the hinge the wall is against the LA.\nBelow it, the open leaflet moves away from the LA.',
            transform=ax.transAxes, **note)

    # ---- 2.6 curtain end in this plane
    ax = fig.add_subplot(gs[1, 1])
    v, gap, gsm = smooth_gap(r)
    ax.plot(v[ok], gap[ok], '.', color='#9aa', ms=4, label='gap at each 0.1 mm step')
    ax.plot(v, gsm, color='#0a8fa3', lw=2.2, label='gap smoothed over 1 mm')
    ax.axhline(lc.GAP_THR, color='#b45309', ls='--', lw=1.2); ax.text(v[-1], lc.GAP_THR, ' 0.2 mm', color='#b45309', fontsize=9, va='center')
    ax.axvline(e, color='#c79a00', lw=2); ax.axvspan(e, e + lc.LOOK * lc.H, color='#ffd400', alpha=0.22, lw=0, label='next 3 mm: the check')
    i0 = int(np.argmin(np.abs(v - e))); i1 = i0 + lc.LOOK
    ax.annotate('', xy=(v[i1], gsm[i1]), xytext=(v[i1], gsm[i0]), arrowprops=dict(arrowstyle='<->', color='#b45309', lw=1.6))
    ax.text(v[i1] + 0.2, gsm[i0] - 0.12, f'rise {gsm[i1] - gsm[i0]:.1f} mm in 3 mm\n(must be more than 0.3)', color='#b45309', fontsize=9.5, va='top')
    ax.text(e + 0.1, 0.62 * ax.get_ylim()[1], f'curtain end\n{e:.1f} mm', color='#8a6a00', fontsize=10.5, va='center', fontweight='bold')
    ax.set_xlim(lc.START, min(v[-1] + 2.2, 16)); ax.set_xlabel('mm below the valve plane'); ax.set_ylabel('gap, wall to LA edge (mm)')
    ax.legend(fontsize=8.5, frameon=False, loc='upper left')
    ax.set_title(f'2.6  Find the curtain end in one plane (plane {off:+d})', **T)
    for sp in ('top', 'right'): ax.spines[sp].set_visible(False)

    # ---- 2.7 combine the five planes
    ax = fig.add_subplot(gs[1, 2])
    for o, en in zip(lc.OFFS, ends):
        if en is None: ax.text(o, ce, 'no end', ha='center', va='center', fontsize=9, color='#888', rotation=90)
        else:
            ax.plot(o, en, 'o', color='#c79a00' if o == off else '#555', ms=12)
            ax.text(o, en - 0.28, f'{en:.1f}', ha='center', va='bottom', fontsize=10)
    ax.axhline(ce, color='#2a9d3a', lw=2.2, ls=':')
    ax.text(7.2, ce + 0.15, f'case result:\n{ce:.1f} mm', color='#2a9d3a', fontsize=11, va='top', ha='right', fontweight='bold')
    found = [x for x in ends if x is not None]
    ax.set_xlim(-5, 7.3); ax.set_ylim(max(found) + 1.6, min(found) - 1.6); ax.set_xticks(lc.OFFS); ax.set_xticklabels([f'{o:+d}' for o in lc.OFFS])
    ax.set_xlabel('plane (mm across the LVOT)'); ax.set_ylabel('curtain end (mm below the valve plane)')
    ax.set_title('2.7  Find the curtain end for the case', **T)
    ax.text(0.03, 0.03, '3 or more planes: use the median.\n2 planes within 1.5 mm: use the mean.\nIn all other conditions: "not found".\n'
            f'This case: {rule}.', transform=ax.transAxes, fontsize=9.5, va='bottom')
    for sp in ('top', 'right'): ax.spines[sp].set_visible(False)

    # ---- 3 short-axis measurements at 2 mm, mid and curtain end
    sub = gs[1, 3].subgridspec(1, 3, wspace=0.06)
    g = np.arange(-22, 22.01, lc.SA_H)
    for kk, (name, d) in enumerate((('2 mm', 2.0), ('mid', ce / 2), ('curtain end', ce))):
        ax = fig.add_subplot(sub[0, kk])
        sa = sample(I, A, c - d * n, e1, e2, g, g, 1); m = lc.shortaxis(L == 3, A, c, n, e1, e2, d)
        lv = sample((L == 3).astype(np.uint8), A, c - d * n, e1, e2, g, g, 0)
        lab, _ = ndi.label(lv); cen = lab[len(g) // 2, len(g) // 2]
        ax.imshow(sa.T, cmap='gray', vmin=vmin, vmax=vmax, extent=[g[0], g[-1], g[-1], g[0]])
        if cen: ax.contour(g, g, (lab == cen).T.astype(float), [0.5], colors=C_END, linewidths=1.6)
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_xlabel(f'{name} ({d:.1f} mm)\narea {m["area_mm2"]:.0f} mm²\nellipticity {m["ellipticity"]:.2f}\n{m["major_mm"]:.1f} × {m["minor_mm"]:.1f} mm', fontsize=9.5)
        if kk == 0: ax.set_title('3  Measure the cross-section', **T)
        if kk == 2: ax.set_title('yellow = LV label', fontsize=9, loc='right')

    fig.suptitle(f'How the rule finds the end of the aorto-mitral curtain, step by step (ImageCAS case {cid}; step numbers as in README_lvot_curtain.md)',
                 fontsize=15, x=0.035, ha='left', fontweight='bold')
    out = f'curtain_steps_case{cid}.png'
    fig.savefig(out, dpi=110); print(out, 'curtain end', ce, rule, 'plane shown', off, 'plane ends', ends)


if __name__ == '__main__':
    main()
