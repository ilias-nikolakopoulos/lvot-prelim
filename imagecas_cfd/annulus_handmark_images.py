#!/usr/bin/env python3
"""
Images for hand marks of the aortic annulus level (option 5 of the valve-plane fix, 2026-11-16).

One image per scan: three long-axis planes, each through one sinus (found from the aorta label 3 mm above the
current 0 mm), with a depth scale in mm relative to the current 0 mm, and short-axis slices from -2 to +5 mm.
The images show the CT only. They do not show the labels or the automatic estimate, so that the marks are independent.

Usage (from this folder, without -I because it imports lvot_curtain):
  uv run --python 3.12 --with nibabel --with scipy --with numpy --with matplotlib --with pandas python annulus_handmark_images.py
Outputs in annulus_handmarks/: mark_{id}.png, handmarks_template.csv, subset.csv
"""
import glob, os
import numpy as np
import nibabel as nib
import pandas as pd
from scipy import ndimage as ndi
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import lvot_curtain as lc

LAB_DIR = 'labels_DTU/ImageCAS-STACOM2025-02-10-2025/segmentations'
OUT = 'annulus_handmarks'
SEED = 2026
N_RANDOM = 16
CFD_DIASTOLIC = [142, 167, 978, 587, 467, 663]
SINUS_LEVEL = -3.0      # mm (above the current 0 mm) where the three sinuses are located
SA_LEVELS = [-2, -1, 0, 1, 2, 3, 4, 5]
PLANE_COL = ['#ffd400', '#35d0e6', '#ff7ab6']


def subset():
    d = pd.read_csv('lvot_out/lvot_curtain_results.csv')
    sysf = sorted(d[d.systolic_flag].id.tolist())
    rest = d[~d.id.isin(sysf + CFD_DIASTOLIC)].id.to_numpy()
    rnd = sorted(np.random.default_rng(SEED).choice(rest, N_RANDOM, replace=False).tolist())
    rows = [dict(id=i, group='systolic-flagged') for i in sysf] + [dict(id=i, group='CFD diastolic') for i in CFD_DIASTOLIC] + [dict(id=i, group='random') for i in rnd]
    return pd.DataFrame(rows).sort_values('id').reset_index(drop=True)


def sinus_angles(lab_sec, g):
    """Angles (deg, in the e1-e2 plane) of the three sinuses: the three largest lobes of the aorta label section."""
    m = lab_sec == 6
    cc, nl = ndi.label(m)
    if nl == 0: return [90, 210, 330], False
    m = cc == (1 + int(np.argmax(ndi.sum(m, cc, range(1, nl + 1)))))
    p = np.argwhere(m) * (g[1] - g[0]) + g[0]; cen = p.mean(0); q = p - cen
    ang = np.degrees(np.arctan2(q[:, 1], q[:, 0])) % 360; r = np.hypot(q[:, 0], q[:, 1])
    prof = np.array([r[(ang >= a) & (ang < a + 5)].max() if ((ang >= a) & (ang < a + 5)).any() else 0 for a in range(0, 360, 5)])
    sm = ndi.gaussian_filter1d(prof, 2, mode='wrap')
    peaks = [i for i in range(72) if sm[i] >= sm[i - 1] and sm[i] > sm[(i + 1) % 72]]
    peaks.sort(key=lambda i: -sm[i]); chosen = []
    for i in peaks:
        if all(min(abs(i - j), 72 - abs(i - j)) * 5 >= 70 for j in chosen): chosen.append(i)
        if len(chosen) == 3: break
    if len(chosen) < 3: return [90, 210, 330], False
    return sorted(c * 5 + 2.5 for c in chosen), True


def make(cid, img):
    ln = nib.load(os.path.join(LAB_DIR, f'{cid}.img.nii.gz')); L = np.asarray(ln.dataobj).astype(np.uint8); A = ln.affine; Ai = np.linalg.inv(A)
    I = np.asarray(nib.load(img).dataobj).astype(np.float32)
    blood, myo = float(np.median(I[L == 3])), float(np.median(I[L == 1])); con = blood - myo
    c, n, e1, e2 = lc.frame(L, A)

    def samp(vol, o, a, b, ga, gb, order):
        U, V = np.meshgrid(ga, gb, indexing='ij'); pts = o[None, None, :] + U[..., None] * a + V[..., None] * b
        ijk = (Ai[:3, :3] @ pts.reshape(-1, 3).T).T + Ai[:3, 3]
        return ndi.map_coordinates(vol, ijk.T, order=order, cval=-1000 if order else 0).reshape(U.shape)
    kw = dict(cmap='gray', vmin=blood - 1.15 * con, vmax=blood + 0.9 * con, interpolation='bilinear')
    g = np.arange(-20, 20.01, 0.1)
    angs, found = sinus_angles(samp(L, c - SINUS_LEVEL * n, e1, e2, g, g, 0), g)
    fig = plt.figure(figsize=(21, 9.4)); gs = fig.add_gridspec(2, 24, height_ratios=[2.0, 1], hspace=0.12, wspace=0.9, left=0.035, right=0.995, top=0.9, bottom=0.01)
    gu = np.arange(-5, 17.01, 0.06); gv = np.arange(-8, 8.01, 0.06)   # the sinus side of each plane
    for k, (th, col, nm) in enumerate(zip(angs, PLANE_COL, 'ABC')):
        ax = fig.add_subplot(gs[0, 8 * k:8 * k + 8]); a = np.cos(np.radians(th)) * e1 + np.sin(np.radians(th)) * e2
        s = ndi.gaussian_filter(samp(I, c, a, -n, gu, gv, 1), 3.5)
        ax.imshow(s.T, extent=[gu[0], gu[-1], gv[-1], gv[0]], aspect='equal', **kw)
        for d in range(-8, 9):
            ax.axhline(d, color='white', lw=0.9 if d == 0 else 0.35, ls='-' if d == 0 else ':', alpha=0.9 if d == 0 else 0.55)
        ax.set_yticks(range(-8, 9)); ax.set_yticklabels([f'{d:+d}' if d else '0' for d in range(-8, 9)], fontsize=11); ax.set_xticks([0]); ax.set_xticklabels(['LVOT axis'], fontsize=9)
        ax.axvline(0, color='white', lw=0.4, ls=':', alpha=0.6)
        ax.set_title(f'Plane {nm}: long axis through sinus {nm}', fontsize=13, color='black', fontweight='bold')
        for sp in ax.spines.values(): sp.set_color(col); sp.set_linewidth(3)
        if k == 0: ax.set_ylabel('mm from the current 0 mm (+ = toward the LV)', fontsize=11)
    gsa = np.arange(-17, 17.01, 0.1)
    for k, d in enumerate(SA_LEVELS):
        ax = fig.add_subplot(gs[1, 3 * k:3 * k + 3]); s = ndi.gaussian_filter(samp(I, c - d * n, e1, e2, gsa, gsa, 1), 2)
        ax.imshow(s.T, extent=[gsa[0], gsa[-1], gsa[-1], gsa[0]], **kw); ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(f'short axis at {d:+d} mm' if d else 'short axis at 0 mm', fontsize=11)
        if d == SA_LEVELS[0]:
            for th, col, nm in zip(angs, PLANE_COL, 'ABC'):
                dx, dy = np.cos(np.radians(th)), np.sin(np.radians(th))
                ax.plot([-16 * dx, 16 * dx], [-16 * dy, 16 * dy], color=col, lw=1.2); ax.text(14.5 * dx, 14.5 * dy, nm, color=col, fontsize=12, fontweight='bold', ha='center', va='center')
    fig.suptitle(f'Case {cid}: find the lowest point of each cusp attachment (the sinus floor) and read its depth on the scale.'
                 + ('' if found else '   [sinuses not found automatically: planes at fixed angles]'), fontsize=14, x=0.03, ha='left', fontweight='bold')
    fig.savefig(os.path.join(OUT, f'mark_{cid}.png'), dpi=80); plt.close(fig)
    return found


if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    sub = subset(); sub.to_csv(os.path.join(OUT, 'subset.csv'), index=False)
    where = {int(os.path.basename(p).split('.')[0]): p for p in glob.glob('ImageCAS_*/*/*.img.nii.gz')}
    for cid in sub.id:
        print(cid, 'sinuses found' if make(int(cid), where[int(cid)]) else 'fixed planes', flush=True)   # one CT in memory at a time
    tp = os.path.join(OUT, 'handmarks_template.csv')
    if not os.path.exists(tp):   # never overwrite marks
        pd.DataFrame(dict(id=sub.id, annulus_depth_mm='', nadir_A_mm='', nadir_B_mm='', nadir_C_mm='', comment='')).to_csv(tp, index=False)
