#!/usr/bin/env python3
"""
LV-cavity version of the LVOT twins: the inflow comes from the patient's own LV cavity instead of a straight tube.

Domain (depth d along the LVOT axis, mm below the valve plane; positive = into the LV):
  d >= 16 mm        : the LV cavity (DTU label 3), identical for both twins
  12 < d < 16 mm    : smooth blend (signed-distance fields) from the twin's 12 mm slice into the LV cavity
  2 <= d <= 12 mm   : the twin's LVOT segment (real lumen or same-area circle), taken from {twin}.stl
  d_out <= d < 2 mm : the same aortic outlet extension as {twin}.stl
So the two twins differ only in the LVOT segment and the 4 mm blend; the LV and the outlet are shared.

Inflow: systolic ejection is driven by wall contraction. In this rigid model it is represented by a uniform
inflow through the endocardium: LV wall deeper than INLET_MIN_D and more than LA_GAP from the LA label = 'inlet'
(flowRateInletVelocity spreads Q uniformly over it). Everything else is 'wall', the aortic cap is 'outlet'.

Outputs in geometry/{id}/: {twin}_lv.stl (mm, solids inlet/outlet/wall), lv_info.json, qc_lv.png

Run from the imagecas folder:
  uv run --python 3.12 --with nibabel --with scipy --with numpy --with scikit-image --with matplotlib python -I cfd/build_lv_twins.py --ids 388,142
Enlarged-cavity sensitivity (systolic scans grown to the median diastolic cavity volume; LVOT unchanged):
  ... build_lv_twins.py --ids 388,576,33,586,253,146 --grow-to 88
"""
import argparse, json, os, sys
import numpy as np
import nibabel as nib
from scipy import ndimage as ndi
from skimage import measure
from matplotlib.path import Path

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from build_twins import LAB_DIR, write_stl  # noqa: E402

H = 0.4             # mm grid
BLEND = (12.0, 16.0)
INLET_MIN_D = 20.0  # mm
LA_GAP = 3.0        # mm
SMOOTH = 0.4        # mm Gaussian on the combined signed-distance field
GROW_TAPER = (16.0, 28.0)   # mm: LV enlargement (--grow-to) ramps in over this depth range
GROW_PAD = 15.0     # mm extra grid margin when growing


def read_rings(stl, c, n, e1, e2):
    """Recover the twin's rings (depth -> (NPTS,2) in u,v) from the wall solid of the existing STL."""
    V, inwall = [], False
    for line in open(stl):
        t = line.split()
        if not t: continue
        if t[0] == 'solid': inwall = t[1] == 'wall'
        elif t[0] == 'vertex' and inwall: V.append([float(x) for x in t[1:4]])
    rel = np.unique(np.round(np.array(V), 6), axis=0) - c
    d, u, v = -(rel @ n), rel @ e1, rel @ e2
    rings = {}
    for dd in np.unique(np.round(d, 3)):
        m = np.abs(d - dd) < 1e-3
        uv = np.stack([u[m], v[m]], 1); cen = uv.mean(0)
        rings[float(dd)] = uv[np.argsort(np.arctan2(uv[:, 1] - cen[1], uv[:, 0] - cen[0]))]
    return rings


def sdf2d(mask, h):
    return (ndi.distance_transform_edt(~mask) - ndi.distance_transform_edt(mask)) * h


def build(cid, twins=('real', 'circle'), apical=False, grow_to=None):
    gdir = os.path.join(HERE, 'geometry', str(cid))
    info = json.load(open(os.path.join(gdir, 'slices.json')))
    fr = {k: np.array(v) for k, v in info['frame'].items()}
    c, n, e1, e2 = fr['c'], fr['n'], fr['e1'], fr['e2']
    ln = nib.load(os.path.join(LAB_DIR, f'{cid}.img.nii.gz'))
    L = np.asarray(ln.dataobj).astype(np.uint8); A = ln.affine; Ai = np.linalg.inv(A)
    LV = (L == 3).astype(np.float32); LA = (L == 2).astype(np.float32)
    # LV voxels in frame coords -> grid extent
    ijk = np.argwhere(L == 3)[::7]
    xyz = ijk @ A[:3, :3].T + A[:3, 3]; rel = xyz - c
    dl, ul, vl = -(rel @ n), rel @ e1, rel @ e2
    d_out = info['start_mm'] - info['outlet_ext_mm']
    # inflow region: whole endocardium below INLET_MIN_D, or (apical sensitivity) only the deepest half of the LV
    d_inlet = INLET_MIN_D + 0.5 * (dl.max() - INLET_MIN_D) if apical else INLET_MIN_D
    sfx = '_lv_apical' if apical else ('_lv_big' if grow_to else '_lv')
    pad = GROW_PAD if grow_to else 0.0
    rings = {tw: read_rings(os.path.join(gdir, f'{tw}.stl'), c, n, e1, e2) for tw in twins}
    allr = np.vstack([r for tw in twins for r in rings[tw].values()])
    gd = np.arange(d_out - 2, dl.max() + 3 + pad, H)
    gu = np.arange(min(ul.min(), allr[:, 0].min()) - 3 - pad, max(ul.max(), allr[:, 0].max()) + 3 + pad, H)
    gv = np.arange(min(vl.min(), allr[:, 1].min()) - 3 - pad, max(vl.max(), allr[:, 1].max()) + 3 + pad, H)
    U, Vg = np.meshgrid(gu, gv, indexing='ij')
    print(f'{cid}: grid {len(gd)}x{len(gu)}x{len(gv)} = {len(gd)*len(gu)*len(gv)/1e6:.1f}M')
    # resample LV and LA labels slice by slice
    lv = np.zeros((len(gd), len(gu), len(gv)), bool); la = np.zeros_like(lv)
    for i, d in enumerate(gd):
        P = (c - d * n) + U[..., None] * e1 + Vg[..., None] * e2
        q = P.reshape(-1, 3) @ Ai[:3, :3].T + Ai[:3, 3]
        lv[i] = (ndi.map_coordinates(LV, q.T, order=1, cval=0) > 0.5).reshape(U.shape)
        la[i] = (ndi.map_coordinates(LA, q.T, order=1, cval=0) > 0.5).reshape(U.shape)
    # largest LV component
    cc, nl = ndi.label(lv); sz = np.bincount(cc.ravel()); sz[0] = 0; lv = cc == sz.argmax()
    phi_lv = (ndi.distance_transform_edt(~lv) - ndi.distance_transform_edt(lv)).astype(np.float32) * H
    la_dist = ndi.distance_transform_edt(~la).astype(np.float32) * H
    lv_label_ml = float(lv.sum() * H ** 3 / 1000); grow_mm = 0.0
    if grow_to:
        # enlarge the cavity to grow_to mL by an outward offset of its signed-distance field, tapered in from 0 at
        # GROW_TAPER[0] to the full offset at GROW_TAPER[1] mm below the valve, so the LVOT segment and blend stay identical
        t = np.clip((gd - GROW_TAPER[0]) / (GROW_TAPER[1] - GROW_TAPER[0]), 0, 1); t = (t * t * (3 - 2 * t))[:, None, None]
        lo_, hi_ = 0.0, GROW_PAD - 1.0
        for _ in range(30):
            mid = 0.5 * (lo_ + hi_)
            if ((phi_lv - mid * t) < 0).sum() * H ** 3 / 1000 < grow_to: lo_ = mid
            else: hi_ = mid
        grow_mm = 0.5 * (lo_ + hi_)
        phi_lv = (phi_lv - grow_mm * t).astype(np.float32); lv = phi_lv < 0
        print(f'  LV grown from {lv_label_ml:.1f} to {lv.sum() * H ** 3 / 1000:.1f} mL (offset {grow_mm:.2f} mm, tapered {GROW_TAPER[0]:g}-{GROW_TAPER[1]:g} mm)')
    # deep LV point for the pressure probe: farthest from the wall, deeper than 25 mm
    tmp = np.where((gd[:, None, None] > 25) & lv, -phi_lv, -1e9)
    ip = np.unravel_index(np.argmax(tmp), tmp.shape)
    p_lv = (c - gd[ip[0]] * n + gu[ip[1]] * e1 + gv[ip[2]] * e2)
    out = dict(id=cid, grid_mm=H, blend_mm=BLEND, inlet_min_d_mm=d_inlet, apical=apical, la_gap_mm=LA_GAP,
               p_probe_lv_world_mm=p_lv.tolist(), lv_volume_ml=float(lv.sum() * H ** 3 / 1000), lv_label_volume_ml=lv_label_ml,
               grow_to_ml=grow_to, grow_offset_mm=grow_mm, grow_taper_mm=GROW_TAPER if grow_to else None, twins={})
    keys = np.array(sorted(rings[twins[0]].keys()))
    pts2 = np.stack([U.ravel(), Vg.ravel()], 1)
    secs = {}
    for tw in twins:
        cache = {}
        def phi_twin(d):
            k = float(keys[np.argmin(np.abs(keys - d))])
            if k not in cache:
                m = Path(rings[tw][k]).contains_points(pts2).reshape(U.shape)
                cache[k] = sdf2d(m, H)
            return cache[k]
        phi = np.empty_like(phi_lv)
        for i, d in enumerate(gd):
            if d <= BLEND[0]:
                phi[i] = phi_twin(min(max(d, keys.min()), BLEND[0]))
            elif d < BLEND[1]:
                t = (d - BLEND[0]) / (BLEND[1] - BLEND[0]); t = t * t * (3 - 2 * t)
                phi[i] = (1 - t) * phi_twin(BLEND[0]) + t * phi_lv[i]
            else:
                phi[i] = phi_lv[i]
        phi = ndi.gaussian_filter(phi, SMOOTH / H)
        phi = np.maximum(phi, (d_out - gd)[:, None, None])           # flat aortic cap at d_out
        inside = phi < 0
        cc, nl = ndi.label(inside)
        i7 = np.argmin(np.abs(gd - 7.0)); s7 = info['slices'][10]['real']
        seed = cc[i7, np.argmin(np.abs(gu - s7['cu'])), np.argmin(np.abs(gv - s7['cv']))]
        assert seed > 0, 'LVOT seed outside the domain'
        phi = np.where((cc != seed) & inside, H, phi)
        secs[tw] = phi
        vtx, faces, _, _ = measure.marching_cubes(phi, 0.0, spacing=(H, H, H))
        dv, uv_, vv = vtx[:, 0] + gd[0], vtx[:, 1] + gu[0], vtx[:, 2] + gv[0]
        W = c - dv[:, None] * n + uv_[:, None] * e1 + vv[:, None] * e2
        tri = W[faces]
        # orient outward (positive signed volume)
        vol = np.einsum('ij,ij->i', tri[:, 0], np.cross(tri[:, 1], tri[:, 2])).sum() / 6
        if vol < 0: faces = faces[:, [0, 2, 1]]; tri = W[faces]
        fd = dv[faces].mean(1)
        nrm = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]); nrm /= np.linalg.norm(nrm, axis=1, keepdims=True) + 1e-12
        is_out = (fd < d_out + 0.5 * H) & (np.abs(nrm @ n) > 0.9)
        fi = np.clip(np.round((fd - gd[0]) / H).astype(int), 0, len(gd) - 1)
        fu = np.clip(np.round((uv_[faces].mean(1) - gu[0]) / H).astype(int), 0, len(gu) - 1)
        fv = np.clip(np.round((vv[faces].mean(1) - gv[0]) / H).astype(int), 0, len(gv) - 1)
        is_in = (~is_out) & (fd > d_inlet) & (la_dist[fi, fu, fv] > LA_GAP)
        area = 0.5 * np.linalg.norm(np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1)
        patches = {'inlet': [tuple(t) for t in tri[is_in]], 'outlet': [tuple(t) for t in tri[is_out]],
                   'wall': [tuple(t) for t in tri[~is_in & ~is_out]]}
        write_stl(os.path.join(gdir, f'{tw}{sfx}.stl'), patches)
        # check: section areas at the echo depths vs slices.json
        chk = {}
        for de in (5.0, 7.5, 10.0):
            i = np.argmin(np.abs(gd - de)); a = float((phi[i] < 0).sum() * H * H)
            ref = [s for s in info['slices'] if abs(s['depth_mm'] - de) < 1e-6][0][tw]['area_mm2']
            chk[str(de)] = dict(area_mm2=a, slices_json_mm2=ref, diff_pct=100 * (a / ref - 1))
        out['twins'][tw] = dict(volume_ml=float(abs(vol) / 1000), n_tri=int(len(faces)), inlet_area_mm2=float(area[is_in].sum()),
                                outlet_area_mm2=float(area[is_out].sum()), section_check=chk)
        print(f'  {tw}: {len(faces)} tri, vol {abs(vol)/1000:.1f} mL, inlet {area[is_in].sum():.0f} mm2, outlet {area[is_out].sum():.0f} mm2, '
              f'area check ' + ', '.join(f'{k}: {v["diff_pct"]:+.1f}%' for k, v in chk.items()))
    json.dump(out, open(os.path.join(gdir, f'lv_info{sfx[3:]}.json'), 'w'), indent=1)
    if not apical: qc(cid, gdir, gd, gu, gv, secs, lv, info, twins, name=f'qc{sfx}.png')


def qc(cid, gdir, gd, gu, gv, secs, lv, info, twins, name='qc_lv.png'):
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    s7 = info['slices'][10]['real']
    iu = np.argmin(np.abs(gu - s7['cu'])); iv = np.argmin(np.abs(gv - s7['cv']))
    fig, axs = plt.subplots(2, 3, figsize=(12, 8))
    for r, tw in enumerate(twins):
        phi = secs[tw]
        for k, (title, sl, ext, xl) in enumerate([
                ('long axis: d vs AP (through LVOT centroid)', phi[:, :, iv], [gu[0], gu[-1], gd[-1], gd[0]], 'AP (mm, + toward LA)'),
                ('long axis: d vs ML (through LVOT centroid)', phi[:, iu, :], [gv[0], gv[-1], gd[-1], gd[0]], 'ML (mm)')]):
            ax = axs[r, k]
            ax.imshow(sl < 0, cmap='Greys', extent=ext, aspect='equal', alpha=0.6)
            lvs = lv[:, :, iv] if k == 0 else lv[:, iu, :]
            ax.contour(np.linspace(ext[0], ext[1], lvs.shape[1]), np.linspace(ext[3], ext[2], lvs.shape[0]), lvs, [0.5], colors='tab:red', linewidths=0.7)
            for y in (2, 12, 16, INLET_MIN_D): ax.axhline(y, color='tab:blue' if y != INLET_MIN_D else 'tab:green', lw=0.6, ls='--')
            ax.set_ylim(min(gd[-1], 110), -40); ax.set_title(f'{tw}: {title}', fontsize=8); ax.set_xlabel(xl, fontsize=8); ax.set_ylabel('d (mm)', fontsize=8)
        ax = axs[r, 2]
        for de, col in ((5.0, 'tab:orange'), (10.0, 'tab:purple'), (14.0, 'tab:blue'), (20.0, 'tab:green')):
            i = np.argmin(np.abs(gd - de)); ax.contour(gv, gu, secs[tw][i], [0], colors=col, linewidths=1)
            ax.plot([], [], color=col, label=f'{de:g} mm')
        ax.set_aspect('equal'); ax.legend(fontsize=7); ax.set_title(f'{tw}: cross-sections (AP up, ML across)', fontsize=8)
    fig.suptitle(f'Case {cid}: LV-cavity twins. Grey = domain, red = LV label, blue dashed = 2/12/16 mm, green = inlet starts ({INLET_MIN_D:g} mm)', fontsize=9)
    fig.tight_layout(); fig.savefig(os.path.join(gdir, name), dpi=130); plt.close(fig)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--ids', default='388,142')
    ap.add_argument('--apical', action='store_true', help='inflow only through the deepest half of the LV (sensitivity)')
    ap.add_argument('--grow-to', type=float, default=None, help='enlarge the LV cavity to this volume (mL); writes {twin}_lv_big.stl, lv_info_big.json')
    a = ap.parse_args()
    for cid in [int(x) for x in a.ids.split(',')]:
        build(cid, apical=a.apical, grow_to=a.grow_to)
