#!/usr/bin/env python3
"""
Build the LVOT CFD twins for one or more ImageCAS cases.

For each case:
  real   = the LVOT lumen (LV + aorta labels), cut in short-axis slices perpendicular to the LVOT axis
  circle = the same slices replaced by circles with the same area and the same centroid
Both twins share the centroid path, the area of every slice, the axial extent and the extensions,
so the only difference is the cross-section shape.

Domain along the LVOT axis (depth d, mm below the aortic valve plane; positive = into the LV):
  outlet ext  : the valve-plane slice repeated for OUT_EXT_D equivalent diameters (into the aorta)
  LVOT segment: d = START ... DEPTH mm (slices every DZ mm)
  inlet ext   : the deepest slice repeated for IN_EXT_D equivalent diameters (into the LV)
Flow goes from inlet (LV side) to outlet (aorta side). Leaflets are not modelled.

Outputs per case in OUT/{id}/:
  real.stl, circle.stl            ASCII STL with solids 'inlet', 'outlet', 'wall' (snappyHexMesh-ready), units mm
  slices.json                     per-slice area, ellipticity, AP/ML widths, echo-style area, centroid
  qc.png                          real contours on CT at several depths + both twins in 3D

Usage:
  python -I build_twins.py --ids 388,142 --out cfd/geometry
Labels: 1 myocardium, 2 LA, 3 LV, 6 aorta.
"""
import argparse, glob, json, os, sys
import numpy as np
import nibabel as nib
from scipy import ndimage as ndi
from skimage import measure

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
import lvot_curtain as lc  # reuse the locked LVOT frame (valve plane, axis, LA direction)

LAB_DIR = os.path.join(ROOT, 'labels_DTU/ImageCAS-STACOM2025-02-10-2025/segmentations')
IMG = {int(os.path.basename(p).split('.')[0]): p for p in glob.glob(os.path.join(ROOT, 'ImageCAS_*/*/*.img.nii.gz'))}
START = 2.0         # mm: first slice (the valve-plane slice cuts the sinuses and coronary ostia)
DEPTH = 12.0        # mm of LVOT below the valve plane (deeper slices open into the LV cavity)
DZ = 0.5            # mm between slices
SA_H = 0.2          # mm in-plane sampling of slices
SA_HALF = 30.0      # mm half-width of slice grid
NPTS = 128          # points per contour
MASK_SIGMA = 0.6    # mm in-plane smoothing of the binary slice before contouring
Z_SIGMA = 1.0       # mm smoothing of contours along the axis (removes voxel staircase)
IN_EXT_D, OUT_EXT_D = 2.0, 5.0   # extensions in equivalent diameters
EXT_DZ = 1.0        # mm between extension rings
ECHO_DEPTHS = (5.0, 7.5, 10.0)   # PW-Doppler / AP-diameter sample depths, mm
LEAFLET_SPLIT = True
CT_SIGMA = 0.5     # mm smoothing of CT before thresholding
THR_FRAC = 0.65    # threshold = myo + THR_FRAC * (blood - myo)
NRAY = 360         # rays per slice
RAY_H = 0.05       # mm step along rays
DARK_RUN = 0.4     # mm of consecutive dark samples that count as the leaflet
DARK_MIN_R = 3.0   # mm: ignore dark samples closer than this to the centre (noise)
MED_RAYS = 25      # rays (25 deg) in the circular median filter
EDGE_TOL = 0.6     # mm: a dark run closer than this to the label edge is the wall blur, not the leaflet


def _central(m):
    cc, nl = ndi.label(m)
    if nl == 0:
        return None
    cen = np.array(m.shape) // 2
    k = cc[tuple(cen)]
    if k == 0:
        pts_ = np.argwhere(cc > 0); k = cc[tuple(pts_[np.argmin(((pts_ - cen) ** 2).sum(1))])]
    return cc == k


def slice_mask(L, I, Ai, c, n, e1, e2, d, g):
    """Central LV+aorta label component in one short-axis slice, and the smoothed CT on the same grid."""
    U, V = np.meshgrid(g, g, indexing='ij')
    pts = (c - d * n) + U[..., None] * e1 + V[..., None] * e2
    ijk = (Ai[:3, :3] @ pts.reshape(-1, 3).T).T + Ai[:3, 3]
    lab = ndi.map_coordinates(L, ijk.T, order=0).reshape(U.shape)
    lab_m = ndi.binary_fill_holes(_central((lab == 3) | (lab == 6)))
    ct = ndi.gaussian_filter(ndi.map_coordinates(I, ijk.T, order=1).reshape(U.shape), CT_SIGMA / SA_H)
    return lab_m, ct


def ray_contour(lab_m, ct, thr, g, origin=None):
    """Star-shaped lumen contour from rays cast from the lumen centre.
    Each ray stops at the label edge or at the first dark run (CT < thr for >= DARK_RUN mm) inside the label:
    the anterior mitral leaflet / curtain, which the DTU LV label crosses. A circular median filter over angle
    removes rays that slip through the gap at the leaflet tip. Returns (NPTS, 2) points (u, v) in mm."""
    ii = np.argwhere(lab_m)
    org = g[0] + ii.mean(0) * SA_H if origin is None else origin
    th = np.linspace(0, 2 * np.pi, NRAY, endpoint=False)
    s = np.arange(0, SA_HALF, RAY_H)
    P = org[None, None, :] + s[None, :, None] * np.stack([np.cos(th), np.sin(th)], 1)[:, None, :]   # (NRAY, ns, 2)
    idx = (P - g[0]) / SA_H
    inl = ndi.map_coordinates(lab_m.astype(np.uint8), idx.reshape(-1, 2).T, order=0).reshape(P.shape[:2]) > 0
    val = ndi.map_coordinates(ct, idx.reshape(-1, 2).T, order=1).reshape(P.shape[:2])
    nrun = int(round(DARK_RUN / RAY_H))
    r = np.empty(NRAY)
    for k in range(NRAY):
        out_ = np.where(~inl[k])[0]; r_lab = s[out_[0]] if len(out_) else s[-1]
        dark = (val[k] < thr) & inl[k] & (s >= DARK_MIN_R)
        run = np.convolve(dark.astype(int), np.ones(nrun, int), 'valid') >= nrun
        hit = np.where(run)[0]
        r_dark = s[hit[0]] if len(hit) else r_lab
        # a dark run just inside the label edge is the blurred wall, not the leaflet: keep the label edge
        r[k] = r_lab if r_lab - r_dark < EDGE_TOL else r_dark
    r = ndi.median_filter(r, size=MED_RAYS, mode='wrap')
    r = ndi.gaussian_filter1d(r, 1.0, mode='wrap')
    thp = np.linspace(0, 2 * np.pi, NPTS, endpoint=False)
    rp = np.interp(thp, np.concatenate([th, [2 * np.pi]]), np.concatenate([r, r[:1]]))
    return org + rp[:, None] * np.stack([np.cos(thp), np.sin(thp)], 1)


def contour_of(mask, g):
    pad = 5
    sm = ndi.gaussian_filter(np.pad(mask, pad).astype(float), MASK_SIGMA / SA_H)
    cs = measure.find_contours(sm, 0.5)
    cs = [c for c in cs if np.allclose(c[0], c[-1])]
    c = max(cs, key=len)
    xy = g[0] + (c - pad) * SA_H               # (u, v) in mm
    return xy


def resample_polar(xy, cen):
    """Resample a closed contour to NPTS points at equal angles about cen; checks star-shape."""
    rel = xy - cen
    ang = np.unwrap(np.arctan2(rel[:, 1], rel[:, 0]))
    if ang[-1] < ang[0]:
        rel, ang = rel[::-1], ang[::-1]
    star = bool(np.all(np.diff(ang) > -1e-9))
    r = np.hypot(rel[:, 0], rel[:, 1])
    ang0 = ang - ang[0]
    th = np.linspace(0, 2 * np.pi, NPTS, endpoint=False)
    # equal-angle sampling (valid for star-shaped contours)
    a = np.mod(ang, 2 * np.pi)
    o = np.argsort(a)
    a_s, r_s = a[o], r[o]
    a_s = np.concatenate([a_s - 2 * np.pi, a_s, a_s + 2 * np.pi]); r_s = np.tile(r_s, 3)
    rr = np.interp(th, a_s, r_s)
    return np.stack([rr * np.cos(th), rr * np.sin(th)], 1) + cen, star


def poly_area(xy):
    x, y = xy[:, 0], xy[:, 1]
    return 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


def poly_centroid(xy):
    x, y = xy[:, 0], xy[:, 1]
    cr = x * np.roll(y, -1) - np.roll(x, -1) * y
    A = cr.sum() / 2
    return np.array([((x + np.roll(x, -1)) * cr).sum(), ((y + np.roll(y, -1)) * cr).sum()]) / (6 * A)


def shape_metrics(xy):
    """Area, ellipticity (sqrt eigen-ratio of the filled region), AP width (along u=e1) and ML width (along v=e2)
    through the centroid."""
    A = poly_area(xy); cen = poly_centroid(xy)
    # second moments of the polygon by rasterizing
    lo, hi = xy.min(0) - 1, xy.max(0) + 1
    h = 0.1
    gx, gy = np.meshgrid(np.arange(lo[0], hi[0], h), np.arange(lo[1], hi[1], h), indexing='ij')
    from matplotlib.path import Path
    inside = Path(xy).contains_points(np.stack([gx.ravel(), gy.ravel()], 1)).reshape(gx.shape)
    q = np.stack([gx[inside], gy[inside]], 1); q -= q.mean(0)
    ev = np.sort(np.linalg.eigvalsh(np.cov(q.T)))[::-1]
    # caliper widths through the centroid along e1 (AP) and e2 (ML)
    def width(axis):
        line = cen[None, :] + np.outer(np.arange(-30, 30, 0.05), axis)
        ins = Path(xy).contains_points(line)
        return float(ins.sum() * 0.05)
    ap, ml = width(np.array([1.0, 0.0])), width(np.array([0.0, 1.0]))
    return dict(area_mm2=float(A), ellipticity=float(np.sqrt(ev[0] / ev[1])), AP_mm=ap, ML_mm=ml,
                echo_area_mm2=float(np.pi * ap ** 2 / 4), cu=float(cen[0]), cv=float(cen[1]))


def to_world(rings_uv, depths, c, n, e1, e2):
    return np.array([[c - d * n + p[0] * e1 + p[1] * e2 for p in ring] for ring, d in zip(rings_uv, depths)])


def tri_surface(W):
    """W: (nrings, NPTS, 3) ordered from inlet to outlet. Returns dict of triangle lists per patch."""
    nr, npnt, _ = W.shape
    wall = []
    for i in range(nr - 1):
        for j in range(npnt):
            a, b = W[i, j], W[i, (j + 1) % npnt]
            c_, d_ = W[i + 1, j], W[i + 1, (j + 1) % npnt]
            wall += [(a, c_, b), (b, c_, d_)]
    def cap(ring, flip):
        cen = ring.mean(0)
        t = [(cen, ring[j], ring[(j + 1) % npnt]) for j in range(npnt)]
        return [(x, z, y) for x, y, z in t] if flip else t
    return dict(wall=wall, inlet=cap(W[0], False), outlet=cap(W[-1], True))


def orient_outward(patches, W):
    """Flip all triangles so normals point out of the lumen (checked with the signed volume)."""
    vol = 0.0
    for tris in patches.values():
        for a, b, c_ in tris:
            vol += np.dot(a, np.cross(b, c_)) / 6
    if vol < 0:
        patches = {k: [(a, c_, b) for a, b, c_ in v] for k, v in patches.items()}
    return patches, abs(vol)


def write_stl(path, patches):
    with open(path, 'w') as fh:
        for name, tris in patches.items():
            fh.write(f'solid {name}\n')
            for a, b, c_ in tris:
                nrm = np.cross(b - a, c_ - a); nn = np.linalg.norm(nrm); nrm = nrm / nn if nn > 0 else nrm
                fh.write(f' facet normal {nrm[0]:.6e} {nrm[1]:.6e} {nrm[2]:.6e}\n  outer loop\n')
                for p in (a, b, c_):
                    fh.write(f'   vertex {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}\n')
                fh.write('  endloop\n endfacet\n')
            fh.write(f'endsolid {name}\n')


def build(cid, out_root):
    out = os.path.join(out_root, str(cid)); os.makedirs(out, exist_ok=True)
    ln = nib.load(os.path.join(LAB_DIR, f'{cid}.img.nii.gz'))
    L = np.asarray(ln.dataobj).astype(np.uint8); A = ln.affine; Ai = np.linalg.inv(A)
    c, n, e1, e2 = lc.frame(L, A)
    I = np.asarray(nib.load(IMG[cid]).dataobj).astype(np.float32)
    bl, my = float(np.median(I[L == 3])), float(np.median(I[L == 1])); thr = my + THR_FRAC * (bl - my)
    g = np.arange(-SA_HALF, SA_HALF + SA_H / 2, SA_H)
    depths = np.arange(START, DEPTH + DZ / 2, DZ)
    rings, stars, edge, cut = [], [], [], []
    for d in depths:
        lab_m, ct = slice_mask(L, I, Ai, c, n, e1, e2, d, g)
        edge.append(bool(lab_m[0].any() or lab_m[-1].any() or lab_m[:, 0].any() or lab_m[:, -1].any()))
        if LEAFLET_SPLIT:
            ring = ray_contour(lab_m, ct, thr, g)
            ring = ray_contour(lab_m, ct, thr, g, origin=poly_centroid(ring))   # re-centre once
            star = True
        else:
            xy = contour_of(lab_m, g); ring, star = resample_polar(xy, poly_centroid(xy))
        cut.append(float(1 - poly_area(ring) / (lab_m.sum() * SA_H ** 2)))
        rings.append(ring); stars.append(star)
    R = np.array(rings)                                     # (nz, NPTS, 2)
    R = ndi.gaussian_filter1d(R, Z_SIGMA / DZ, axis=0, mode='nearest')   # smooth along the axis
    # per-slice metrics of the smoothed real contours; circles with the same area and centroid
    met = [shape_metrics(r) for r in R]
    th = np.linspace(0, 2 * np.pi, NPTS, endpoint=False)
    C = np.array([np.array([m_['cu'], m_['cv']]) + np.sqrt(m_['area_mm2'] / np.pi) * np.stack([np.cos(th), np.sin(th)], 1)
                  for m_ in met])
    deq0 = 2 * np.sqrt(met[0]['area_mm2'] / np.pi); deqN = 2 * np.sqrt(met[-1]['area_mm2'] / np.pi)
    d_out = START - np.arange(EXT_DZ, OUT_EXT_D * deq0 + EXT_DZ / 2, EXT_DZ)[::-1]   # negative depth = into aorta
    d_in = DEPTH + np.arange(EXT_DZ, IN_EXT_D * deqN + EXT_DZ / 2, EXT_DZ)
    all_d = np.concatenate([d_in[::-1], depths[::-1], d_out[::-1]])       # inlet (deep LV) -> outlet (aorta)
    res = {}
    for name, S in (('real', R), ('circle', C)):
        ring_seq = [S[-1]] * len(d_in) + list(S[::-1]) + [S[0]] * len(d_out)
        W = to_world(ring_seq, all_d, c, n, e1, e2)
        patches, vol = orient_outward(tri_surface(W), W)
        write_stl(os.path.join(out, f'{name}.stl'), patches)
        res[name] = dict(volume_mm3=float(vol))
    cm = [shape_metrics(r) for r in C]
    info = dict(id=cid, frame=dict(c=c.tolist(), n=n.tolist(), e1=e1.tolist(), e2=e2.tolist()),
                depth_mm=DEPTH, dz_mm=DZ, inlet_ext_mm=float(d_in[-1] - DEPTH), outlet_ext_mm=float(START - d_out[0]), start_mm=START,
                inlet_Deq_mm=float(deqN), outlet_Deq_mm=float(deq0), all_slices_star_shaped=bool(all(stars)), non_star_depths_mm=[float(d) for d, s in zip(depths, stars) if not s], slices_touching_grid_edge=int(sum(edge)),
                twins=res, leaflet_split=dict(on=LEAFLET_SPLIT, thr_HU=thr, blood_HU=bl, myo_HU=my), slices=[dict(depth_mm=float(d), label_area_cut_frac=cf, real=m_, circle=dict(area_mm2=cm_['area_mm2'], ellipticity=cm_['ellipticity']))
                                   for d, m_, cm_, cf in zip(depths, met, cm, cut)])
    json.dump(info, open(os.path.join(out, 'slices.json'), 'w'), indent=1)
    qc(cid, out, L, A, c, n, e1, e2, depths, R, C, all_d, d_in, d_out)
    return info


def qc(cid, out, L, A, c, n, e1, e2, depths, R, C, all_d, d_in, d_out):
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    I = np.asarray(nib.load(IMG[cid]).dataobj).astype(np.float32)
    Ai = np.linalg.inv(A); g = np.arange(-SA_HALF, SA_HALF + 0.01, 0.25); U, V = np.meshgrid(g, g, indexing='ij')
    bl, my = np.median(I[L == 3]), np.median(I[L == 1])
    show = [2, 4, 6, 8, 10, 12]
    f = plt.figure(figsize=(16, 7))
    for k, d in enumerate(show):
        ax = f.add_subplot(2, 6, k + 1)
        pts = (c - d * n) + U[..., None] * e1 + V[..., None] * e2
        ijk = (Ai[:3, :3] @ pts.reshape(-1, 3).T).T + Ai[:3, 3]
        im = ndi.map_coordinates(I, ijk.T, order=1).reshape(U.shape)
        lb = ndi.map_coordinates(L, ijk.T, order=0).reshape(U.shape)
        ax.imshow(im.T, cmap='gray', vmin=my - 150, vmax=bl + 150, origin='lower', extent=[g[0], g[-1], g[0], g[-1]])
        ax.contour(g, g, ((lb == 3) | (lb == 6)).T.astype(float), [0.5], colors='magenta', linewidths=0.6, linestyles='dotted')
        i = int(round((d - START) / DZ))
        for S, col in ((R, 'cyan'), (C, 'orange')):
            p = np.vstack([S[i], S[i][:1]]); ax.plot(p[:, 0], p[:, 1], color=col, lw=1)
        ax.set_title(f'{d:g} mm below valve', fontsize=9); ax.set_xticks([]); ax.set_yticks([])
        ax.set_xlim(-18, 18); ax.set_ylim(-18, 18)
    for k, (S, nm) in enumerate(((R, 'real'), (C, 'circle'))):
        ax = f.add_subplot(2, 2, 3 + k, projection='3d')
        seq = [S[-1]] * len(d_in) + list(S[::-1]) + [S[0]] * len(d_out)
        for ring, d in list(zip(seq, all_d))[::3]:
            p = np.vstack([ring, ring[:1]])
            ax.plot(p[:, 0], p[:, 1], -d, color='cyan' if nm == 'real' else 'orange', lw=0.5)
        ax.set_title(f'{nm} twin (u = toward LA/AP, v = ML, axis up = toward aorta)', fontsize=9)
        ax.set_box_aspect((1, 1, 2.5))
    f.suptitle(f'case {cid}: real lumen (cyan), same-area circle (orange), LV+aorta label (magenta dotted)', fontsize=11)
    plt.tight_layout(); plt.savefig(os.path.join(out, 'qc.png'), dpi=70); plt.close(f)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--ids', required=True); ap.add_argument('--out', default=os.path.join(HERE, 'geometry'))
    a = ap.parse_args()
    for cid in [int(x) for x in a.ids.split(',')]:
        info = build(cid, a.out)
        s = {x['depth_mm']: x['real'] for x in info['slices']}
        print(cid, 'star-shaped' if info['all_slices_star_shaped'] else f"NOT STAR-SHAPED at {info['non_star_depths_mm']}",
              ' '.join(f"d{d:g}: A={s[d]['area_mm2']:.0f} e={s[d]['ellipticity']:.2f} AP={s[d]['AP_mm']:.1f} ML={s[d]['ML_mm']:.1f}"
                       for d in (2.0, 7.5, 12.0)), '| edge slices:', info['slices_touching_grid_edge'],
              '| label cut at 5/10 mm: %.0f%% / %.0f%%' % tuple(100 * x['label_area_cut_frac'] for x in info['slices'] if x['depth_mm'] in (5.0, 10.0)), flush=True)
