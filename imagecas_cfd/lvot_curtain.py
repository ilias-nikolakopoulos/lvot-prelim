#!/usr/bin/env python3
"""
LVOT measurement on ImageCAS CCTA with the DTU whole-heart labels.
The LVOT ends where the aorto-mitral curtain ends (anterior mitral leaflet hinge).
Rule agreed 2026-11-11; see README_lvot_curtain.md.

Usage:
  python3 -I lvot_curtain.py --img-dir IMG_DIR [IMG_DIR2 ...] --lab-dir LAB_DIR --out OUT_DIR [--ids 1,2,3] [--workers 4]
  IMG_DIR holds {id}.img.nii.gz (ImageCAS CT); LAB_DIR holds {id}.img.nii.gz (DTU labels).
  The script skips cases that already have OUT_DIR/cases/{id}.json, so you can stop and restart it.
Labels: 1 myocardium, 2 LA, 3 LV, 6 aorta.
"""
import argparse, glob, json, os, sys, traceback
import numpy as np
import warnings; warnings.filterwarnings("ignore")
import nibabel as nib
from scipy import ndimage as ndi
from scipy.signal import peak_prominences

# ---------------- parameters (locked) ----------------
H = 0.1                      # in-plane sampling of long-axis planes, mm
GU = np.arange(-8, 18 + H / 2, H)      # mm along e1 (toward LA)
GV = np.arange(-8, 16 + H / 2, H)      # mm below valve plane
OFFS = [-4, -2, 0, 2, 4]     # parallel long-axis planes, mm
SMOOTH_SIGMA = 3             # px (0.3 mm) image smoothing for the trace
WIN_FRONT, WIN_BACK = 60, 15 # px: search 6 mm in front of LA edge, 1.5 mm behind
START = 2.0                  # mm: start trace 2 mm below valve plane
STEP = 3                     # px: max sideways move per 0.1 mm step (0.3 mm)
PROM = 30                    # HU: min dip contrast within +-1.5 mm
PROM_WLEN = 31               # px window for contrast
MAXMISS = 10                 # rows: stop after 1 mm of rejected points
GAP_SMOOTH = 10              # rows: 1 mm smoothing of gap
GAP_THR = 0.2                # mm
LOOK = 30                    # rows: 3 mm confirmation window
RISE = 0.3                   # mm: gap must rise by more than this within LOOK
AGREE2 = 1.5                 # mm: two planes must agree within this
SA_H = 0.25                  # short-axis grid, mm
SA_HALF = 30                 # short-axis half-width, mm
SYSTOLE_RATIO = 0.6          # LV/myocardium volume ratio below this -> flag as systolic


def frame(L, A):
    """Valve plane centre c, LVOT axis n (into aorta), e1 toward LA, e2."""
    LV, AO, LA = L == 3, L == 6, L == 2
    itf = np.argwhere(LV & ndi.binary_dilation(AO, structure=ndi.generate_binary_structure(3, 1)))
    if len(itf) < 20:
        raise RuntimeError('no LV-aorta interface')
    P = (A[:3, :3] @ itf.T).T + A[:3, 3]
    c = P.mean(0)
    w, v = np.linalg.eigh(np.cov((P - c).T))
    n = v[:, 0]
    aop = (A[:3, :3] @ np.argwhere(AO).T).T + A[:3, 3]
    near = aop[np.linalg.norm(aop - c, axis=1) < 15]
    if (near.mean(0) - c) @ n < 0:
        n = -n
    lap = (A[:3, :3] @ np.argwhere(LA)[::50].T).T + A[:3, 3]
    t = lap.mean(0) - c
    t -= (t @ n) * n
    e1 = t / np.linalg.norm(t)
    e2 = np.cross(n, e1)
    return c, n, e1, e2


def longaxis(vol, A, c, n, e1, e2, off, order):
    U, V = np.meshgrid(GU, GV, indexing='ij')
    Ai = np.linalg.inv(A)
    pts = (c + off * e2)[None, None, :] + U[..., None] * e1 + (-V)[..., None] * n
    ijk = (Ai[:3, :3] @ pts.reshape(-1, 3).T).T + Ai[:3, 3]
    return ndi.map_coordinates(vol, ijk.T, order=order, cval=-1000 if order else 0).reshape(U.shape)


def trace(s, Ls):
    """Rows: depth, wall x, LA edge x, gap (nan if rejected), accepted flag."""
    sm = ndi.gaussian_filter(s, SMOOTH_SIGMA)
    out, prev, miss = [], None, 0
    for j, v in enumerate(GV):
        if v < START:
            continue
        la = np.where(Ls[:, j] == 2)[0]
        if len(la) == 0:
            out.append((v, np.nan, np.nan, np.nan, 0)); miss += 1
            if prev is not None and miss > MAXMISS: break
            continue
        xla = la.min()
        lo, hi = max(xla - WIN_FRONT, 0), xla + WIN_BACK
        if prev is not None:
            lo, hi = max(lo, prev - STEP), min(hi, prev + STEP + 1)
        if hi <= lo:
            out.append((v, np.nan, GU[xla], np.nan, 0)); miss += 1
            if miss > MAXMISS: break
            continue
        x = lo + int(np.argmin(sm[lo:hi, j]))
        prom = peak_prominences(-sm[:, j], [x], wlen=PROM_WLEN)[0][0]
        if prom >= PROM and lo < x < hi - 1:
            prev, miss = x, 0
            out.append((v, GU[x], GU[xla], (xla - x) * H, 1))
        else:
            miss += 1
            out.append((v, GU[x], GU[xla], np.nan, 0))
            if miss > MAXMISS: break
    return np.array(out) if out else np.zeros((0, 5))


def curtain_end_plane(r):
    if len(r) == 0: return None
    v, gap = r[:, 0], r[:, 3]
    good = np.isfinite(gap)
    if good.sum() < 30: return None
    g = ndi.uniform_filter1d(np.interp(v, v[good], gap[good]), GAP_SMOOTH)
    last = v[good].max()
    for i in np.where(v >= START)[0]:
        j = i + LOOK
        if j >= len(v) or v[j] > last: break
        if g[i] > GAP_THR and np.all(np.diff(g[i:j:5]) >= -0.05) and g[j] > g[i] + RISE:
            return float(v[i])
    return None


def combine(ends):
    p = [e for e in ends if e is not None]
    if len(p) >= 3: return float(np.median(p)), f'{len(p)} planes'
    if len(p) == 2 and abs(p[0] - p[1]) <= AGREE2: return float(np.mean(p)), '2 planes agree'
    return None, f'not found ({len(p)} plane(s)' + (f': {p[0]:.1f} vs {p[1]:.1f}' if len(p) == 2 else '') + ')'


def shortaxis(LV, A, c, n, e1, e2, d):
    Ai = np.linalg.inv(A)
    g = np.arange(-SA_HALF, SA_HALF + SA_H / 2, SA_H)
    U, V = np.meshgrid(g, g, indexing='ij')
    o = c - d * n
    pts = o[None, None, :] + U[..., None] * e1 + V[..., None] * e2
    ijk = (Ai[:3, :3] @ pts.reshape(-1, 3).T).T + Ai[:3, 3]
    m = ndi.map_coordinates(LV.astype(np.uint8), ijk.T, order=0).reshape(U.shape) > 0
    lab, nl = ndi.label(m)
    if nl == 0: return None
    cen = np.array(m.shape) // 2
    if lab[tuple(cen)] > 0: k = lab[tuple(cen)]
    else:
        cc = np.argwhere(lab > 0); k = lab[tuple(cc[np.argmin(((cc - cen) ** 2).sum(1))])]
    R = lab == k
    q = np.argwhere(R) * SA_H; q = q - q.mean(0)
    ev = np.sort(np.linalg.eigvalsh(np.cov(q.T)))[::-1]
    area = R.sum() * SA_H ** 2
    return dict(area_mm2=float(area), ellipticity=float(np.sqrt(ev[0] / ev[1])), Deq_mm=float(2 * np.sqrt(area / np.pi)),
                major_mm=float(4 * np.sqrt(ev[0])), minor_mm=float(4 * np.sqrt(ev[1])))


def qc_png(path, cid, planes, res):
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    f, axs = plt.subplots(1, 5, figsize=(20, 4.8))
    for ax, (off, (s, Ls, r, e)) in zip(axs, planes.items()):
        ax.imshow(s.T, cmap='gray', vmin=res['_myo'] - 150, vmax=res['_blood'] + 150, origin='upper', extent=[GU[0], GU[-1], GV[-1], GV[0]])
        ax.contour(GU, GV, (Ls == 2).T.astype(float), [0.5], colors='magenta', linewidths=0.7)
        if len(r):
            ok = r[:, 4] == 1
            ax.plot(r[ok, 1], r[ok, 0], 'c.', ms=1.2); ax.plot(r[~ok, 1], r[~ok, 0], 'rx', ms=3)
        ax.axhline(0, color='white', ls='--', lw=0.8)
        if e is not None: ax.axhline(e, color='yellow', ls='--', lw=1.3)
        if res['curtain_end_mm'] is not None: ax.axhline(res['curtain_end_mm'], color='lime', ls=':', lw=1.3)
        ax.set_title(f'plane {off:+d}: ' + (f'{e:.1f} mm' if e is not None else 'no end'), fontsize=10)
        ax.set_xlim(-2, 18); ax.set_ylim(15, -5); ax.set_xticks([])
    ce = res['curtain_end_mm']
    plt.suptitle(f"case {cid}: curtain end {ce:.1f} mm ({res['curtain_rule']})" if ce is not None else f"case {cid}: {res['curtain_rule']}"
                 + ('  [SYSTOLIC FLAG]' if res['systolic_flag'] else ''), fontsize=12)
    plt.tight_layout(rect=(0, 0, 1, 0.92)); plt.savefig(path, dpi=60); plt.close(f)


def run_case(cid, img_dir, lab_dir, out_dir):
    jpath = os.path.join(out_dir, 'cases', f'{cid}.json')
    if os.path.exists(jpath): return 'skip'
    res = dict(id=cid)
    try:
        ln = nib.load(os.path.join(lab_dir, f'{cid}.img.nii.gz'))
        L = np.asarray(ln.dataobj).astype(np.uint8); A = ln.affine
        im = nib.load(os.path.join(img_dir, f'{cid}.img.nii.gz'))
        if not np.allclose(im.affine, A, atol=1e-3): raise RuntimeError('image and label affines differ')
        I = np.asarray(im.dataobj).astype(np.float32)
        vx = abs(np.linalg.det(A[:3, :3])) / 1000
        res.update(LV_ml=float((L == 3).sum() * vx), myo_ml=float((L == 1).sum() * vx), LA_ml=float((L == 2).sum() * vx))
        res['LV_myo_ratio'] = res['LV_ml'] / res['myo_ml'] if res['myo_ml'] > 0 else None
        res['systolic_flag'] = bool(res['LV_myo_ratio'] is not None and res['LV_myo_ratio'] < SYSTOLE_RATIO)
        res['_blood'], res['_myo'] = float(np.median(I[L == 3])), float(np.median(I[L == 1]))
        c, n, e1, e2 = frame(L, A)
        planes = {}
        for off in OFFS:
            Ls = longaxis(L, A, c, n, e1, e2, off, 0); s = longaxis(I, A, c, n, e1, e2, off, 1)
            r = trace(s, Ls); planes[off] = (s, Ls, r, curtain_end_plane(r))
        ends = [planes[o][3] for o in OFFS]
        res['plane_ends_mm'] = dict(zip([str(o) for o in OFFS], ends))
        ce, rule = combine(ends)
        if res['systolic_flag'] and ce is None: rule += '; systolic: mark by hand'
        res['curtain_end_mm'], res['curtain_rule'] = ce, rule
        LV = L == 3
        levels = {'2mm': 2.0}
        if ce is not None:
            levels.update(mid=ce / 2, curtain_end=ce)
        for k, d in levels.items():
            m = shortaxis(LV, A, c, n, e1, e2, d)
            if m:
                for q, val in m.items(): res[f'{q}_{k}'] = val
            res[f'depth_{k}'] = d
        os.makedirs(os.path.join(out_dir, 'qc'), exist_ok=True)
        qc_png(os.path.join(out_dir, 'qc', f'{cid}.png'), cid, planes, res)
    except Exception as e:
        res['error'] = repr(e); res['trace'] = traceback.format_exc(limit=2)
    res = {k: v for k, v in res.items() if not k.startswith('_')}
    os.makedirs(os.path.join(out_dir, 'cases'), exist_ok=True)
    with open(jpath + '.tmp', 'w') as fh: json.dump(res, fh)
    os.replace(jpath + '.tmp', jpath)
    return res.get('curtain_end_mm', res.get('error'))


def _wrap(args):
    return run_case(*args)


def collect(out_dir):
    import csv
    rows = [json.load(open(p)) for p in glob.glob(os.path.join(out_dir, 'cases', '*.json'))]
    rows.sort(key=lambda r: int(r['id']))
    for r in rows:
        for o, v in (r.pop('plane_ends_mm', None) or {}).items(): r[f'plane_{o}_end_mm'] = v
        r.pop('trace', None)
    keys = []
    for r in rows:
        for k in r:
            if k not in keys: keys.append(k)
    with open(os.path.join(out_dir, 'lvot_curtain_results.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=keys); w.writeheader(); w.writerows(rows)
    found = sum(r.get('curtain_end_mm') is not None for r in rows)
    print(f'{len(rows)} cases; curtain end found in {found}; errors {sum("error" in r for r in rows)}')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--img-dir', required=True, nargs='+', help='one or more image folders'); ap.add_argument('--lab-dir', required=True); ap.add_argument('--out', required=True)
    ap.add_argument('--ids', default=None, help='comma-separated ids; default = all labels that have an image')
    ap.add_argument('--workers', type=int, default=2)
    a = ap.parse_args()
    where = {}
    for d in a.img_dir:
        for p in glob.glob(os.path.join(d, '*.img.nii.gz')): where[int(os.path.basename(p).split('.')[0])] = d
    if a.ids: ids = [int(x) for x in a.ids.split(',')]
    else:
        ids = sorted(i for i in where if os.path.exists(os.path.join(a.lab_dir, f'{i}.img.nii.gz')))
    ids = [i for i in ids if i in where]
    ids = [i for i in ids if i not in (90, 141)]   # invalid DTU labels
    os.makedirs(a.out, exist_ok=True)
    if a.workers > 1:
        from multiprocessing import Pool
        with Pool(a.workers) as p:
            for cid, r in zip(ids, p.imap(_wrap, [(i, where[i], a.lab_dir, a.out) for i in ids])):
                print(cid, r, flush=True)
    else:
        for i in ids: print(i, run_case(i, where[i], a.lab_dir, a.out), flush=True)
    collect(a.out)
