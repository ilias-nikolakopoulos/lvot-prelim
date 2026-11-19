#!/usr/bin/env python3
"""
Automatic per-scan annulus level for the ImageCAS LVOT pipeline (option 3 of the valve-plane fix, 2026-11-16).

Problem: lvot_curtain.py puts 0 mm at the LV-aorta label interface. This interface crosses the closed cusps and
the sinus floors, so it is above the true annulus (the level of the cusp nadirs). All depths are shifted by this offset.

Method: in the frame of lvot_curtain.py, step down from the current 0 mm. At each level cut a short-axis section and
test if it is a single lumen:
  1. aorta label inside the central LV+aorta label section: less than AO_MAX of its area;
  2. the CT blood region (half-max threshold, connected to the lumen centre) has no lobes: solidity >= SOL_MIN;
  3. the CT blood region is not much larger than the label section: excess <= EXC_MAX (a sinus floor in or beside the lumen).
The annulus level is the first level where all tests are true and stay true for the next HOLD mm.
Offset = annulus level minus current 0 mm (positive = deeper, toward the LV).

Thresholds (study on 20 scans, 2026-11-16):
  - Test 3 sets the level in most scans. The excess falls smoothly with depth (0.1-0.4 at 0 mm, about 0 at 4 mm), so
    there is no sharp step and the level depends on EXC_MAX: median offset 1.6 / 1.7 / 2.2 mm for 0.15 / 0.10 / 0.06.
    The primary value is 0.10. The offsets for 0.05 and 0.15 and for the label test alone are also saved, so that the
    hand marks (option 5) can show which one is correct. Clean levels (4-8 mm) have excess -0.03 (median), 0.02 (99th pct).
  - Test 2 changes few scans (0.94 to 0.965 moves 3 of 20 by 0.2-0.4 mm). SOL_MIN = 0.94 is the 5th percentile of the
    clean levels (4-8 mm), where the open mitral side also lowers the solidity.

This script does not change lvot_curtain.py or its results. Outputs go to a new folder.

Usage (from this folder, without -I because it imports lvot_curtain):
  uv run --python 3.12 --with nibabel --with scipy --with numpy python annulus_offset.py --out annulus_out [--ids 1,2,3] [--profiles]
Outputs: annulus_out/cases/{id}.json (resumable), annulus_out/annulus_offset.csv
"""
import argparse, csv, glob, json, os
import numpy as np
import nibabel as nib
from scipy import ndimage as ndi
from scipy.spatial import ConvexHull
import lvot_curtain as lc

LAB_DIR = 'labels_DTU/ImageCAS-STACOM2025-02-10-2025/segmentations'
# ---------------- parameters ----------------
H = 0.25                   # mm, short-axis grid
HALF = 30.0                # mm, half-width of the grid
DSTEP = 0.2                # mm between levels
DMAX = 8.0                 # mm, deepest level tested
CT_SIGMA = 0.5             # mm, CT smoothing before the threshold
THR_FRAC = 0.5             # blood threshold = myo + THR_FRAC * (blood - myo)
REACH = 4.0                # mm, the blood region is searched this far outside the label section
AO_MAX = 0.01              # test 1
SOL_MIN = 0.94             # test 2
EXC_MAX = 0.10             # test 3 (primary)
EXC_ALT = (0.05, 0.15)     # also saved, for the comparison with hand marks
HOLD = 1.0                 # mm, the tests must stay true over this distance


def central(m, cen=None):
    cc, nl = ndi.label(m)
    if nl == 0: return np.zeros_like(m)
    cen = np.array(m.shape) // 2 if cen is None else cen
    k = cc[tuple(cen)]
    if k == 0:
        p = np.argwhere(cc > 0); k = cc[tuple(p[np.argmin(((p - cen) ** 2).sum(1))])]
    return cc == k


def solidity(m):
    p = np.argwhere(m)
    if len(p) < 20: return np.nan
    return float(len(p) / ConvexHull(p).volume) if len(p) else np.nan   # pixel count / hull area (2-D 'volume')


def profile(L, I, A):
    """Per-level measures from the current 0 mm down to DMAX."""
    Ai = np.linalg.inv(A)
    blood, myo = float(np.median(I[L == 3])), float(np.median(I[L == 1])); thr = myo + THR_FRAC * (blood - myo)
    c, n, e1, e2 = lc.frame(L, A)
    g = np.arange(-HALF, HALF + H / 2, H); U, V = np.meshgrid(g, g, indexing='ij')
    rows = []
    for d in np.arange(0, DMAX + DSTEP / 2, DSTEP):
        pts = (c - d * n) + U[..., None] * e1 + V[..., None] * e2
        ijk = (Ai[:3, :3] @ pts.reshape(-1, 3).T).T + Ai[:3, 3]
        lab = ndi.map_coordinates(L, ijk.T, order=0).reshape(U.shape)
        ct = ndi.gaussian_filter(ndi.map_coordinates(I, ijk.T, order=1).reshape(U.shape), CT_SIGMA / H)
        la = ndi.binary_fill_holes(central((lab == 3) | (lab == 6)))
        if la.sum() < 50:
            rows.append(dict(d=float(d), ao_frac=np.nan, sol_blood=np.nan, sol_label=np.nan, excess=np.nan, A_label=0.0, A_blood=0.0)); continue
        cen = np.round(np.argwhere(la).mean(0)).astype(int)
        allowed = np.isin(lab, [0, 1, 3, 6]) & ndi.binary_dilation(la, iterations=int(REACH / H))
        bl = ndi.binary_fill_holes(central((ct >= thr) & allowed, cen)) & allowed
        rows.append(dict(d=float(d), ao_frac=float(((lab == 6) & la).sum() / la.sum()), sol_blood=solidity(bl), sol_label=solidity(la),
                         excess=float(bl.sum() / la.sum() - 1), A_label=float(la.sum() * H * H), A_blood=float(bl.sum() * H * H)))
    return rows


def annulus_level(rows, exc_max=EXC_MAX, sol_min=SOL_MIN):
    ok = np.array([(r['ao_frac'] < AO_MAX) and (r['sol_blood'] >= sol_min) and (r['excess'] <= exc_max) for r in rows])
    k = int(round(HOLD / DSTEP))
    for i in range(len(rows) - k):
        if ok[i:i + k + 1].all(): return rows[i]['d']
    return None


def run_case(cid, img, out, keep_profile):
    jp = os.path.join(out, 'cases', f'{cid}.json')
    if os.path.exists(jp): return 'skip'
    res = dict(id=cid)
    try:
        ln = nib.load(os.path.join(LAB_DIR, f'{cid}.img.nii.gz')); L = np.asarray(ln.dataobj).astype(np.uint8); A = ln.affine
        I = np.asarray(nib.load(img).dataobj).astype(np.float32)
        rows = profile(L, I, A); lvl = annulus_level(rows)
        res['annulus_offset_mm'] = lvl
        for e in EXC_ALT: res[f'offset_exc{int(round(100 * e)):02d}_mm'] = annulus_level(rows, exc_max=e)
        res['offset_label_only_mm'] = annulus_level(rows, exc_max=1e9, sol_min=0.0)
        for name, test in (('ao', lambda r: r['ao_frac'] < AO_MAX), ('sol', lambda r: r['sol_blood'] >= SOL_MIN), ('exc', lambda r: r['excess'] <= EXC_MAX)):
            first = next((r['d'] for r in rows if test(r)), None); res[f'first_ok_{name}_mm'] = first   # which test sets the level
        r2 = min(rows, key=lambda r: abs(r['d'] - 2.0)); res.update(sol_blood_at_2mm=r2['sol_blood'], excess_at_2mm=r2['excess'], ao_frac_at_2mm=r2['ao_frac'])
        if keep_profile: res['profile'] = rows
    except Exception as e:
        res['error'] = repr(e)
    os.makedirs(os.path.join(out, 'cases'), exist_ok=True)
    with open(jp + '.tmp', 'w') as fh: json.dump(res, fh)
    os.replace(jp + '.tmp', jp)
    return res.get('annulus_offset_mm', res.get('error'))


def collect(out):
    rows = [json.load(open(p)) for p in glob.glob(os.path.join(out, 'cases', '*.json'))]
    rows.sort(key=lambda r: int(r['id']))
    for r in rows: r.pop('profile', None)
    keys = []
    for r in rows:
        for k in r:
            if k not in keys: keys.append(k)
    with open(os.path.join(out, 'annulus_offset.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=keys); w.writeheader(); w.writerows(rows)
    v = [r['annulus_offset_mm'] for r in rows if r.get('annulus_offset_mm') is not None]
    print(f'{len(rows)} cases; level found in {len(v)}; errors {sum("error" in r for r in rows)}; median offset {np.median(v) if v else float("nan"):.2f} mm')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', required=True); ap.add_argument('--ids', default=None); ap.add_argument('--profiles', action='store_true', help='keep the per-level profile in each json')
    a = ap.parse_args()
    where = {int(os.path.basename(p).split('.')[0]): p for p in glob.glob('ImageCAS_*/*/*.img.nii.gz')}
    ids = [int(x) for x in a.ids.split(',')] if a.ids else sorted(i for i in where if os.path.exists(os.path.join(LAB_DIR, f'{i}.img.nii.gz')))
    ids = [i for i in ids if i in where and i not in (90, 141)]   # invalid DTU labels
    os.makedirs(a.out, exist_ok=True)
    for i in ids: print(i, run_case(i, where[i], a.out, a.profiles), flush=True)   # one CT in memory at a time
    collect(a.out)
