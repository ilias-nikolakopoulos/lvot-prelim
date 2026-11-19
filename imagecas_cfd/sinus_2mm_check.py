#!/usr/bin/env python3
"""
Flag scans whose 2 mm LVOT plane still cuts an aortic sinus (rule agreed 2026-11-18).

Why: 0 mm in lvot_curtain.py is the LV-aorta label boundary. It sits about 1-2 mm above the annulus, so in some scans
the 2 mm plane is still at sinus level. These scans are excluded from the CT analysis. The annulus correction
(annulus_offset.py checked against MPR hand marks) is tabled as a possible later fix.

Rule, on the 2 mm short-axis plane of lvot_curtain.py (labels only, no CT needed):
  'inside': aorta label inside the central LV+aorta section, area >= INSIDE_MIN_FRAC of the section;
  'beside': an aorta-label region of at least BESIDE_MIN_MM2 that comes within BESIDE_REACH_MM of the section.
A scan is excluded if 'inside' or 'beside' is true.

This script does not change lvot_curtain.py or its results.

Usage (from this folder, without -I because it imports lvot_curtain):
  uv run --python 3.12 --with nibabel --with scipy --with numpy python sinus_2mm_check.py --out sinus_2mm_out [--ids 1,2,3]
Outputs: sinus_2mm_out/cases/{id}.json (resumable), sinus_2mm_out/sinus_2mm_flags.csv
"""
import argparse, csv, glob, json, os
import numpy as np
import nibabel as nib
from scipy import ndimage as ndi
import lvot_curtain as lc

LAB_DIR = 'labels_DTU/ImageCAS-STACOM2025-02-10-2025/segmentations'
# ---------------- parameters ----------------
DEPTH = 2.0              # mm below the current 0 mm (the primary LVOT level)
H = 0.25                 # mm, short-axis grid
HALF = 30.0              # mm, half-width of the grid
INSIDE_MIN_FRAC = 0.01   # 'inside': aorta label >= 1% of the section (the smallest clear sinus tips in the review: 1-2%)
BESIDE_MIN_MM2 = 5.0     # 'beside': ignore aorta-label specks smaller than this
BESIDE_REACH_MM = 3.0    # 'beside': the region must come this close to the section


def check(L, A):
    Ai = np.linalg.inv(A)
    c, n, e1, e2 = lc.frame(L, A)
    g = np.arange(-HALF, HALF + H / 2, H); U, V = np.meshgrid(g, g, indexing='ij')
    pts = (c - DEPTH * n) + U[..., None] * e1 + V[..., None] * e2
    ijk = (Ai[:3, :3] @ pts.reshape(-1, 3).T).T + Ai[:3, 3]
    lab = ndi.map_coordinates(L, ijk.T, order=0).reshape(U.shape)
    m = (lab == 3) | (lab == 6)
    cc, nl = ndi.label(m)
    if nl == 0: raise RuntimeError('no LV or aorta label at 2 mm')
    cen = np.array(m.shape) // 2; k = cc[tuple(cen)]
    if k == 0:
        p = np.argwhere(cc > 0); k = cc[tuple(p[np.argmin(((p - cen) ** 2).sum(1))])]
    sec = ndi.binary_fill_holes(cc == k)
    ao = lab == 6
    inside = float((ao & sec).sum() * H * H); sec_a = float(sec.sum() * H * H)
    near = ndi.binary_dilation(sec, iterations=int(round(BESIDE_REACH_MM / H))) & ~sec
    cc2, n2 = ndi.label(ao & ~sec); beside = 0.0
    for j in range(1, n2 + 1):
        comp = cc2 == j; a = float(comp.sum() * H * H)
        if a >= BESIDE_MIN_MM2 and (comp & near).any(): beside += a
    f_in = inside / sec_a >= INSIDE_MIN_FRAC; f_be = beside > 0
    return dict(section_mm2=sec_a, ao_inside_mm2=inside, ao_inside_frac=inside / sec_a, ao_beside_mm2=beside,
                sinus_inside=bool(f_in), sinus_beside=bool(f_be), exclude=bool(f_in or f_be))


def run_case(cid, out):
    jp = os.path.join(out, 'cases', f'{cid}.json')
    if os.path.exists(jp): return 'skip'
    res = dict(id=cid)
    try:
        ln = nib.load(os.path.join(LAB_DIR, f'{cid}.img.nii.gz'))
        res.update(check(np.asarray(ln.dataobj).astype(np.uint8), ln.affine))
    except Exception as e:
        res['error'] = repr(e)
    os.makedirs(os.path.join(out, 'cases'), exist_ok=True)
    with open(jp + '.tmp', 'w') as fh: json.dump(res, fh)
    os.replace(jp + '.tmp', jp)
    return res.get('exclude', res.get('error'))


def collect(out):
    rows = sorted((json.load(open(p)) for p in glob.glob(os.path.join(out, 'cases', '*.json'))), key=lambda r: int(r['id']))
    keys = []
    for r in rows:
        for k in r:
            if k not in keys: keys.append(k)
    with open(os.path.join(out, 'sinus_2mm_flags.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=keys); w.writeheader(); w.writerows(rows)
    ex = [r for r in rows if r.get('exclude')]
    print(f"{len(rows)} cases; exclude {len(ex)} (inside {sum(r.get('sinus_inside', False) for r in rows)}, "
          f"beside {sum(r.get('sinus_beside', False) for r in rows)}); errors {sum('error' in r for r in rows)}")


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', required=True); ap.add_argument('--ids', default=None)
    a = ap.parse_args()
    have = {int(os.path.basename(p).split('.')[0]) for p in glob.glob('ImageCAS_*/*/*.img.nii.gz')}
    ids = [int(x) for x in a.ids.split(',')] if a.ids else sorted(i for i in have if os.path.exists(os.path.join(LAB_DIR, f'{i}.img.nii.gz')))
    ids = [i for i in ids if i not in (90, 141)]   # invalid DTU labels
    os.makedirs(a.out, exist_ok=True)
    for i in ids: print(i, run_case(i, a.out), flush=True)   # one label volume in memory at a time
    collect(a.out)
