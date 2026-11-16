#!/usr/bin/env python3
"""Check twin STLs (watertight, volume) and tabulate the geometric (area) part of the echo error.
Usage: python -I cfd/summarize_twins.py   (run from the imagecas folder)"""
import json, os, numpy as np, pandas as pd, trimesh
HERE = os.path.dirname(os.path.abspath(__file__))
sel = pd.read_csv(os.path.join(HERE, 'case_selection.csv'))
rows = []
for r in sel.itertuples():
    g = os.path.join(HERE, 'geometry', str(r.id)); info = json.load(open(os.path.join(g, 'slices.json')))
    chk = {}
    for tw in ('real', 'circle'):
        m = trimesh.load(os.path.join(g, f'{tw}.stl'), force='mesh'); m.merge_vertices()
        chk[tw] = (m.is_watertight, m.volume)
    s = {round(x['depth_mm'], 1): x['real'] for x in info['slices']}
    row = dict(pair=r.pair, arm=r.arm, id=r.id, watertight=chk['real'][0] and chk['circle'][0],
               vol_diff_pct=100 * (chk['real'][1] / chk['circle'][1] - 1))
    for d in (5.0, 7.5, 10.0):
        row[f'ell_{d:g}'] = s[d]['ellipticity']; row[f'area_{d:g}'] = s[d]['area_mm2']
        row[f'AP_{d:g}'] = s[d]['AP_mm']; row[f'ML_{d:g}'] = s[d]['ML_mm']
        row[f'echo_area_err_pct_{d:g}'] = 100 * (s[d]['echo_area_mm2'] / s[d]['area_mm2'] - 1)
    rows.append(row)
t = pd.DataFrame(rows); t.to_csv(os.path.join(HERE, 'twin_geometry_summary.csv'), index=False)
pd.set_option('display.width', 220)
print(t[['pair', 'arm', 'id', 'watertight', 'vol_diff_pct', 'ell_7.5', 'area_7.5', 'AP_7.5', 'ML_7.5',
         'echo_area_err_pct_5', 'echo_area_err_pct_7.5', 'echo_area_err_pct_10']].round(2).to_string(index=False))
print('\nmedian echo-area error by arm (%):')
print(t.groupby('arm')[['echo_area_err_pct_5', 'echo_area_err_pct_7.5', 'echo_area_err_pct_10', 'ell_7.5']].median().round(2))
