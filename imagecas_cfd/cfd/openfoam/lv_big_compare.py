"""Enlarged-LV sensitivity: does the systolic over-read fall to the diastolic level when the systolic LV is grown to 88 mL?

Compares, for the real twin, the spectral max read-out (echo_spectrum_<tag>.csv, beam angle 0) and SV_ratio of
  lv_lam_h1 systolic      (patient systolic LV, as segmented)
  lv_big_lam_h1           (the same systolic scans, LV grown to 88 mL by build_lv_twins.py --grow-to 88; LVOT untouched)
  lv_lam_h1 diastolic     (the paired diastolic scans)
Run from the imagecas folder:
  uv run --python 3.12 --with numpy --with pandas --with matplotlib python -I cfd/openfoam/lv_big_compare.py
Output: cfd/qc/lv_big_compare.png, printed table.
"""
import json, os
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
CFD = os.path.dirname(HERE)


def main():
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    sel = pd.read_csv(os.path.join(CFD, 'case_selection.csv'))
    s = pd.concat([pd.read_csv(os.path.join(CFD, f'echo_spectrum_{t}.csv')) for t in ('lv_lam_h1', 'lv_big_lam_h1')])
    s = s[(s.angle_deg == 0) & (s.twin == 'real') & (s.readout == 'max')]
    sys_ids = sel[sel.arm == 'systolic'].id.tolist()
    vol = {}
    for cid in sel.id:
        vol[cid] = json.load(open(os.path.join(CFD, 'geometry', str(cid), 'lv_info.json')))['lv_volume_ml']
    rows = []
    for cid in sys_ids:
        dia = int(sel[(sel.pair == sel[sel.id == cid].pair.iloc[0]) & (sel.arm == 'diastolic')].id.iloc[0])
        for lab, tag, i in (('systolic LV', 'lv_lam_h1', cid), ('systolic, LV grown', 'lv_big_lam_h1', cid), ('diastolic LV', 'lv_lam_h1', dia)):
            q = s[(s.tag == tag) & (s.id == i)]
            big = json.load(open(os.path.join(CFD, 'geometry', str(i), 'lv_info_big.json'))) if tag == 'lv_big_lam_h1' else None
            rows.append(dict(pair_sys_id=cid, condition=lab, id=i, lv_ml=big['grow_to_ml'] if big else vol[i],
                             max_ratio=q.profile_ratio.median(), SV_ratio_max=q.SV_ratio.median()))
    df = pd.DataFrame(rows)
    pd.set_option('display.width', 200)
    print(df.round(3).to_string(index=False))
    print(df.groupby('condition', sort=False)[['lv_ml', 'max_ratio', 'SV_ratio_max']].median().round(3))

    fig, axs = plt.subplots(1, 2, figsize=(12, 4.6), gridspec_kw=dict(wspace=0.3))
    conds = ['systolic LV', 'systolic, LV grown', 'diastolic LV']
    cols = ['#d1495b', '#e08a2c', '#2a78d6']
    for ax, col, yl in ((axs[0], 'max_ratio', 'max read-out ÷ true mean velocity'), (axs[1], 'SV_ratio_max', 'echo SV ÷ true SV (max read-out)')):
        labs = []
        for cid, g in df.groupby('pair_sys_id'):
            ax.plot(range(3), g.set_index('condition').loc[conds, col], '-', color='#bbbbbb', lw=1, zorder=1)
            labs.append((g[g.condition == conds[2]][col].iloc[0], f'{cid}/{g[g.condition == conds[2]].id.iloc[0]}'))
        labs.sort(); gap = 0.035 * np.ptp(df[col]); prev = -np.inf
        for y, t in labs:       # push labels apart so they do not overlap
            y = max(y, prev + gap); prev = y
            ax.text(2.08, y, t, fontsize=7, va='center', color='#666')
        for k, c in enumerate(conds):
            y = df[df.condition == c][col]
            ax.scatter([k] * len(y), y, s=28, color=cols[k], zorder=2)
            ax.plot([k - 0.25, k + 0.25], [y.median()] * 2, color='k', lw=2, zorder=3)
        ax.axhline(1, color='k', lw=0.7)
        ax.set_xticks(range(3), [f'{c}\n(median {df[df.condition == c].lv_ml.median():.0f} mL)' for c in conds], fontsize=8.5)
        ax.set_ylabel(yl)
        for sp in ('top', 'right'):
            ax.spines[sp].set_visible(False)
    axs[0].set_title('A  Velocity over-read (spectral max, real LVOT)', fontsize=10, loc='left')
    axs[1].set_title('B  Echo SV error (AP-diameter area × max)', fontsize=10, loc='left')
    fig.suptitle('Enlarged-LV sensitivity: same systolic LVOT, LV cavity grown to 88 mL (steady peak flow, patient LV inflow)', fontsize=11)
    fig.text(0.5, -0.04, 'Lines join one pair (systolic scan / its diastolic scan, labelled sys/dia). Black bar = median. '
             'Each point = median over the 5, 7.5 and 10 mm gates, beam along the LVOT axis.', ha='center', fontsize=8)
    out = os.path.join(CFD, 'qc', 'lv_big_compare.png'); fig.savefig(out, dpi=135, bbox_inches='tight'); plt.close(fig)
    print('wrote', os.path.relpath(out, CFD))


if __name__ == '__main__':
    main()
