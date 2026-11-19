"""
CT LVOT summary without the scans whose 2 mm plane still cuts an aortic sinus (version 2, 2026-11-18).

Inputs (unchanged): lvot_out/lvot_curtain_results.csv and sinus_2mm_out/sinus_2mm_flags.csv.
Outputs (new; the version-1 figure and table are kept):
  lvot_out/lvot_curtain_summary_v2.png    the same six panels as plot_summary.py, on the kept scans
  lvot_out/summary_v1_vs_v2.csv           key numbers, all scans (v1) vs kept scans (v2)

Level names: '2 mm' is 2 mm below the LV-aorta label boundary (not 2 mm below the annulus; see the limitations).

Run from the imagecas folder:
  uv run --python 3.12 --with pandas --with numpy --with matplotlib python -I lvot_out/summary_sinus_excluded.py
"""
import os
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
res = pd.read_csv(os.path.join(HERE, 'lvot_curtain_results.csv'))
fl = pd.read_csv(os.path.join(ROOT, 'sinus_2mm_out', 'sinus_2mm_flags.csv'))
assert set(fl.id) == set(res.id), 'sinus check incomplete'
assert 'error' not in fl.columns or fl['error'].isna().all(), 'sinus check has errors'
d = res.merge(fl[['id', 'exclude', 'sinus_inside', 'sinus_beside']], on='id', how='left')
keep = d[~d.exclude.astype(bool)].copy()
LEV = [('2mm', '2 mm below label boundary'), ('mid', 'mid-curtain'), ('curtain_end', 'curtain end')]


def numbers(x, tag):
    f = x[x.curtain_end_mm.notna()]
    row = dict(version=tag, n_scans=len(x), n_curtain_found=len(f), pct_found=100 * len(f) / len(x),
               n_systolic_flag=int(x.systolic_flag.sum()), curtain_end_median_mm=f.curtain_end_mm.median(),
               curtain_end_q1_mm=f.curtain_end_mm.quantile(.25), curtain_end_q3_mm=f.curtain_end_mm.quantile(.75))
    for k, _ in LEV:
        s = x if k == '2mm' else f        # 2 mm exists in every scan; mid and curtain end only where the end was found
        row[f'ellipticity_{k}_median'] = s[f'ellipticity_{k}'].median()
        row[f'area_{k}_median_mm2'] = s[f'area_mm2_{k}'].median()
        row[f'major_{k}_median_mm'] = s[f'major_mm_{k}'].median(); row[f'minor_{k}_median_mm'] = s[f'minor_mm_{k}'].median()
        # echo-style circle on the minor (AP) axis vs the measured area, %
        row[f'echo_circle_err_{k}_median_pct'] = (100 * (np.pi / 4 * s[f'minor_mm_{k}'] ** 2 / s[f'area_mm2_{k}'] - 1)).median()
        row[f'ellipticity_{k}_found_only_median'] = f[f'ellipticity_{k}'].median()
    return row


tab = pd.DataFrame([numbers(d, 'v1 all scans'), numbers(keep, 'v2 sinus at 2 mm excluded')])
tab.loc[len(tab)] = dict(version='excluded scans', n_scans=int(d.exclude.sum()),
                         n_systolic_flag=int(d[d.exclude.astype(bool)].systolic_flag.sum()))
tab.to_csv(os.path.join(HERE, 'summary_v1_vs_v2.csv'), index=False)
print(tab.T.to_string())

# same figure as plot_summary.py, on the kept scans: run it with the data swapped in and the output renamed
src = open(os.path.join(HERE, 'plot_summary.py')).read()
src = src.replace('d=pd.read_csv("lvot_out/lvot_curtain_results.csv")', 'd=KEEP')
src = src.replace("'lvot_out/lvot_curtain_summary.png'", "'lvot_out/lvot_curtain_summary_v2.png'")
src = src.replace("f'{len(d)} cases (90, 141 excluded).", f"f'{{len(d)}} cases (90, 141 and {int(d.exclude.sum())} scans with a sinus at 2 mm excluded). 2 mm = 2 mm below the LV-aorta label boundary.")
src = src.replace("xl=['2 mm','mid','curtain\\nend']", "xl=['2 mm*','mid','curtain\\nend']")
src = src.replace("xlabel='mm below aortic valve plane'", "xlabel='mm below the LV-aorta label boundary'")
src = src.replace("2 mm below the LV-aorta label boundary. Panels D", "2 mm below the LV-aorta label boundary.\\nPanels D")
src = src.replace("fig.text(0.01,0.935,", "fig.text(0.01,0.912,")
src = src.replace("plt.tight_layout(rect=(0,0,1,0.93))", "plt.tight_layout(rect=(0,0,1,0.895))")
assert src.count('KEEP') == 1 and 'summary_v2.png' in src and 'sinus at 2 mm excluded' in src and "'2 mm*'" in src and 'LV-aorta label boundary\'' in src and 'rect=(0,0,1,0.895)' in src, 'plot_summary.py changed; update the replacements'
cwd = os.getcwd(); os.chdir(ROOT)
try:
    exec(compile(src, 'plot_summary_v2', 'exec'), dict(KEEP=keep.drop(columns=['exclude', 'sinus_inside', 'sinus_beside'])))
finally:
    os.chdir(cwd)
print('wrote lvot_out/lvot_curtain_summary_v2.png and lvot_out/summary_v1_vs_v2.csv')
