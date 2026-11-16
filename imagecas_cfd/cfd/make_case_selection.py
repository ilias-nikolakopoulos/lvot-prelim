#!/usr/bin/env python3
"""Write cfd/case_selection.csv: the 6 systolic/diastolic CFD pairs chosen by hand on 2026-11-12
(see cfd/README_case_selection.md for the criteria and the QC exclusions).
Run from the imagecas folder:
  uv run --python 3.12 --with pandas python -I cfd/make_case_selection.py
"""
import pandas as pd
d = pd.read_csv('lvot_out/lvot_curtain_results.csv').set_index('id')
pairs = [(388, 142), (576, 167), (33, 978), (586, 587), (253, 467), (146, 663)]   # (systolic, diastolic)
note = {146: 'lower contrast; kept for ellipticity range', 576: 'very bright contrast, clean',
        586: 'high noise (210 HU) but clean labels', 142: 'small aorta label islands at +2 mm, not in LVOT lumen'}
rows = []
for k, (s, dd) in enumerate(pairs, 1):
    for arm, i in [('systolic', s), ('diastolic', dd)]:
        r = d.loc[i]
        rows.append(dict(pair=k, arm=arm, id=i, LV_myo_ratio=round(r.LV_myo_ratio, 2), LV_ml=round(r.LV_ml, 1),
                         area_mm2_2mm=round(r.area_mm2_2mm, 1), ellipticity_2mm=round(r.ellipticity_2mm, 3),
                         curtain_end_mm=r.curtain_end_mm if arm == 'diastolic' else None, note=note.get(i, '')))
pd.DataFrame(rows).to_csv('cfd/case_selection.csv', index=False)
print(pd.DataFrame(rows).to_string(index=False))
