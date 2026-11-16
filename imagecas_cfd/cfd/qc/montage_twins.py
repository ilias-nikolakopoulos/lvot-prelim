#!/usr/bin/env python3
"""Stack the per-case twin QC images (cfd/geometry/{id}/qc.png, top half = contours) into
cfd/qc/twins_contours_systolic.png and twins_contours_diastolic.png.
Run from the imagecas folder after build_twins.py:
  uv run --python 3.12 --with matplotlib --with numpy python -I cfd/qc/montage_twins.py
"""
import matplotlib.image as mi, numpy as np
ids = [388, 576, 33, 586, 253, 146, 142, 167, 978, 587, 467, 663]   # 6 systolic, then their 6 diastolic partners
ims = [mi.imread(f'cfd/geometry/{i}/qc.png')[..., :3] for i in ids]; h = ims[0].shape[0] // 2 + 10
mi.imsave('cfd/qc/twins_contours_systolic.png', np.vstack([x[:h] for x in ims[:6]]))
mi.imsave('cfd/qc/twins_contours_diastolic.png', np.vstack([x[:h] for x in ims[6:]]))
print('wrote cfd/qc/twins_contours_{systolic,diastolic}.png')
