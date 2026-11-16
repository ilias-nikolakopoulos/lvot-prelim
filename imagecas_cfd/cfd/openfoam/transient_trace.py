"""Time traces of the PW sample-volume velocity in a transient run (probes every time step).
Run from the imagecas folder:
  uv run --python 3.12 --with numpy --with pandas --with matplotlib python -I cfd/openfoam/transient_trace.py --tag lv_apical_lam_tr_h1
"""
import argparse, glob, json, os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from postprocess import read_probe, CFD  # noqa: E402
ap = argparse.ArgumentParser(); ap.add_argument('--tag', default='lv_apical_lam_tr_h1'); args = ap.parse_args()
cases = sorted(glob.glob(os.path.join(CFD, 'runs', args.tag, '*_*')))
fig, axs = plt.subplots(len(cases), 1, figsize=(10, 3.2 * len(cases)), sharex=True, squeeze=False)
cols = {5.0: '#1f5fbf', 7.5: '#d07a00', 10.0: '#2a9d4b'}
for ax, case in zip(axs[:, 0], cases):
    cj = json.load(open(os.path.join(case, 'case.json')))
    info = json.load(open(os.path.join(CFD, 'geometry', str(cj['id']), 'slices.json')))
    sl = {round(s['depth_mm'], 2): s for s in info['slices']}
    it, U = read_probe(case, 'probesU', 'U'); U = U.reshape(len(it), -1, 3) @ np.array(cj['n'])
    for d in cj['echo_depths_mm']:
        idx = [i for i, m in enumerate(cj['probe_meta']) if m['depth_mm'] == d]
        s = sl[round(d, 2)]; tw = cj['twin']
        A = s['real']['area_mm2'] if tw == 'real' else s['circle']['area_mm2']
        ap_mm = s['real']['AP_mm'] if tw == 'real' else 2 * np.sqrt(A / np.pi)
        ar = np.pi * ap_mm ** 2 / 4 / A
        prof = U[:, idx].mean(1) / (cj['Q_m3s'] / (A * 1e-6))
        w = it >= cj['t_avg_s']
        ax.plot(it, prof, color=cols[d], lw=1, label=f'{d:g} mm: mean {prof[w].mean():.3f}, SD {100*prof[w].std():.1f}%  (echo SV/true {ar*prof[w].mean():.2f})')
    ax.axvspan(cj['t_avg_s'], it[-1], color='0.92', zorder=0)
    ax.axhline(1, color='k', lw=0.6)
    ax.set_ylabel('gate velocity / mean'); ax.set_title(f"case {cj['id']} {cj['twin']} ({args.tag}): PW-gate velocity over time at constant peak flow", fontsize=10)
    ax.legend(fontsize=8, loc='lower right'); ax.grid(alpha=0.3)
axs[-1, 0].set_xlabel('time (s); grey = averaging window')
out = os.path.join(CFD, 'qc', f'transient_trace_{args.tag}.png'); fig.tight_layout(); fig.savefig(out, dpi=140); print('wrote', os.path.relpath(out, CFD))
