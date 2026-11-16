"""Pulsatile run (lv_lam_pulse_h1): emulated PW-Doppler VTI over one systole vs the true VTI.

The run (make_cases.py --pulse) ejects one systole from rest: raised-cosine flow rise to Q = 400 mL/s at T_PEAK = 0.11 s,
fall to 0 at T_SYS = 0.33 s (SV 66 mL). echoGateFO samples U in the emulated PW gate (echo_spectrum.py geometry: 4 mm
range gate, 3 mm FWHM beam, at the slice centroid, beam along the LVOT axis) every 1 ms at 5 / 7.5 / 10 mm.

For each case, depth and time step:
  spectrum of U.n in the gate -> mean, modal, max (97.5th pct)   (echo_spectrum.spectrum)
  v_true(t) = Q(t) / A_true                                       (true spatial mean velocity)
  ring(t)   = mean U.n over the centroid + 1.5 mm ring probes     (postProcessing/probesU, every time step; cross-check)
VTI_x = integral of x(t) dt from 0 to T_SYS (trapezoid, v = 0 at t = 0)
VTI_ratio = VTI_readout / VTI_true;  SV_ratio = area_ratio x VTI_ratio  (area_ratio from cfd_results_lv_lam_h1.csv)
The steady values (echo_spectrum_lv_lam_h1.csv, angle 0) are the peak-flow-only estimate it is compared with.

Run from the imagecas folder:
  uv run --python 3.12 --with numpy --with scipy --with pandas --with matplotlib python -I cfd/openfoam/pulse_vti.py
Outputs: cfd/pulse_vti_lv_lam_pulse_h1.csv, cfd/pulse_vti_timeseries_lv_lam_pulse_h1.csv, cfd/qc/pulse_vti_lv_lam_pulse_h1.png
"""
import glob, json, os, sys
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from echo_spectrum import gate_points, world, frame, spectrum, DEPTHS, MAX_PCT   # noqa: E402
from make_cases import pulse_q, T_SYS, T_PEAK                                    # noqa: E402
from postprocess import read_probe                                               # noqa: E402

CFD = os.path.dirname(HERE)
TAG = 'lv_lam_pulse_h1'
STEADY = 'lv_lam_h1'
COLS = {'mean': '#2a78d6', 'modal': '#d1495b', 'max': '#7a3db8', 'ring': '#2a9d4b'}


def ring_series(case, cj, d):
    t, U = read_probe(case, 'probesU', 'U')
    U = U.reshape(len(t), -1, 3) @ np.array(cj['n'])
    idx = [i for i, m in enumerate(cj['probe_meta']) if m['depth_mm'] == d and abs(U[-1, i]) < 1e9]
    return np.asarray(t, float), U[:, idx].mean(1)


def run_case(case, res):
    cj = json.load(open(os.path.join(case, 'case.json')))
    fr, sl = frame(cj['id'])
    rows = []
    for d in DEPTHS:
        A, U, V, w = gate_points(fr, sl[d]['real'])
        tree = cKDTree(world(fr, sl[d]['real'], d, A, U, V))
        r0 = res.loc[(cj['id'], cj['twin'], d)]
        a_true = r0['A_true_mm2'] * 1e-6
        name = 'g' + str(d).replace('.', 'p') + '_U.xy'
        for f in glob.glob(os.path.join(case, 'postProcessing', 'echoGate', '*', name)):
            t = float(os.path.basename(os.path.dirname(f)))
            D = np.loadtxt(f)
            dist, idx = tree.query(D[:, :3]); ok = dist < 1e-6
            _, _, modal, vmax, mean = spectrum(D[ok, 3:6] @ fr['n'], w[idx[ok]])
            rows.append(dict(id=cj['id'], twin=cj['twin'], depth_mm=d, t=t, v_true=float(pulse_q(t)) / a_true,
                             mean=mean, modal=modal, max=vmax, area_ratio=r0['area_ratio'], frac_in_mesh=ok.mean()))
    ts = pd.DataFrame(rows).sort_values(['depth_mm', 't'])
    # ring probes (every solver step), interpolated onto the gate times
    for d in DEPTHS:
        tr, vr = ring_series(case, cj, d)
        k = ts.depth_mm == d
        ts.loc[k, 'ring'] = np.interp(ts.loc[k, 't'], tr, vr)
    return cj, ts


def vti(t, v):
    t = np.r_[0.0, t]; v = np.r_[0.0, v]
    return np.trapezoid(v, t)


def main():
    res = pd.read_csv(os.path.join(CFD, f'cfd_results_{STEADY}.csv')).set_index(['id', 'twin', 'depth_mm'])
    sel = pd.read_csv(os.path.join(CFD, 'case_selection.csv'))[['id', 'pair', 'arm']]
    steady = pd.read_csv(os.path.join(CFD, f'echo_spectrum_{STEADY}.csv'))
    steady = steady[steady.angle_deg == 0].pivot_table(index=['id', 'twin', 'depth_mm'], columns='readout', values='profile_ratio')
    allts, out = [], []
    for case in sorted(glob.glob(os.path.join(CFD, 'runs', TAG, '*_*'))):
        if not os.path.isfile(os.path.join(case, 'solve.done')):
            continue
        cj, ts = run_case(case, res)
        allts.append(ts)
        for d, g in ts.groupby('depth_mm'):
            g = g[g.t <= T_SYS + 1e-9]
            vt = vti(g.t.values, g.v_true.values)
            row = dict(id=cj['id'], twin=cj['twin'], depth_mm=d, n_times=len(g), t_last=g.t.max(),
                       VTI_true_cm=100 * vt, area_ratio=g.area_ratio.iloc[0], frac_in_mesh_min=g.frac_in_mesh.min())
            pk = g.iloc[(g.t - T_PEAK).abs().argmin()]
            for ro in ('mean', 'modal', 'max', 'ring'):
                r = vti(g.t.values, g[ro].values) / vt
                row[f'VTI_ratio_{ro}'] = r
                row[f'SV_ratio_{ro}'] = row['area_ratio'] * r
                row[f'peak_ratio_{ro}'] = pk[ro] / pk.v_true
                if ro != 'ring':
                    row[f'steady_ratio_{ro}'] = steady.loc[(cj['id'], cj['twin'], d), ro]
            out.append(row)
        print('done', os.path.basename(case))
    df = pd.DataFrame(out).merge(sel, on='id')
    ts = pd.concat(allts)
    df.to_csv(os.path.join(CFD, f'pulse_vti_{TAG}.csv'), index=False)
    ts.to_csv(os.path.join(CFD, f'pulse_vti_timeseries_{TAG}.csv'), index=False)
    pd.set_option('display.width', 250)
    show = ['id', 'arm', 'twin', 'depth_mm', 'n_times', 'area_ratio'] + [f'{p}_{ro}' for ro in ('max', 'modal', 'mean') for p in ('steady_ratio', 'peak_ratio', 'VTI_ratio')] + ['VTI_ratio_ring', 'SV_ratio_max', 'SV_ratio_mean']
    print(df[show].round(3).to_string(index=False))
    figure(df, ts)


def figure(df, ts):
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    fig, axs = plt.subplots(2, 3, figsize=(15, 9), gridspec_kw=dict(hspace=0.45, wspace=0.28))
    d = 7.5
    for k, cid in enumerate((142, 388)):
        ax = axs[0, k]; arm = df[df.id == cid].arm.iloc[0]
        g = ts[(ts.id == cid) & (ts.twin == 'real') & (ts.depth_mm == d)]
        ax.plot(g.t * 1e3, g.v_true, 'k--', lw=1.5, label='true mean (Q/A)')
        for ro in ('mean', 'modal', 'max'):
            ax.plot(g.t * 1e3, g[ro], color=COLS[ro], lw=1.4, label=f'gate {ro}' + (f' ({MAX_PCT:g}th pct)' if ro == 'max' else ''))
        ax.set_xlabel('time from start of ejection (ms)'); ax.set_ylabel('velocity along LVOT axis (m/s)')
        ax.set_title(f'{"AB"[k]}  Case {cid} ({arm}), real LVOT, {d:g} mm', fontsize=9.5, loc='left')
        ax.legend(fontsize=7.5, frameon=False)
    ax = axs[0, 2]
    for cid, ls in ((142, '-'), (388, '--')):
        for tw, c in (('real', '#7a3db8'), ('circle', '#9a9a9a')):
            g = ts[(ts.id == cid) & (ts.twin == tw) & (ts.depth_mm == d)]
            g = g[g.v_true > 0.15 * g.v_true.max()]
            ax.plot(g.t * 1e3, g['max'] / g.v_true, ls=ls, color=c, lw=1.4, label=f'{cid} {tw}')
    ax.axvline(T_PEAK * 1e3, color='k', lw=0.6, ls=':'); ax.axhline(1, color='k', lw=0.6)
    ax.set_xlabel('time (ms; shown where Q > 15% of peak)'); ax.set_ylabel('gate max ÷ true mean')
    ax.set_title(f'C  Over-read through systole (max read-out, {d:g} mm)', fontsize=9.5, loc='left'); ax.legend(fontsize=7.5, frameon=False)
    # bottom: steady peak ratio vs pulsatile VTI ratio, per read-out, real twin, dots = depths
    for k, (cid) in enumerate((142, 388)):
        ax = axs[1, k]; q = df[(df.id == cid) & (df.twin == 'real')]
        xt, xl = [], []
        for j, ro in enumerate(('mean', 'modal', 'max')):
            for i, (col, lab, a) in enumerate(((f'steady_ratio_{ro}', 'steady', 0.35), (f'VTI_ratio_{ro}', 'pulse', 0.8))):
                x = j * 3 + i
                ax.bar(x, q[col].median(), color=COLS[ro], alpha=a, width=0.8)
                ax.plot([x] * len(q), q[col], 'o', ms=3.5, color='k')
                xt.append(x); xl.append(f'{ro}\n{lab}')
        ax.axhline(1, color='k', lw=0.8); ax.set_ylim(0.9, 1.45)
        ax.set_xticks(xt, xl, fontsize=7.5); ax.set_ylabel('read-out ÷ true (velocity or VTI)')
        ax.set_title(f'{"DE"[k]}  Case {cid} ({q.arm.iloc[0]}): steady peak vs pulsatile VTI\n(real LVOT; dots = 5 / 7.5 / 10 mm)', fontsize=9.5, loc='left')
    ax = axs[1, 2]
    z = df[df.twin == 'real']
    lab, area, vtis, sv = [], [], [], []
    for cid in (142, 388):
        for ro in ('modal', 'max'):
            q = z[z.id == cid]
            lab.append(f'{cid} ({q.arm.iloc[0][:3]}), {ro}'); area.append(100 * (q.area_ratio.median() - 1))
            vtis.append(100 * (q[f'VTI_ratio_{ro}'].median() - 1)); sv.append(100 * (q[f'SV_ratio_{ro}'].median() - 1))
    lab.append('Huh 2022\n(echo vs 4D flow)'); area.append(-25.5); vtis.append(29.4); sv.append(3.8)
    y = np.arange(len(lab))[::-1]
    for off, vals, c, nm in ((-0.25, area, '#8a8a8a', 'area'), (0, vtis, '#e08a2c', 'VTI'), (0.25, sv, '#2a9d4b', 'SV')):
        ax.barh(y + off, vals, height=0.24, color=c, label=nm)
    ax.axvline(0, color='k', lw=0.8); ax.set_yticks(y, lab, fontsize=8); ax.set_xlabel('% error of echo vs truth (median over depths)')
    ax.legend(fontsize=8, frameon=False, loc='upper left', bbox_to_anchor=(0, -0.16), ncol=3)
    ax.set_title('F  Pulsatile error split vs Huh et al. 2022', fontsize=9.5, loc='left')
    for a in axs.flat:
        for sp in ('top', 'right'):
            a.spines[sp].set_visible(False)
    fig.suptitle(f'Pulsatile systole (one beat from rest, peak Q 400 mL/s at {T_PEAK * 1e3:.0f} ms, ejection {T_SYS * 1e3:.0f} ms, SV 66 mL), '
                 'patient LV inflow: emulated PW-Doppler VTI', fontsize=11.5)
    fig.text(0.5, 0.01, 'Gate 4 mm long, beam FWHM 3 mm, at the LVOT centroid, beam along the LVOT axis, sampled every 1 ms. VTI = time integral of the '
             'read-out over ejection; VTI ratio = VTI read-out ÷ VTI of Q/A.\nSteady = the same read-out at constant peak flow (lv_lam_h1). '
             'SV ratio = area ratio (AP-diameter circle ÷ true area) × VTI ratio.', ha='center', fontsize=8)
    out = os.path.join(CFD, 'qc', f'pulse_vti_{TAG}.png'); fig.savefig(out, dpi=135, bbox_inches='tight'); plt.close(fig)
    print('wrote', os.path.relpath(out, CFD))


if __name__ == '__main__':
    main()
