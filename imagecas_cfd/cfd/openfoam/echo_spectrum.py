"""Emulate the PW-Doppler spectral read-out in the LVOT: modal velocity and maximum within the gate, vs the spatial mean.

Gate model (one per echo depth 5 / 7.5 / 10 mm, centred on the slice centroid, beam along the LVOT axis n):
  - range gate: GATE_LEN_MM long along the beam, uniform weighting (rectangular gate)
  - beam profile: Gaussian across the beam, FWHM BEAM_FWHM_MM, truncated at BEAM_R_MM
  - points every GRID_MM inside that cylinder; U sampled with OpenFOAM `sets` (cellPoint interpolation)
Spectrum: histogram of the beam-direction velocity U.b, weighted by the gate/beam weight of each point (blood scatters
uniformly, so Doppler power is proportional to sampled volume). Bins BIN_CM_S wide.
Read-outs (all as profile_ratio = read-out / true mean velocity Q/A):
  mean    weighted mean (what the CFD tables used so far, but over this larger gate)
  modal   peak of the spectrum (3-bin smoothed), = the densest band a sonographer traces
  max     weighted 97.5th percentile (the outer edge of the spectrum; one point cannot set it)
SV_ratio = area_ratio (from cfd_results_<tag>.csv) x profile_ratio. Beam angle: b = cos(a) n + sin(a) e1 for a = 0, +-10,
+-20 deg (tilt in the A5C plane, + = toward the LA / mitral side); no angle correction, as in clinical practice.

Steps (from the imagecas folder):
  uv run --python 3.12 --with numpy python -I cfd/openfoam/echo_spectrum.py write --tags lam_h1,lv_lam_h1
  bash cfd/openfoam/echo_spectrum_sample.sh lam_h1 lv_lam_h1     # postProcess in one container, ~10 s per case
  uv run --python 3.12 --with numpy --with pandas --with matplotlib python -I cfd/openfoam/echo_spectrum.py analyse --tags lam_h1,lv_lam_h1
Outputs: cfd/echo_spectrum_<tag>.csv, cfd/echo_spectrum_summary.csv, cfd/qc/echo_spectrum_*.png
"""
import argparse, glob, json, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CFD = os.path.dirname(HERE)
GATE_LEN_MM = 4.0
BEAM_FWHM_MM = 3.0
BEAM_R_MM = 2.5
GRID_MM = 0.25
BIN_CM_S = 2.0
MAX_PCT = 97.5
ANGLES = (0, 10, -10, 20, -20)
DEPTHS = (5.0, 7.5, 10.0)


def gate_points(fr, s):
    """Local gate grid (along, u, v) in mm and the weights, plus world points (m)."""
    g = np.arange(-BEAM_R_MM, BEAM_R_MM + 1e-9, GRID_MM)
    a = np.arange(-GATE_LEN_MM / 2, GATE_LEN_MM / 2 + 1e-9, GRID_MM)
    A, U, V = np.meshgrid(a, g, g, indexing='ij')
    r = np.hypot(U, V); k = r <= BEAM_R_MM
    A, U, V, r = A[k], U[k], V[k], r[k]
    sig = BEAM_FWHM_MM / 2.3548
    w = np.exp(-0.5 * (r / sig) ** 2)
    return A, U, V, w


def world(fr, s, d, A, U, V):
    c = fr['c'] - d * fr['n'] + s['cu'] * fr['e1'] + s['cv'] * fr['e2']
    return (c + np.outer(-A, fr['n']) + np.outer(U, fr['e1']) + np.outer(V, fr['e2'])) * 1e-3


def frame(cid):
    info = json.load(open(os.path.join(CFD, 'geometry', str(cid), 'slices.json')))
    fr = {k: np.array(v) for k, v in info['frame'].items()}
    sl = {round(x['depth_mm'], 2): x for x in info['slices']}
    return fr, sl


def write(tags, runtime_dt=None):
    for tag in tags:
        for case in sorted(glob.glob(os.path.join(CFD, 'runs', tag, '*_*'))):
            if not runtime_dt and not os.path.isfile(os.path.join(case, 'solve.done')):
                continue
            cj = json.load(open(os.path.join(case, 'case.json')))
            fr, sl = frame(cj['id'])
            sets = []
            for d in DEPTHS:
                A, U, V, w = gate_points(fr, sl[d]['real'])
                P = world(fr, sl[d]['real'], d, A, U, V)
                name = 'g' + str(d).replace('.', 'p')
                sets.append(f'        {name} {{ type cloud; axis xyz; points (\n'
                            + ''.join('            (%.9g %.9g %.9g)\n' % tuple(p) for p in P) + '        ); }\n')
            body = ("echoGate\n{\n    type sets; libs (sampling); setFormat raw; interpolationScheme cellPoint; fields (U);\n"
                    + (f"    writeControl adjustableRunTime; writeInterval {runtime_dt:g};\n" if runtime_dt else '')
                    + "    sets\n    {\n" + ''.join(sets) + "    }\n}\n")
            if runtime_dt:      # included by the controlDict of a pulsatile run (sampled during the run)
                open(os.path.join(case, 'system', 'echoGateFO'), 'w').write(body)
            else:               # post-processing of a finished run
                open(os.path.join(case, 'system', 'sampleEchoGate'), 'w').write(
                    "FoamFile { version 2.0; format ascii; class dictionary; object sampleEchoGate; }\nfunctions\n{\n" + body + "}\n")
        print('wrote system/sampleEchoGate for', tag)


def spectrum(vb, w):
    """Histogram of beam velocity (m/s) -> bin centres (m/s), power, modal, max (pct), mean."""
    b = BIN_CM_S / 100
    edges = np.arange(np.floor(vb.min() / b) * b, vb.max() + 2 * b, b)
    h, e = np.histogram(vb, bins=edges, weights=w)
    hs = np.convolve(h, np.ones(3) / 3, mode='same')
    cen = 0.5 * (e[1:] + e[:-1])
    o = np.argsort(vb); cw = np.cumsum(w[o]) / w.sum()
    vmax = vb[o][np.searchsorted(cw, MAX_PCT / 100)]
    return cen, h / h.max(), cen[np.argmax(hs)], vmax, np.sum(w * vb) / w.sum()


def load_gate(case, d):
    f = sorted(glob.glob(os.path.join(case, 'postProcessing', 'echoGate', '*', 'g' + str(d).replace('.', 'p') + '_U.xy')))
    if not f:
        return None
    return np.loadtxt(f[-1])


def analyse(tags):
    import pandas as pd
    sel = pd.read_csv(os.path.join(CFD, 'case_selection.csv'))[['id', 'pair', 'arm']]
    allrows = []
    for tag in tags:
        res = pd.read_csv(os.path.join(CFD, f'cfd_results_{tag}.csv')).set_index(['id', 'twin', 'depth_mm'])
        rows = []
        for case in sorted(glob.glob(os.path.join(CFD, 'runs', tag, '*_*'))):
            cj = json.load(open(os.path.join(case, 'case.json')))
            fr, sl = frame(cj['id'])
            for d in DEPTHS:
                D = load_gate(case, d)
                if D is None:
                    continue
                A, U, V, w = gate_points(fr, sl[d]['real'])
                P = world(fr, sl[d]['real'], d, A, U, V)
                # match sampled points back to the gate grid (points outside the mesh are dropped by OpenFOAM)
                from scipy.spatial import cKDTree
                dist, idx = cKDTree(P).query(D[:, :3])
                ok = dist < 1e-6
                Uw, ww = D[ok, 3:6], w[idx[ok]]
                r0 = res.loc[(cj['id'], cj['twin'], d)]
                vmean = r0['v_mean']
                for ang in ANGLES:
                    a = np.radians(ang); b = np.cos(a) * fr['n'] + np.sin(a) * fr['e1']
                    cen, pw, modal, vmax, mean = spectrum(Uw @ b, ww)
                    for ro, v in (('mean', mean), ('modal', modal), ('max', vmax)):
                        rows.append(dict(tag=tag, id=cj['id'], twin=cj['twin'], depth_mm=d, angle_deg=ang, readout=ro,
                                         v=v, v_mean=vmean, profile_ratio=v / vmean, area_ratio=r0['area_ratio'],
                                         SV_ratio=r0['area_ratio'] * v / vmean, frac_in_mesh=ok.mean(),
                                         ring_profile_ratio=r0['profile_ratio']))
        df = pd.DataFrame(rows).merge(sel, on='id')
        df.to_csv(os.path.join(CFD, f'echo_spectrum_{tag}.csv'), index=False)
        print(f'wrote echo_spectrum_{tag}.csv ({len(df)} rows); points inside mesh: min {df.frac_in_mesh.min():.3f}')
        allrows.append(df)
    df = pd.concat(allrows)
    z = df[df.angle_deg == 0]
    summ = z.groupby(['tag', 'arm', 'twin', 'depth_mm', 'readout'])[['profile_ratio', 'SV_ratio']].median().round(3).unstack('readout')
    summ.to_csv(os.path.join(CFD, 'echo_spectrum_summary.csv'))
    pd.set_option('display.width', 220)
    print(summ.to_string())
    ang = df[(df.twin == 'real')].groupby(['tag', 'readout', 'angle_deg'])['SV_ratio'].median().unstack('angle_deg').round(3)
    print('\nreal twin, median SV_ratio over scans and depths, by beam angle (deg):'); print(ang.to_string())
    return df


def figure(df):
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(15, 9.5))
    gs = fig.add_gridspec(2, 3, height_ratios=[1, 1.05], hspace=0.42, wspace=0.28)
    # top: example spectra (real twin, 7.5 mm, LV inflow): three cases, 142 (diastolic), 388 and 33 (both systolic)
    for k, cid in enumerate((142, 388, 33)):
        ax = fig.add_subplot(gs[0, k])
        case = os.path.join(CFD, 'runs', 'lv_lam_h1', f'{cid}_real')
        fr, sl = frame(cid); d = 7.5
        D = load_gate(case, d); A, U, V, w = gate_points(fr, sl[d]['real']); P = world(fr, sl[d]['real'], d, A, U, V)
        from scipy.spatial import cKDTree
        dist, idx = cKDTree(P).query(D[:, :3]); ok = dist < 1e-6
        vmean = df[(df.tag == 'lv_lam_h1') & (df.id == cid) & (df.twin == 'real') & (df.depth_mm == d)].v_mean.iloc[0]
        cen, pw, modal, vmax, mean = spectrum(D[ok, 3:6] @ fr['n'], w[idx[ok]])
        ax.fill_between(cen / vmean, pw, step='mid', color='#9db4cc', alpha=0.8, lw=0)
        for x, c, lab in ((1.0, 'k', 'true mean (Q/A)'), (mean / vmean, '#2a78d6', 'gate mean'), (modal / vmean, '#d1495b', 'modal'),
                          (vmax / vmean, '#7a3db8', f'max ({MAX_PCT:g}th pct)')):
            ax.axvline(x, color=c, lw=1.6, ls='--' if c == 'k' else '-', label=f'{lab}: {x:.2f}')
        ax.set_xlim(0.9, 1.75); ax.set_ylim(0, 1.6); ax.set_yticks([])
        ax.set_xlabel('beam velocity ÷ true mean velocity'); ax.legend(fontsize=7.5, frameon=False, loc='upper right')
        arm = 'diastolic' if cid == 142 else 'systolic'
        ax.set_title(f'{chr(65 + k)}  Case {cid} ({arm}), real LVOT, LV inflow, 7.5 mm', fontsize=9.5, loc='left')
        for sp in ('top', 'right', 'left'):
            ax.spines[sp].set_visible(False)
    # bottom: SV_ratio by read-out, real twin, angle 0, per tag x arm (dots = scans, bars = median), all depths pooled
    z = df[(df.angle_deg == 0) & (df.twin == 'real')]
    cols = {'mean': '#2a78d6', 'modal': '#d1495b', 'max': '#7a3db8'}
    for k, (tag, ttl) in enumerate((('lam_h1', 'Straight-tube inflow'), ('lv_lam_h1', 'Patient LV inflow'))):
        ax = fig.add_subplot(gs[1, k])
        xt, xl = [], []
        for i, arm in enumerate(('diastolic', 'systolic')):
            for j, ro in enumerate(('mean', 'modal', 'max')):
                x = i * 4 + j
                y = z[(z.tag == tag) & (z.arm == arm) & (z.readout == ro)].SV_ratio
                ax.bar(x, y.median(), color=cols[ro], alpha=0.35, width=0.8)
                ax.plot(x + np.random.default_rng(0).uniform(-0.25, 0.25, len(y)), y, 'o', ms=3, color=cols[ro])
                xt.append(x); xl.append(f'{ro}\n{arm}' if j == 1 else ro)
        ax.axhline(1, color='k', lw=0.8)
        ax.set_xticks(xt, xl, fontsize=8); ax.set_ylim(0.4, 1.6)
        ax.set_ylabel('echo SV ÷ true SV (real LVOT)')
        ax.set_title(f'{chr(68 + k)}  {ttl}: echo SV by read-out\n(real LVOT; dots = scans x depths)', fontsize=9.5, loc='left')
        for sp in ('top', 'right'):
            ax.spines[sp].set_visible(False)
    # Huh 2022 reference panel
    ax = fig.add_subplot(gs[1, 2])
    zz = df[(df.angle_deg == 0) & (df.twin == 'real') & (df.tag == 'lv_lam_h1')]
    lab, area, vti, sv = [], [], [], []
    for ro in ('mean', 'modal', 'max'):
        q = zz[zz.readout == ro]
        lab.append(f'CFD LV, {ro}'); area.append(100 * (q.area_ratio.median() - 1)); vti.append(100 * (q.profile_ratio.median() - 1)); sv.append(100 * (q.SV_ratio.median() - 1))
    lab.append('Huh 2022\n(echo vs 4D flow)'); area.append(-25.5); vti.append(29.4); sv.append(3.8)
    y = np.arange(len(lab))[::-1]
    for off, vals, c, nm in ((-0.25, area, '#8a8a8a', 'area'), (0, vti, '#e08a2c', 'velocity / VTI'), (0.25, sv, '#2a9d4b', 'SV')):
        ax.barh(y + off, vals, height=0.24, color=c, label=nm)
    ax.axvline(0, color='k', lw=0.8); ax.set_yticks(y, lab, fontsize=8); ax.set_xlabel('% error of echo vs reference (median)')
    ax.legend(fontsize=8, frameon=False, loc='upper left', bbox_to_anchor=(0.0, -0.18), ncol=3)
    ax.set_title('F  Error split (median) vs Huh et al. 2022', fontsize=9.5, loc='left')
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)
    fig.suptitle(f'Emulated PW-Doppler spectral read-out (gate {GATE_LEN_MM:g} mm long, beam FWHM {BEAM_FWHM_MM:g} mm, at the LVOT centroid, beam along the LVOT axis)', fontsize=11.5)
    fig.text(0.5, 0.01, 'Spectrum = volume-weighted histogram of beam-direction velocity in the gate (bins 2 cm/s). Steady peak systole (Q = 400 mL/s): '
             'a peak-velocity ratio, not a full-systole VTI.\nHuh et al., Int J Cardiovasc Imaging 2022;38:1815–23: area −24.5 to −26.5%, VTI +28.5 to +30.2%, SV +3.8%.',
             ha='center', fontsize=8)
    out = os.path.join(CFD, 'qc', 'echo_spectrum_readout.png'); fig.savefig(out, dpi=135, bbox_inches='tight'); plt.close(fig)
    print('wrote', os.path.relpath(out, CFD))


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('step', choices=['write', 'analyse'])
    ap.add_argument('--tags', default='lam_h1,lv_lam_h1')
    ap.add_argument('--runtime', type=float, default=None, help='write system/echoGateFO sampled every RUNTIME s during a (pulsatile) run')
    a = ap.parse_args(); tags = a.tags.split(',')
    if a.step == 'write':
        write(tags, a.runtime)
    else:
        figure(analyse(tags))
