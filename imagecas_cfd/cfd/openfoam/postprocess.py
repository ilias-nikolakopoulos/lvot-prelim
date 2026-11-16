"""Echo stroke-volume error from the twin CFD runs.

For each case and echo depth d:
  A_true   twin slice area (slices.json)
  A_echo   pi*AP^2/4 (real: AP caliper through the centroid; circle: its diameter, so A_echo = A_true)
  v_mean   Q / A_true (true mean velocity)
  v_sv     PW-Doppler sample volume: mean of U.n over the centroid + 6 points on a 1.5 mm ring
  v_c      U.n at the centroid alone
  area_ratio    = A_echo / A_true
  profile_ratio = v_sv / v_mean
  SV_ratio      = area_ratio * profile_ratio   (echo SV / true SV; 1 = no error)
  dp_Pa    pressure drop between the inlet-extension probe and the aortic probe (rho = 1060 kg/m^3)
  drift_pct     change of v_sv over the last 50 iterations (convergence check)
Transient runs (case.json transient=true, pimpleFoam): v_sv, v_c and the ratios are time means over
t >= t_avg_s; v_sv_std_pct and v_sv_p2p_pct give the fluctuation (std and max-min, % of v_mean),
and drift_pct is the change of the running mean over the last half of the averaging window.

Run from the imagecas folder:
  uv run --python 3.12 --with numpy --with pandas python -I cfd/openfoam/postprocess.py --tag lam_h1
"""
import argparse, glob, json, os, re
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
CFD = os.path.dirname(HERE)
RHO = 1060.0
NUM = r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?'


def read_probe(case, fo, field):
    rows = {}
    for f in sorted(glob.glob(os.path.join(case, 'postProcessing', fo, '*', field))):
        for line in open(f):
            if line.startswith('#') or not line.strip():
                continue
            t, rest = line.split(None, 1)
            vals = [float(x) for x in re.findall(NUM, rest)]
            rows[float(t)] = vals
    if not rows:
        return None, None
    it = np.array(sorted(rows))
    return it, np.array([rows[t] for t in it])


def one_case(case, sel):
    cj = json.load(open(os.path.join(case, 'case.json')))
    info = json.load(open(os.path.join(CFD, 'geometry', str(cj['id']), 'slices.json')))
    sl = {round(s['depth_mm'], 2): s for s in info['slices']}
    n = np.array(cj['n'])
    it, U = read_probe(case, 'probesU', 'U')
    itp, P = read_probe(case, 'probesP', 'p')
    out = []
    if it is None:
        return out
    U = U.reshape(len(it), -1, 3)
    un = U @ n
    bad = np.abs(U).max(axis=(0, 2)) > 1e10            # probe outside the mesh
    meta = cj['probe_meta']
    last = it[-1]
    k200 = np.searchsorted(it, last - 50)
    dp = (P[-1, 0] - P[-1, 1]) * RHO if P is not None else np.nan
    tr = bool(cj.get('transient'))
    w = it >= cj['t_avg_s'] if tr else None
    for d in cj['echo_depths_mm']:
        idx = [i for i, m in enumerate(meta) if m['depth_mm'] == d and not bad[i]]
        c_idx = [i for i in idx if meta[i]['k'] == 0]
        s = sl[round(d, 2)]
        if cj['twin'] == 'real':
            A, ap = s['real']['area_mm2'], s['real']['AP_mm']
        else:
            A = s['circle']['area_mm2']; ap = 2 * np.sqrt(A / np.pi)
        A_echo = np.pi * ap ** 2 / 4
        v_mean = cj['Q_m3s'] / (A * 1e-6)
        v_sv_hist = un[:, idx].mean(axis=1)
        v_sv = v_sv_hist[-1]
        extra = {}
        if tr:
            h = v_sv_hist[w]; v_sv = h.mean()
            half = h[: len(h) // 2].mean()
            extra = dict(v_sv_std_pct=100 * h.std() / v_mean, v_sv_p2p_pct=100 * (h.max() - h.min()) / v_mean,
                         t_avg_window_s=float(it[w][-1] - it[w][0]))
        out.append(dict(id=cj['id'], twin=cj['twin'], tag=cj['tag'], depth_mm=d, iterations=int(last),
                        A_true_mm2=A, A_echo_mm2=A_echo, v_mean=v_mean, v_sv=v_sv,
                        v_c=(un[w, c_idx[0]].mean() if tr else un[-1, c_idx[0]]) if c_idx else np.nan, n_sv_points=len(idx),
                        area_ratio=A_echo / A, profile_ratio=v_sv / v_mean, SV_ratio=A_echo / A * v_sv / v_mean,
                        dp_Pa=dp, drift_pct=100 * (v_sv - half) / v_mean if tr else 100 * (v_sv - v_sv_hist[k200]) / v_mean,
                        **extra))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--tag', default='lam_h1')
    args = ap.parse_args()
    sel = pd.read_csv(os.path.join(CFD, 'case_selection.csv'))
    rows = []
    for case in sorted(glob.glob(os.path.join(CFD, 'runs', args.tag, '*_*'))):
        rows += one_case(case, sel)
    if not rows:
        print('no probe output yet'); return
    df = pd.DataFrame(rows).merge(sel[['id', 'pair', 'arm']], on='id', how='left')
    df = df.sort_values(['pair', 'arm', 'id', 'twin', 'depth_mm'])
    out = os.path.join(CFD, f'cfd_results_{args.tag}.csv')
    df.to_csv(out, index=False)
    print(f'wrote {os.path.relpath(out, CFD)}  ({df.id.nunique()} ids, {len(df)} rows)')
    pd.set_option('display.width', 200)
    piv = df.pivot_table(index=['pair', 'arm', 'id', 'depth_mm'], columns='twin',
                         values=['profile_ratio', 'SV_ratio'], aggfunc='first').round(3)
    print(piv.to_string())
    print('\nmedian SV_ratio by arm, twin and depth:')
    print(df.groupby(['arm', 'twin', 'depth_mm']).SV_ratio.median().unstack().round(3).to_string())
    print('\nmax |drift|: %.2f%% of v_mean' % df.drift_pct.abs().max())
    if 'v_sv_std_pct' in df:
        print('\nsample-volume velocity fluctuation (% of v_mean):')
        print(df[['id', 'twin', 'depth_mm', 'v_sv_std_pct', 'v_sv_p2p_pct', 't_avg_window_s']].round(2).to_string(index=False))


if __name__ == '__main__':
    main()
