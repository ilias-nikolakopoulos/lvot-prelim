"""Velocity across the LVOT cross-section: peak-velocity displacement + figures.

Reads the echo-plane velocity (postProcessing/echoPlanes/<latest>/U_d*.raw) of every case in runs/<tag>,
projects it onto the LVOT frame (u along e1 = AP toward the LA, v along e2 = mediolateral, w = U.n = axial),
and interpolates w onto a 0.1 mm grid.

Peak displacement = position of the maximum axial velocity relative to the slice centroid
(the echo sample-volume position). Reported in mm, split into AP and ML, and as a fraction of the
equivalent radius sqrt(A/pi).

Outputs:
  cfd/peak_displacement_<tag>.csv
  cfd/qc/velocity_sections_<tag>_d<depth>.png    all 12 cases, real vs circle, at one depth
  cfd/qc/velocity_sections_<tag>_<id>.png        one case, all depths, real vs circle

Run from the imagecas folder:
  uv run --python 3.12 --with numpy --with scipy --with pandas --with matplotlib python -I cfd/openfoam/plane_figures.py --tag lam_h1
"""
import argparse, glob, json, os
import numpy as np
import pandas as pd
from scipy.interpolate import griddata
from scipy.ndimage import gaussian_filter

HERE = os.path.dirname(os.path.abspath(__file__))
CFD = os.path.dirname(HERE)
H = 0.1          # mm grid
MEANC = '#ff2bd6'  # magenta: the mean-velocity isoline (w = v_mean)
FIELD = 'U'      # 'UMean' for the time-averaged velocity of transient runs (--field UMean)


def load_plane(case, depth):
    tdirs = sorted(glob.glob(os.path.join(case, 'postProcessing', 'echoPlanes', '*')), key=lambda p: float(os.path.basename(p)))
    tag = 'd' + str(depth).replace('.', 'p')
    f = os.path.join(tdirs[-1], f'{FIELD}_{tag}.raw')
    return np.loadtxt(f, comments='#')


def section(case, depth, info, cj):
    fr = {k: np.array(v) for k, v in info['frame'].items()}
    D = load_plane(case, depth)
    rel = D[:, :3] * 1e3 - fr['c']
    u, v = rel @ fr['e1'], rel @ fr['e2']
    w = D[:, 3:6] @ fr['n']
    t1, t2 = D[:, 3:6] @ fr['e1'], D[:, 3:6] @ fr['e2']      # in-plane (secondary) components
    s = {round(x['depth_mm'], 2): x for x in info['slices']}[round(depth, 2)]
    cu, cv = s['real']['cu'], s['real']['cv']
    A = s['real' if cj['twin'] == 'real' else 'circle']['area_mm2']
    gu = np.arange(np.floor(u.min()), np.ceil(u.max()) + H, H); gv = np.arange(np.floor(v.min()), np.ceil(v.max()) + H, H)
    U_, V_ = np.meshgrid(gu, gv, indexing='ij')
    W = griddata((u, v), w, (U_, V_), method='linear')
    T1 = griddata((u, v), t1, (U_, V_), method='linear'); T2 = griddata((u, v), t2, (U_, V_), method='linear')
    vmean = cj['Q_m3s'] / (A * 1e-6)
    # peak of a lightly smoothed field (0.3 mm) to avoid single-face noise
    Wf = np.where(np.isnan(W), 0, W); M = (~np.isnan(W)).astype(float)
    Ws = gaussian_filter(Wf, 3) / np.maximum(gaussian_filter(M, 3), 1e-9); Ws[np.isnan(W)] = np.nan
    i, j = np.unravel_index(np.nanargmax(Ws), Ws.shape)
    pu, pv = gu[i], gv[j]
    req = np.sqrt(A / np.pi)
    ok = ~np.isnan(W); wpos = np.clip(W[ok], 0, None)
    gcu, gcv = U_[ok].mean(), V_[ok].mean()                                   # geometric centroid of the section
    fu, fv = (U_[ok] * wpos).sum() / wpos.sum(), (V_[ok] * wpos).sum() / wpos.sum()   # flux-weighted flow centre
    # 4D-flow MRI style metrics (Sigovan et al., JMRI 2011): normalized flow displacement and flow angle
    deq = 2 * np.sqrt(A / np.pi)
    mt1, mt2 = (np.nan_to_num(T1[ok]) * wpos).sum() / wpos.sum(), (np.nan_to_num(T2[ok]) * wpos).sum() / wpos.sum()
    mw = (W[ok] * wpos).sum() / wpos.sum()
    flow_angle = float(np.degrees(np.arctan2(np.hypot(mt1, mt2), mw)))
    sec = float(np.nanmean(np.hypot(T1[ok], T2[ok])) / vmean)
    return dict(u=gu, v=gv, W=W / vmean, cu=cu, cv=cv, pu=pu, pv=pv, fu=fu, fv=fv, req=req,
                row=dict(id=cj['id'], twin=cj['twin'], depth_mm=depth, peak_over_mean=np.nanmax(Ws) / vmean,
                         disp_mm=float(np.hypot(pu - cu, pv - cv)), disp_AP_mm=float(pu - cu), disp_ML_mm=float(pv - cv),
                         disp_frac_req=float(np.hypot(pu - cu, pv - cv) / req), req_mm=req,
                         flow_centre_off_mm=float(np.hypot(fu - gcu, fv - gcv)), flow_centre_AP_mm=float(fu - gcu),
                         flow_centre_ML_mm=float(fv - gcv), flow_displacement=float(np.hypot(fu - gcu, fv - gcv) / deq),
                         flow_angle_deg=flow_angle, flow_angle_dir_AP=float(mt1), flow_angle_dir_ML=float(mt2), secondary_over_mean=sec, area_frac_above_1p1=float((W[ok] / vmean > 1.1).mean()),
                         core_cv_pct=float(100 * np.std(W[ok][W[ok] / vmean > 0.9]) / np.mean(W[ok][W[ok] / vmean > 0.9])),
                         sv_centroid_offset_mm=float(np.hypot(cu - gcu, cv - gcv))))


def draw(ax, S, vmax, title=None):
    import matplotlib.pyplot as plt
    im = ax.pcolormesh(S['v'], S['u'], S['W'], cmap='viridis', vmin=0, vmax=vmax, shading='auto')
    Wc = np.nan_to_num(S['W'], nan=-1)
    ax.contour(S['v'], S['u'], Wc, levels=[0.25, 0.5, 0.75, 1.25], colors='w', linewidths=0.4, alpha=0.6)
    ax.contour(S['v'], S['u'], Wc, levels=[1.0], colors=MEANC, linewidths=1.3)                # w = v_mean (VTI-correct)
    ax.add_patch(plt.Circle((S['cv'], S['cu']), 1.5, fill=False, color='w', lw=1.0, ls='--'))   # Doppler sample volume
    ax.plot(S['cv'], S['cu'], '+', color='w', ms=7, mew=1.2)
    ax.plot(S['pv'], S['pu'], 'x', color='#ff5a5f', ms=6, mew=1.5)
    ax.plot(S['fv'], S['fu'], 'o', mfc='none', mec='#ffb000', ms=6, mew=1.5)                      # flux-weighted flow centre
    ax.set_aspect('equal'); ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_visible(False)
    if title: ax.set_title(title, fontsize=8)
    return im


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--tag', default='lam_h1'); ap.add_argument('--one', default='978')
    ap.add_argument('--field', default='U', help='U, or UMean for the time average of a transient run')
    args = ap.parse_args()
    global FIELD; FIELD = args.field
    name = args.tag + ('_mean' if args.field != 'U' else '')
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    sel = pd.read_csv(os.path.join(CFD, 'case_selection.csv'))
    secs, rows = {}, []
    for case in sorted(glob.glob(os.path.join(CFD, 'runs', args.tag, '*_*'))):
        cj = json.load(open(os.path.join(case, 'case.json')))
        if not os.path.isdir(os.path.join(case, 'postProcessing', 'echoPlanes')):
            continue
        info = json.load(open(os.path.join(CFD, 'geometry', str(cj['id']), 'slices.json')))
        for d in cj['echo_depths_mm']:
            S = section(case, d, info, cj); secs[(cj['id'], cj['twin'], d)] = S; rows.append(S['row'])
    df = pd.DataFrame(rows).merge(sel[['id', 'pair', 'arm']], on='id').sort_values(['pair', 'arm', 'twin', 'depth_mm'])
    out = os.path.join(CFD, f'peak_displacement_{name}.csv'); df.to_csv(out, index=False)
    print(f'wrote {os.path.relpath(out, CFD)}')
    print(df.groupby(['twin', 'arm', 'depth_mm'])[['peak_over_mean', 'flow_centre_off_mm', 'flow_displacement', 'flow_angle_deg', 'secondary_over_mean', 'core_cv_pct']].median().round(2).to_string())
    os.makedirs(os.path.join(CFD, 'qc'), exist_ok=True)
    VMAX = 1.6
    note = ('Axial velocity / true mean velocity.  Vertical = AP (up = toward LA, the axis echo measures), horizontal = ML.\n'
            '+ = slice centroid (sample-volume centre), dashed = 1.5 mm sample volume, x = velocity maximum, o = flux-weighted flow centre.\nMagenta line = mean velocity (w = v_mean); white contours 0.25, 0.5, 0.75, 1.25.')
    # 1) all cases at each depth
    pairs = sel.sort_values('pair')
    for d in (5.0, 7.5, 10.0):
        fig, axs = plt.subplots(6, 4, figsize=(8.4, 13.2))
        for r, p in enumerate(sorted(pairs.pair.unique())):
            ids = {a: int(pairs[(pairs.pair == p) & (pairs.arm == a)].id.iloc[0]) for a in ('systolic', 'diastolic')}
            for k, (arm, tw) in enumerate([('systolic', 'real'), ('systolic', 'circle'), ('diastolic', 'real'), ('diastolic', 'circle')]):
                key = (ids[arm], tw, d); ax = axs[r, k]
                if key not in secs: ax.axis('off'); continue
                S = secs[key]; im = draw(ax, S, VMAX)
                ax.set_title(f'{ids[arm]}  max {S["row"]["peak_over_mean"]:.2f}×, centre off {S["row"]["flow_centre_off_mm"]:.1f} mm', fontsize=7)
                if k == 0: ax.set_ylabel(f'pair {p}', fontsize=8)
        for k, t in enumerate(['systolic · real', 'systolic · circle', 'diastolic · real', 'diastolic · circle']):
            fig.text(0.14 + k * 0.215, 0.965, t, ha='center', fontsize=9, weight='bold')
        fig.suptitle(f'LVOT axial velocity at {d:g} mm below the valve ({name})', y=0.995, fontsize=11)
        fig.text(0.5, 0.012, note, ha='center', fontsize=7)
        cb = fig.add_axes([0.93, 0.3, 0.015, 0.4]); fig.colorbar(im, cax=cb).set_label('w / v_mean', fontsize=8)
        fig.subplots_adjust(left=0.04, right=0.91, top=0.95, bottom=0.05, wspace=0.08, hspace=0.25)
        f = os.path.join(CFD, 'qc', f'velocity_sections_{name}_d{str(d).replace(".", "p")}.png'); fig.savefig(f, dpi=150); plt.close(fig)
        print('wrote', os.path.relpath(f, CFD))
    # 2) one case: all depths, real vs circle
    one = int(args.one)
    fig, axs = plt.subplots(3, 3, figsize=(10, 10.5), gridspec_kw=dict(height_ratios=[1, 1, 0.8]))
    for k, d in enumerate((5.0, 7.5, 10.0)):
        for r, tw in enumerate(('real', 'circle')):
            S = secs[(one, tw, d)]; im = draw(axs[r, k], S, VMAX)
            R = S['row']
            axs[r, k].set_title(f'{tw}, {d:g} mm\nmax {R["peak_over_mean"]:.2f}× mean, flow centre {R["flow_centre_off_mm"]:.1f} mm off', fontsize=8)
        ax = axs[2, k]
        for tw, ls in (('real', '-'), ('circle', '--')):
            S = secs[(one, tw, d)]
            iu = np.argmin(np.abs(S['u'] - S['cu'])); iv = np.argmin(np.abs(S['v'] - S['cv']))
            ax.plot(S['u'] - S['cu'], S['W'][:, iv], ls, color='#1f6fb2', lw=1.4, label=f'{tw}: AP')
            ax.plot(S['v'] - S['cv'], S['W'][iu, :], ls, color='#d1495b', lw=1.4, label=f'{tw}: ML')
        ax.axhline(1, color='0.5', lw=0.6); ax.axvspan(-1.5, 1.5, color='0.85', zorder=0)
        ax.set_xlabel('distance from centroid (mm)', fontsize=8); ax.set_ylim(0, 1.6); ax.tick_params(labelsize=7)
        ax.set_title(f'profiles through the centroid, {d:g} mm', fontsize=8)
        if k == 0: ax.set_ylabel('w / v_mean', fontsize=8); ax.legend(fontsize=6.5, frameon=False, loc='lower center', ncol=2)
        for sp in ('top', 'right'): ax.spines[sp].set_visible(False)
    fig.suptitle(f'Case {one}: LVOT axial velocity, real vs same-area circle ({name})', fontsize=11)
    fig.text(0.5, 0.008, note + '\nBottom: AP (blue) and ML (red) lines through the centroid; grey band = sample volume; 1 = true mean.', ha='center', fontsize=7)
    cb = fig.add_axes([0.93, 0.45, 0.015, 0.4]); fig.colorbar(im, cax=cb).set_label('w / v_mean', fontsize=8)
    fig.subplots_adjust(left=0.06, right=0.91, top=0.93, bottom=0.09, wspace=0.12, hspace=0.3)
    f = os.path.join(CFD, 'qc', f'velocity_sections_{name}_{one}.png'); fig.savefig(f, dpi=150); plt.close(fig)
    print('wrote', os.path.relpath(f, CFD))


if __name__ == '__main__':
    main()
