"""Flow pattern on the long-axis (A5C-like) plane: in-plane streamlines + speed, for several inflow variants of one case.

The plane contains the LVOT axis n and the AP axis e1 (septum <-> anterior mitral leaflet / LA) and passes through the
7.5 mm slice centroid, so it cuts the LVOT lengthwise and, in the LV-inflow runs, the LV cavity below it.

Step 1 (writes system/sampleLongAxis in each case):
  uv run --python 3.12 --with numpy python -I cfd/openfoam/long_axis_flow.py write --id 388 --twin real
Step 2 (in the container, for each case; use -fields '(UMean)' for a transient run):
  docker run --rm --user "$(id -u):$(id -g)" -v "$PWD/<case>:/case" -w /case opencfd/openfoam-run:2512 \
      -c "cd /case && postProcess -dict system/sampleLongAxis -latestTime -fields '(U)'"
Step 3 (figure -> cfd/qc/long_axis_flow_<id>_<twin>.png):
  uv run --python 3.12 --with numpy --with scipy --with matplotlib python -I cfd/openfoam/long_axis_flow.py plot --id 388 --twin real
"""
import argparse, glob, json, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CFD = os.path.dirname(HERE)
# (run tag, field, panel title)
RUNS = [('lam_h1', 'U', 'A  Straight tube, plug inflow'),
        ('lv_lam_h1', 'U', 'B  Patient LV cavity, whole-wall inflow'),
        ('lv_apical_lam_tr_h1', 'UMean', 'C  Patient LV cavity, apical inflow\n(transient, time-averaged 0.2–0.5 s)')]
DEPTH = 7.5


def frame(cid):
    info = json.load(open(os.path.join(CFD, 'geometry', str(cid), 'slices.json')))
    fr = {k: np.array(v) for k, v in info['frame'].items()}
    s = {round(x['depth_mm'], 2): x for x in info['slices']}[DEPTH]
    p0 = fr['c'] - DEPTH * fr['n'] + s['real']['cu'] * fr['e1'] + s['real']['cv'] * fr['e2']   # mm
    return fr, p0


def write(cid, twin):
    fr, p0 = frame(cid)
    vec = lambda p: '(%.9g %.9g %.9g)' % tuple(p)
    for tag, field, _ in RUNS:
        case = os.path.join(CFD, 'runs', tag, f'{cid}_{twin}')
        if not os.path.isdir(case):
            print('missing', case); continue
        txt = ("FoamFile { version 2.0; format ascii; class dictionary; object sampleLongAxis; }\n"
               "functions\n{\nsampleLongAxis\n{\n"
               f"  type surfaces; libs (sampling); surfaceFormat raw; fields ({field}); interpolationScheme cellPoint;\n"
               "  surfaces\n  {\n"
               f"    longAxis {{ type cuttingPlane; planeType pointAndNormal; pointAndNormalDict {{ point {vec(p0 * 1e-3)}; normal {vec(fr['e2'])}; }} interpolate false; }}\n"
               "  }\n}\n}\n")
        open(os.path.join(case, 'system', 'sampleLongAxis'), 'w').write(txt)
        print('wrote', os.path.relpath(os.path.join(case, 'system', 'sampleLongAxis'), CFD), field)


def load(case, field):
    d = sorted(glob.glob(os.path.join(case, 'postProcessing', 'sampleLongAxis', '*')), key=lambda p: float(os.path.basename(p)))[-1]
    return np.loadtxt(os.path.join(d, f'{field}_longAxis.raw'), comments='#')


def plot(cid, twin):
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    from scipy.interpolate import griddata
    from scipy.spatial import cKDTree
    fr, p0 = frame(cid)
    Q = json.load(open(os.path.join(CFD, 'runs', 'lam_h1', f'{cid}_{twin}', 'case.json')))['Q_m3s']
    s = {round(x['depth_mm'], 2): x for x in json.load(open(os.path.join(CFD, 'geometry', str(cid), 'slices.json')))['slices']}[DEPTH]
    vmean = Q / (s[twin]['area_mm2'] * 1e-6)
    data = []
    for tag, field, title in RUNS:
        D = load(os.path.join(CFD, 'runs', tag, f'{cid}_{twin}'), field)
        rel = D[:, :3] * 1e3 - fr['c']
        x, depth = rel @ fr['e1'], -(rel @ fr['n'])          # AP position (mm), depth below the valve (mm)
        U = D[:, 3:6]
        data.append((title, x, depth, U @ fr['e1'], -(U @ fr['n']), np.linalg.norm(U, axis=1)))
    ymin = min(d[2].min() for d in data); ymax = max(d[2].max() for d in data)
    xmin = min(d[1].min() for d in data); xmax = max(d[1].max() for d in data)
    ylo, yhi = -12, min(ymax, 46)
    h = 0.25
    gx = np.arange(np.floor(xmin), np.ceil(xmax) + h, h); gy = np.arange(ylo, yhi + h, h)
    GX, GY = np.meshgrid(gx, gy)
    fig, axs = plt.subplots(1, 3, figsize=(15, 8.6), sharey=True)
    for ax, (title, x, y, ux, uy, sp) in zip(axs, data):
        k = (y >= ylo - 2) & (y <= yhi + 2)
        pts = np.c_[x[k], y[k]]
        inside = cKDTree(pts).query(np.c_[GX.ravel(), GY.ravel()])[0].reshape(GX.shape) < 1.1
        S = griddata(pts, sp[k] / vmean, (GX, GY), method='linear'); S[~inside] = np.nan
        UX = griddata(pts, ux[k], (GX, GY), method='linear'); UY = griddata(pts, uy[k], (GX, GY), method='linear')
        UX[~inside] = np.nan; UY[~inside] = np.nan
        im = ax.pcolormesh(gx, gy, S, cmap='magma', vmin=0, vmax=1.6, shading='auto')
        ax.streamplot(gx, gy, np.nan_to_num(UX), np.nan_to_num(UY), density=1.3, color='w', linewidth=0.7, arrowsize=0.8,
                      broken_streamlines=False, minlength=0.2)
        ax.contour(gx, gy, S, levels=[1.0], colors='#ff2bd6', linewidths=1.2)
        for d, ls in ((0, '-'), (5, ':'), (7.5, ':'), (10, ':')):
            ax.axhline(d, color='0.75', lw=0.6, ls=ls)
        ax.set_title(title, fontsize=10.5, loc='left')
        ax.set_aspect('equal'); ax.set_xlim(gx[0], gx[-1]); ax.set_ylim(yhi, ylo)        # apex at the bottom, aorta at the top
        ax.set_xlabel('AP position (mm)   ← septum  |  mitral / LA →')
        ax.set_facecolor('0.92')
    axs[0].set_ylabel('depth below the aortic valve (mm)')
    for d in (0, 5, 7.5, 10):
        axs[0].text(gx[0] + 0.5, d - 0.4, 'valve' if d == 0 else f'{d:g} mm', fontsize=7, color='0.35')
    cb = fig.add_axes([0.3, 0.115, 0.4, 0.018]); fig.colorbar(im, cax=cb, orientation='horizontal').set_label('speed ÷ true mean LVOT velocity at 7.5 mm')
    fig.suptitle(f'Case {cid} {twin} LVOT: flow pattern on the long-axis (A5C-like) plane, steady peak systole, Q = 400 mL/s', fontsize=12)
    fig.text(0.5, 0.012, 'White lines = in-plane streamlines (direction of flow within this plane). Magenta line = speed equal to the LVOT mean velocity. '
             'Grey lines = valve plane and the 5 / 7.5 / 10 mm echo depths.\nFlow goes up (from the LV/inlet toward the aorta). No diastolic vortex: '
             'each run models ejection only, starting from rest; there is no mitral inflow.', ha='center', fontsize=8.5)
    fig.subplots_adjust(left=0.06, right=0.98, top=0.92, bottom=0.2, wspace=0.08)
    out = os.path.join(CFD, 'qc', f'long_axis_flow_{cid}_{twin}.png'); fig.savefig(out, dpi=140); plt.close(fig)
    print('wrote', os.path.relpath(out, CFD))


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('step', choices=['write', 'plot'])
    ap.add_argument('--id', type=int, default=388); ap.add_argument('--twin', default='real')
    a = ap.parse_args()
    write(a.id, a.twin) if a.step == 'write' else plot(a.id, a.twin)
