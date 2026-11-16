"""Pathlines from the tracer-particle run (make_particle_case.py) and a 3-D figure.

The vtkCloud snapshots (postProcessing/tracerPositions/*.vtp, every ~1 ms) store each particle's position and velocity
but no ID, and particle order can change between snapshots in a parallel run. So particles are linked from one snapshot
to the next by predicting x + U*dt and solving a one-to-one assignment (Hungarian algorithm); a link is accepted only if
the predicted and found positions are within MAX_JUMP_MM. Particles that leave through the outlet simply end.

Step 1 (build the pathlines, uv Python):
  uv run --python 3.12 --with numpy --with scipy --with vtk python -I cfd/openfoam/render_pathlines.py build
  -> <case>/postProcessing/pathlines.vtp (lines with speed_rel, t_ms, seed_depth_mm), summary printed
Step 2 (render with ParaView):
  /Applications/ParaView-6.2.0.app/Contents/bin/pvbatch cfd/openfoam/render_pathlines.py render
Step 3 (montage):
  uv run --python 3.12 --with numpy --with matplotlib python -I cfd/openfoam/render_pathlines.py montage
  -> cfd/qc/pathlines_<tag>_<case>.png
"""
import glob, json, os, sys

CFD = os.path.join(os.getcwd(), 'cfd')
TAG, CASE, STL = 'lv_apical_lam_tr_h1_particles', '388_real', 'real_lv_apical.stl'
MAX_JUMP_MM = 0.5
VMAX = 1.6
CDIR = os.path.join(CFD, 'runs', TAG, CASE)
OUT_VTP = os.path.join(CDIR, 'postProcessing', 'pathlines.vtp')
PNGS = [os.path.join(CFD, 'qc', f'pathlines_{TAG}_{CASE}_{k}.png') for k in ('speed', 'time', 'seed')]


def frame_vmean():
    cj = json.load(open(os.path.join(CDIR, 'case.json')))
    info = json.load(open(os.path.join(CFD, 'geometry', str(cj['id']), 'slices.json')))
    s = {round(x['depth_mm'], 2): x for x in info['slices']}[7.5]
    return cj, info['frame'], cj['Q_m3s'] / (s['real']['area_mm2'] * 1e-6)


def build():
    import numpy as np, vtk
    from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
    from scipy.optimize import linear_sum_assignment
    cj, fr, vmean = frame_vmean()
    fr = {k: np.array(v) for k, v in fr.items()}
    series = json.load(open(os.path.join(CDIR, 'postProcessing', 'tracerPositions', 'kinematicCloud.vtp.series')))['files']
    snaps = []
    r = vtk.vtkXMLPolyDataReader()
    for f in series:
        r.SetFileName(os.path.join(CDIR, 'postProcessing', 'tracerPositions', f['name'])); r.Update(); d = r.GetOutput()
        if d.GetNumberOfPoints() == 0:
            snaps.append((f['time'], np.zeros((0, 3)), np.zeros((0, 3)))); continue
        x = vtk_to_numpy(d.GetPoints().GetData()).astype(float) * 1e3                  # mm
        u = vtk_to_numpy(d.GetPointData().GetArray('U')).astype(float)                # m/s
        snaps.append((f['time'], x, u))
    t0 = snaps[0][0]
    # tracks: list of lists of (t, x, u); 'open' maps index in the current snapshot -> track id
    tracks = [[(snaps[0][0], snaps[0][1][i], snaps[0][2][i])] for i in range(len(snaps[0][1]))]
    cur = list(range(len(tracks)))
    jumps, n_links, n_lost = [], 0, 0
    for (ta, xa, ua), (tb, xb, ub) in zip(snaps[:-1], snaps[1:]):
        if len(xa) == 0 or len(xb) == 0:
            cur = []; continue
        dt = tb - ta
        pred = xa + 0.5 * (ua + 0) * dt * 1e3                                     # first guess
        # refine the prediction with the velocity at the found point (two-pass midpoint)
        D = np.linalg.norm(pred[:, None, :] - xb[None, :, :], axis=2)
        ia, ib = linear_sum_assignment(np.where(D < 5, D, 1e6))
        pred2 = pred.copy(); pred2[ia] = xa[ia] + 0.5 * (ua[ia] + ub[ib]) * dt * 1e3
        D = np.linalg.norm(pred2[:, None, :] - xb[None, :, :], axis=2)
        ia, ib = linear_sum_assignment(np.where(D < 5, D, 1e6))
        ok = D[ia, ib] < MAX_JUMP_MM
        nxt = [None] * len(xb)
        for a_, b_, good in zip(ia, ib, ok):
            if good and a_ < len(cur) and cur[a_] is not None:
                tracks[cur[a_]].append((tb, xb[b_], ub[b_])); nxt[b_] = cur[a_]; jumps.append(D[a_, b_]); n_links += 1
            else:
                n_lost += 1
        for b_ in range(len(xb)):                    # an unmatched point starts a new (short) track
            if nxt[b_] is None:
                tracks.append([(tb, xb[b_], ub[b_])]); nxt[b_] = len(tracks) - 1
        cur = nxt
    jumps = np.array(jumps)
    print(f'{len(snaps)} snapshots, {t0:.4f}-{snaps[-1][0]:.4f} s; {len(snaps[0][1])} particles at start, {len(snaps[-1][1])} at end')
    print(f'links accepted {n_links}, rejected {n_lost}; prediction error median {np.median(jumps):.3f} mm, 99th pct {np.percentile(jumps, 99):.3f} mm')
    # summary: seed depth vs ejected (left through the valve plane) by the end
    full = [tr for tr in tracks if abs(tr[0][0] - t0) < 1e-9]
    seed_d = np.array([-(tr[0][1] - fr['c']) @ fr['n'] for tr in full])
    end_d = np.array([-(tr[-1][1] - fr['c']) @ fr['n'] for tr in full])
    ejected = end_d < 0                                        # last seen above the valve plane, or left through the outlet
    t_valve = []
    for tr in full:
        dd = np.array([-(p[1] - fr['c']) @ fr['n'] for p in tr]); k = np.where(dd < 0)[0]
        t_valve.append((tr[k[0]][0] - t0) * 1e3 if len(k) else np.nan)
    t_valve = np.array(t_valve)
    print('seed depth band (mm)  n   ejected by end   median time to valve (ms)')
    for lo, hi in ((2, 10), (10, 20), (20, 30), (30, 45)):
        k = (seed_d >= lo) & (seed_d < hi)
        print(f'  {lo:2d}-{hi:2d}              {k.sum():4d}   {100 * ejected[k].mean():5.1f}%          {np.nanmedian(t_valve[k]) if np.isfinite(t_valve[k]).any() else float("nan"):6.1f}')
    # write lines
    pts, spd, tms, sdep, lines = [], [], [], [], vtk.vtkCellArray()
    for tr in tracks:
        if len(tr) < 3:
            continue
        d0 = -(tr[0][1] - fr['c']) @ fr['n']
        ids = []
        for (t, x, u) in tr:
            ids.append(len(pts)); pts.append(x); spd.append(np.linalg.norm(u) / vmean); tms.append((t - t0) * 1e3); sdep.append(d0)
        pl = vtk.vtkPolyLine(); pl.GetPointIds().SetNumberOfIds(len(ids))
        for k, i in enumerate(ids):
            pl.GetPointIds().SetId(k, i)
        lines.InsertNextCell(pl)
    P = vtk.vtkPoints(); P.SetData(numpy_to_vtk(np.array(pts), deep=True))
    poly = vtk.vtkPolyData(); poly.SetPoints(P); poly.SetLines(lines)
    for name, arr in (('speed_rel', spd), ('t_ms', tms), ('seed_depth_mm', sdep)):
        a = numpy_to_vtk(np.array(arr, float), deep=True); a.SetName(name); poly.GetPointData().AddArray(a)
    w = vtk.vtkXMLPolyDataWriter(); w.SetFileName(OUT_VTP); w.SetInputData(poly); w.Write()
    json.dump(dict(t0=t0, t_end=snaps[-1][0], n_tracks=lines.GetNumberOfCells(), link_err_median_mm=float(np.median(jumps)),
                   link_err_p99_mm=float(np.percentile(jumps, 99)), links_rejected=n_lost,
                   ejected_pct_by_seed_band={f'{lo}-{hi}': float(100 * ejected[(seed_d >= lo) & (seed_d < hi)].mean())
                                             for lo, hi in ((2, 10), (10, 20), (20, 30), (30, 45))}),
              open(os.path.join(CDIR, 'postProcessing', 'pathlines_summary.json'), 'w'), indent=1)
    print('wrote', os.path.relpath(OUT_VTP, CFD))


def render():
    from paraview.simple import (XMLPolyDataReader, STLReader, Clip, Show, CreateView, ColorBy, GetColorTransferFunction,
                                 SaveScreenshot, Delete, GetSources)
    cj, fr, vmean = frame_vmean()
    top = [fr['c'][i] + 12 * fr['n'][i] for i in range(3)]
    for (arr, rng, presets), png in zip((('speed_rel', (0, VMAX), ('Inferno', 'Inferno (matplotlib)')),
                                         ('t_ms', (0, 150), ('Viridis', 'Viridis (matplotlib)')),
                                         ('seed_depth_mm', (2, 45), ('Turbo', 'Jet', 'Rainbow Uniform'))), PNGS):
        tr = XMLPolyDataReader(FileName=[OUT_VTP])
        tr = Clip(Input=tr); tr.ClipType = 'Plane'; tr.ClipType.Origin = top; tr.ClipType.Normal = fr['n']; tr.Invert = 1
        wall = STLReader(FileNames=[os.path.join(CFD, 'geometry', str(cj['id']), STL)])
        wall = Clip(Input=wall); wall.ClipType = 'Plane'; wall.ClipType.Origin = top; wall.ClipType.Normal = fr['n']; wall.Invert = 1
        v = CreateView('RenderView'); v.ViewSize = [900, 1100]; v.Background = [1, 1, 1]; v.OrientationAxesVisibility = 0
        try:
            v.UseColorPaletteForBackground = 0
        except Exception:
            pass
        dw = Show(wall, v); dw.Opacity = 0.12; dw.DiffuseColor = [0.55, 0.55, 0.6]; dw.AmbientColor = [0.55, 0.55, 0.6]
        dt = Show(tr, v); ColorBy(dt, ('POINTS', arr)); dt.LineWidth = 1.4
        lut = GetColorTransferFunction(arr)
        for name in presets:
            try:
                lut.ApplyPreset(name, True); break
            except RuntimeError:
                pass
        lut.RescaleTransferFunction(*rng); dt.SetScalarBarVisibility(v, False)
        e1, e2, n = fr['e1'], fr['e2'], fr['n']
        cx = [e2[1] * n[2] - e2[2] * n[1], e2[2] * n[0] - e2[0] * n[2], e2[0] * n[1] - e2[1] * n[0]]
        sgn = 1 if sum(cx[i] * e1[i] for i in range(3)) > 0 else -1          # septum left, mitral / LA right
        c = list(fr['c']); cam = v.GetActiveCamera()
        cam.SetFocalPoint(*c); cam.SetPosition(*[c[i] - 300 * sgn * e2[i] for i in range(3)]); cam.SetViewUp(*n)
        v.ResetCamera(False); cam.Zoom(1.15)
        SaveScreenshot(png, v, ImageResolution=[900, 1100], TransparentBackground=0)
        print('wrote', png)
        for src in list(GetSources().values()):
            Delete(src)
        Delete(v)


def montage():
    import numpy as np
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    import matplotlib.image as mi
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import Normalize
    summ = json.load(open(os.path.join(CDIR, 'postProcessing', 'pathlines_summary.json')))
    fig, axs = plt.subplots(1, 3, figsize=(16, 7.4))
    titles = ('A  Speed along each path', 'B  Time since release', 'C  Where each parcel started (depth below the valve)')
    cmaps = (('inferno', 0, VMAX, 'speed ÷ true mean LVOT velocity at 7.5 mm'), ('viridis', 0, 150, 'ms since release'),
             ('turbo', 2, 45, 'release depth (mm below the valve)'))
    for ax, png, ttl, (cm, lo, hi, lab) in zip(axs, PNGS, titles, cmaps):
        im = mi.imread(png)[..., :3]; m = (im < 0.97).any(2); r_, c_ = np.where(m)
        ax.imshow(im[max(r_.min() - 20, 0):r_.max() + 20, max(c_.min() - 20, 0):c_.max() + 20]); ax.set_axis_off()
        ax.set_title(ttl, fontsize=10.5, loc='left')
        ax.text(0.0, 0.5, 'septum', transform=ax.transAxes, rotation=90, va='center', fontsize=8, color='0.4')
        ax.text(1.0, 0.5, 'mitral / LA', transform=ax.transAxes, rotation=90, va='center', ha='right', fontsize=8, color='0.4')
        cax = ax.inset_axes([0.15, -0.07, 0.7, 0.025])
        fig.colorbar(ScalarMappable(Normalize(lo, hi), cm), cax=cax, orientation='horizontal').set_label(lab, fontsize=8.5)
    ej = summ['ejected_pct_by_seed_band']
    fig.suptitle(f'Case 388, apical inflow, unsteady flow at peak systole: pathlines of {summ["n_tracks"]} tracer parcels over '
                 f'{1e3 * (summ["t_end"] - summ["t0"]):.0f} ms', fontsize=12)
    fig.text(0.5, 0.02, 'Parcels released at random positions 2–45 mm below the valve at t = 0. Share ejected through the valve by the end: '
             + ', '.join(f'{k} mm: {v:.0f}%' for k, v in ej.items())
             + '.\nConstant peak flow (400 mL/s), ejection only: no diastolic filling or vortex. Linked between 1 ms snapshots '
             f'(median prediction error {summ["link_err_median_mm"]:.3f} mm).', ha='center', fontsize=8.5)
    fig.subplots_adjust(left=0.02, right=0.98, top=0.9, bottom=0.16, wspace=0.05)
    out = os.path.join(CFD, 'qc', f'pathlines_{TAG}_{CASE}.png'); fig.savefig(out, dpi=130); plt.close(fig)
    print('wrote', os.path.relpath(out, CFD))


if __name__ == '__main__':
    {'build': build, 'render': render, 'montage': montage}[sys.argv[1]]()
