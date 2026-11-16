"""Render the 3-D streamlines from streamlines_3d.py with ParaView (offscreen), one PNG per case.

Run from the imagecas folder with ParaView's batch Python:
  /Applications/ParaView-6.2.0.app/Contents/bin/pvbatch cfd/openfoam/render_streamlines.py
Then assemble the panels:
  uv run --python 3.12 --with numpy --with matplotlib python -I cfd/openfoam/render_streamlines.py --montage
Output: cfd/qc/streamlines3d_<tag>_<case>.png (per case) and cfd/qc/streamlines3d_388_33.png (montage)
"""
import glob, json, os, sys

CFD = os.path.join(os.getcwd(), 'cfd')
PANELS = [('lam_h1', '388_real', 'U', 'real.stl', 'A  388, straight tube, plug inflow'),
          ('lv_lam_h1', '388_real', 'U', 'real_lv.stl', 'B  388, patient LV, whole-wall inflow'),
          ('lv_apical_lam_tr_h1', '388_real', 'UMean', 'real_lv_apical.stl', 'C  388, patient LV, apical inflow (time-averaged)'),
          ('lv_lam_h1', '33_real', 'U', 'real_lv.stl', 'D  33 (small end-systolic LV), whole-wall inflow')]
VMAX = 1.6


def png(tag, case):
    return os.path.join(CFD, 'qc', f'streamlines3d_{tag}_{case}.png')


def render():
    from paraview.simple import (XMLPolyDataReader, STLReader, Transform, Calculator, Clip, Show, Hide, GetActiveViewOrCreate,
                                 ColorBy, GetColorTransferFunction, SaveScreenshot, Delete, ResetSession,
                                 CreateView, GetSources, Render)
    for tag, case, field, stl, _ in PANELS:
        cdir = os.path.join(CFD, 'runs', tag, case)
        cj = json.load(open(os.path.join(cdir, 'case.json')))
        cid = cj['id']
        info = json.load(open(os.path.join(CFD, 'geometry', str(cid), 'slices.json')))
        fr = info['frame']
        s = {round(x['depth_mm'], 2): x for x in info['slices']}[7.5]
        vmean = cj['Q_m3s'] / (s['real']['area_mm2'] * 1e-6)
        trk = sorted(glob.glob(os.path.join(cdir, 'postProcessing', 'sets', 'streamLines3D', '*', 'track0.vtp')))[-1]
        tr = XMLPolyDataReader(FileName=[trk])
        sc = Transform(Input=tr); sc.Transform.Scale = [1000, 1000, 1000]            # m -> mm, same units as the STL
        sc.TransformAllInputVectors = 0                                              # keep velocities in m/s
        top = [fr['c'][i] + 12 * fr['n'][i] for i in range(3)]                       # crop the aortic extension 12 mm above the valve
        sc = Clip(Input=sc); sc.ClipType = 'Plane'; sc.ClipType.Origin = top; sc.ClipType.Normal = fr['n']; sc.Invert = 1
        calc = Calculator(Input=sc); calc.ResultArrayName = 'speed_rel'; calc.Function = f'mag({field})/{vmean}'
        wall = STLReader(FileNames=[os.path.join(CFD, 'geometry', str(cid), stl)])
        wall = Clip(Input=wall); wall.ClipType = 'Plane'; wall.ClipType.Origin = top; wall.ClipType.Normal = fr['n']; wall.Invert = 1
        if 'lv' not in tag:   # tube: show the same 12 mm above to 30 mm below the valve
            bot = [fr['c'][i] - 30 * fr['n'][i] for i in range(3)]
            for k, src in (('w', wall), ('t', sc)):
                cl = Clip(Input=src); cl.ClipType = 'Plane'; cl.ClipType.Origin = bot; cl.ClipType.Normal = fr['n']; cl.Invert = 0
                if k == 'w': wall = cl
                else: sc = cl
            calc = Calculator(Input=sc); calc.ResultArrayName = 'speed_rel'; calc.Function = f'mag({field})/{vmean}'
        v = CreateView('RenderView'); v.ViewSize = [900, 1100]; v.Background = [1, 1, 1]
        v.OrientationAxesVisibility = 0
        try:
            v.UseColorPaletteForBackground = 0
        except Exception:
            pass
        dw = Show(wall, v); dw.Opacity = 0.12; dw.DiffuseColor = [0.55, 0.55, 0.6]; dw.AmbientColor = [0.55, 0.55, 0.6]
        dt = Show(calc, v); ColorBy(dt, ('POINTS', 'speed_rel')); dt.LineWidth = 1.6
        lut = GetColorTransferFunction('speed_rel')
        for name in ('Inferno', 'Inferno (matplotlib)', 'Magma (matplotlib)', 'Viridis (matplotlib)'):
            try:
                lut.ApplyPreset(name, True); break
            except RuntimeError:
                pass
        lut.RescaleTransferFunction(0.0, VMAX)
        dt.SetScalarBarVisibility(v, False)
        # A5C-like view: look along e2 (mediolateral), flow axis n pointing up the screen, septum on the left
        c = list(fr['c'])
        cam = v.GetActiveCamera()
        cam.SetFocalPoint(*c)
        # choose the camera side so that screen-right = +e1 (mitral / LA) and the septum is on the left
        e1, e2, n = fr['e1'], fr['e2'], fr['n']
        cx = [e2[1] * n[2] - e2[2] * n[1], e2[2] * n[0] - e2[0] * n[2], e2[0] * n[1] - e2[1] * n[0]]
        sgn = 1 if sum(cx[i] * e1[i] for i in range(3)) > 0 else -1          # right = cross(view dir, up)
        cam.SetPosition(*[c[i] - 300 * sgn * e2[i] for i in range(3)])
        cam.SetViewUp(*fr['n'])
        v.ResetCamera(False)
        cam.Zoom(1.15)
        v.CameraParallelProjection = 0
        SaveScreenshot(png(tag, case), v, ImageResolution=[900, 1100], TransparentBackground=0)
        print('wrote', png(tag, case))
        for src in list(GetSources().values()):
            Delete(src)
        Delete(v)


def montage():
    import numpy as np
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    import matplotlib.image as mi
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import Normalize
    fig, axs = plt.subplots(1, 4, figsize=(18, 6.4))
    for ax, (tag, case, _, _, title) in zip(axs, PANELS):
        im = mi.imread(png(tag, case))[..., :3]; m = (im < 0.97).any(2); r_, c_ = np.where(m)
        im = im[max(r_.min() - 20, 0):r_.max() + 20, max(c_.min() - 20, 0):c_.max() + 20]
        ax.imshow(im); ax.set_axis_off(); ax.set_title(title, fontsize=10, loc='left')
        ax.text(0.02, 0.5, 'septum', transform=ax.transAxes, rotation=90, va='center', fontsize=8, color='0.4')
        ax.text(0.98, 0.5, 'mitral / LA', transform=ax.transAxes, rotation=90, va='center', ha='right', fontsize=8, color='0.4')
    cax = fig.add_axes([0.35, 0.1, 0.3, 0.022])
    fig.colorbar(ScalarMappable(Normalize(0, VMAX), 'inferno'), cax=cax, orientation='horizontal').set_label('speed ÷ true mean LVOT velocity at 7.5 mm')
    fig.suptitle('3-D streamlines seeded on the inflow surface, viewed like an apical 5-chamber view (aorta at the top)', fontsize=12)
    fig.text(0.5, 0.015, 'Each line follows the steady (or time-averaged) velocity from where blood enters. Grey shell = blood-pool surface. '
             'Ejection only, from rest: no diastolic vortex.', ha='center', fontsize=8.5)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.88, bottom=0.19, wspace=0.02)
    out = os.path.join(CFD, 'qc', 'streamlines3d_388_33.png'); fig.savefig(out, dpi=130); plt.close(fig)
    print('wrote', os.path.relpath(out, CFD))


if __name__ == '__main__':
    montage() if '--montage' in sys.argv else render()
