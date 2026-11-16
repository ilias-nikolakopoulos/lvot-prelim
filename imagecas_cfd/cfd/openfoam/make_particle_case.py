"""Pathlines: continue a finished transient run with massless-like tracer particles (OpenFOAM icoUncoupledKinematicCloud).

Copies runs/<src_tag>/<case> (mesh, the latest time folder, system) to runs/<src_tag>_particles/<case>, then:
  - controlDict: restart from the latest time, run T_RUN more seconds, keep all writes, add the particle cloud
  - constant/kinematicCloudProperties: N tracers (20 um, blood density -> Stokes time ~6 us << time step), seeded at
    random positions inside the LV cavity and LVOT (2-45 mm below the valve) at the restart time; walls rebound,
    outlet lets them escape; a vtkCloud function object saves every particle position every POS_EVERY steps
  - constant/kinematicCloudPositions: the seed points (picked inside the twin STL)
The original run is not changed.

Run from the imagecas folder:
  uv run --python 3.12 --with numpy --with scipy --with trimesh --with vtk python -I cfd/openfoam/make_particle_case.py --t-run 0.15
  caffeinate -i cfd/openfoam/run_all.sh lv_apical_lam_tr_h1_particles > cfd/runs/lv_apical_lam_tr_h1_particles.log 2>&1
Then render with render_pathlines.py (reads postProcessing/tracerPositions/*.vtp).
"""
import argparse, glob, json, os, shutil
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CFD = os.path.dirname(HERE)
RHO = 1060.0
POS_EVERY = 12       # time steps between saved particle positions (dt ~ 0.085 ms -> ~1 ms)

CLOUD_FO = """
    // tracer particles (pathlines), added by make_particle_case.py
    tracers
    {
        type            icoUncoupledKinematicCloud;
        libs            (lagrangianFunctionObjects);
        cloud           kinematicCloud;
        U               U;
    }
    // particle positions every POS_EVERY time steps (~1 ms), joined into pathlines by origId in render_pathlines.py
    tracerPositions
    {
        type            vtkCloud;
        libs            (lagrangianFunctionObjects);
        clouds          (kinematicCloud);
        fields          (U d);
        writeControl    timeStep;
        writeInterval   POS_EVERY;
    }
"""

PROPS = """FoamFile {{ version 2.0; format ascii; class dictionary; object kinematicCloudProperties; }}

solution
{{
    active          true;
    coupled         false;
    transient       yes;
    cellValueSourceCorrection off;
    maxCo           0.3;
    interpolationSchemes
    {{
        rho         cell;
        U           cellPoint;
        mu          cell;
    }}
    integrationSchemes
    {{
        U           analytical;
    }}
    sourceTerms
    {{
        schemes {{ U explicit 1; }}
    }}
}}

constantProperties
{{
    rho0            {rho};
    youngsModulus   1e9;
    poissonsRatio   0.35;
}}

subModels
{{
    particleForces
    {{
        sphereDrag;
    }}

    injectionModels
    {{
        seeds
        {{
            type            manualInjection;
            massTotal       0;
            parcelBasisType fixed;
            nParticle       1;
            SOI             {t0};
            positionsFile   "kinematicCloudPositions";
            U0              (0 0 0);
            sizeDistribution
            {{
                type        fixedValue;
                fixedValueDistribution {{ value 2e-5; }}
            }}
        }}
    }}

    dispersionModel         none;
    patchInteractionModel   localInteraction;
    localInteractionCoeffs
    {{
        patches
        (
            lumen_wall   {{ type rebound; e 1; mu 0; }}
            lumen_inlet  {{ type rebound; e 1; mu 0; }}
            lumen_outlet {{ type escape; }}
        );
    }}
    heatTransferModel       none;
    surfaceFilmModel        none;
    stochasticCollisionModel none;
    collisionModel          none;

    cloudFunctions
    {{
    }}
}}
"""


def latest_time(case):
    ts = [p for p in glob.glob(os.path.join(case, '*')) if os.path.basename(p).replace('.', '', 1).isdigit() and os.path.basename(p) != '0']
    return max(ts, key=lambda p: float(os.path.basename(p)))


def inside(m, p):
    """Points inside a closed surface, with VTK's compiled vtkSelectEnclosedPoints (trimesh's ray test is too slow here)."""
    import vtk
    from vtk.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray
    surf = vtk.vtkPolyData(); pts = vtk.vtkPoints(); pts.SetData(numpy_to_vtk(np.ascontiguousarray(m.vertices), deep=True)); surf.SetPoints(pts)
    f = np.hstack([np.full((len(m.faces), 1), 3), m.faces]).astype(np.int64).ravel()
    cells = vtk.vtkCellArray(); cells.SetCells(len(m.faces), numpy_to_vtkIdTypeArray(f, deep=True)); surf.SetPolys(cells)
    q = vtk.vtkPolyData(); qp = vtk.vtkPoints(); qp.SetData(numpy_to_vtk(np.ascontiguousarray(p), deep=True)); q.SetPoints(qp)
    sel = vtk.vtkSelectEnclosedPoints(); sel.SetInputData(q); sel.SetSurfaceData(surf); sel.SetTolerance(1e-6); sel.Update()
    return np.array([sel.IsInside(i) for i in range(len(p))], bool)


def seed(dst, cj, a):
    """Random seed points inside the twin surface, d_min-d_max mm below the valve, >= 0.6 mm from the wall."""
    import trimesh
    from scipy.spatial import cKDTree
    info = json.load(open(os.path.join(CFD, 'geometry', str(cj['id']), 'slices.json')))
    fr = {k: np.array(v) for k, v in info['frame'].items()}
    m = trimesh.load(os.path.join(CFD, 'geometry', str(cj['id']), a.stl), force='mesh')
    surf = cKDTree(trimesh.sample.sample_surface_even(m, 400000, seed=1)[0])
    rng = np.random.default_rng(1)
    lo, hi = m.bounds; pts = []
    while len(pts) < a.n:
        p = rng.uniform(lo, hi, size=(20000, 3))
        d = -(p - fr['c']) @ fr['n']
        p = p[(d > a.d_min) & (d < a.d_max)]
        p = p[surf.query(p)[0] > 0.6]
        p = p[inside(m, p)]
        pts.extend(p)
    pts = np.array(pts[:a.n]) * 1e-3
    with open(os.path.join(dst, 'constant', 'kinematicCloudPositions'), 'w') as fh:
        fh.write('FoamFile { version 2.0; format ascii; class vectorField; object kinematicCloudPositions; }\n\n')
        fh.write(f'{len(pts)}\n(\n' + ''.join('(%.9g %.9g %.9g)\n' % tuple(q) for q in pts) + ')\n')
    print(f'seeded {len(pts)} tracers')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--src-tag', default='lv_apical_lam_tr_h1'); ap.add_argument('--case', default='388_real')
    ap.add_argument('--t-run', type=float, default=0.15, help='seconds of physical time to track')
    ap.add_argument('--n', type=int, default=1500, help='number of tracer particles')
    ap.add_argument('--d-min', type=float, default=2.0); ap.add_argument('--d-max', type=float, default=45.0)
    ap.add_argument('--stl', default='real_lv_apical.stl')
    ap.add_argument('--reseed', action='store_true', help='only rewrite the seed positions of an existing particle case')
    a = ap.parse_args()
    src = os.path.join(CFD, 'runs', a.src_tag, a.case)
    dst = os.path.join(CFD, 'runs', a.src_tag + '_particles', a.case)
    if a.reseed:
        seed(dst, json.load(open(os.path.join(dst, 'case.json'))), a); return
    if os.path.exists(dst):
        raise SystemExit(f'{dst} exists; remove it first')
    cj = json.load(open(os.path.join(src, 'case.json')))
    os.makedirs(dst)
    tdir = latest_time(src); t0 = float(os.path.basename(tdir))
    shutil.copytree(os.path.join(src, 'constant'), os.path.join(dst, 'constant'))
    shutil.copytree(os.path.join(src, 'system'), os.path.join(dst, 'system'))
    os.makedirs(os.path.join(dst, os.path.basename(tdir)))
    for f in ('U', 'p', 'phi', 'U_0', 'phi_0', 'uniform'):           # flow state only; the field averages are not needed
        s = os.path.join(tdir, f)
        (shutil.copytree if os.path.isdir(s) else shutil.copy2)(s, os.path.join(dst, os.path.basename(tdir), f))
    shutil.copytree(os.path.join(src, '0.orig'), os.path.join(dst, '0'))   # Allrun copies 0 -> 0.orig; needs a 0 folder
    for f in ('Allrun', 'mesh.done'):
        shutil.copy2(os.path.join(src, f), dst)
    cj.update(particles=dict(t0=t0, t_run=a.t_run, n=a.n, seed_depth_mm=[a.d_min, a.d_max], source=f'{a.src_tag}/{a.case}'))
    json.dump(cj, open(os.path.join(dst, 'case.json'), 'w'), indent=1)

    # controlDict: restart, run t_run more, keep every write, drop field averaging, add the cloud
    cd = open(os.path.join(dst, 'system', 'controlDict')).read()
    cd = cd.replace('endTime         0.5;', f'endTime         {t0 + a.t_run:g};')
    cd = cd.replace('purgeWrite      2;', 'purgeWrite      0;').replace('writeInterval   0.05;', 'writeInterval   0.025;')
    i = cd.index('    // time-averaged velocity'); j = cd.index('    // PW-Doppler sample volumes')
    cd = cd[:i] + CLOUD_FO.replace('POS_EVERY', str(POS_EVERY)) + cd[j:]
    cd = cd.replace('fields          (U p UMean);', 'fields          (U p);')
    open(os.path.join(dst, 'system', 'controlDict'), 'w').write(cd)
    assert f'endTime         {t0 + a.t_run:g};' in cd and 'tracers' in cd

    # transport: the cloud needs the fluid density
    tp = os.path.join(dst, 'constant', 'transportProperties')
    open(tp, 'a').write(f'rhoInf          [1 -3 0 0 0 0 0] {RHO:g};\nrho             [1 -3 0 0 0 0 0] {RHO:g};\n')
    open(os.path.join(dst, 'constant', 'g'), 'w').write(
        'FoamFile { version 2.0; format ascii; class uniformDimensionedVectorField; object g; }\n'
        'dimensions [0 1 -2 0 0 0 0];\nvalue (0 0 0);\n')
    open(os.path.join(dst, 'constant', 'kinematicCloudProperties'), 'w').write(PROPS.format(rho=RHO, t0=t0))

    seed(dst, cj, a)
    print('wrote', os.path.relpath(dst, CFD), f'restart t={t0:g}, end t={t0 + a.t_run:g}')


if __name__ == '__main__':
    main()
