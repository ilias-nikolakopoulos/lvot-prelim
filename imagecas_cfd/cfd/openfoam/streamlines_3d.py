"""3-D streamlines of finished CFD runs (OpenFOAM streamLine function object, seeded on the inlet patch).

Step 1 (writes system/streamLines3D in each case listed in CASES):
  uv run --python 3.12 python -I cfd/openfoam/streamlines_3d.py
Step 2 (in the container, per case; the field is U, or UMean for the transient run):
  docker run --rm --user "$(id -u):$(id -g)" -v "$PWD/<case>:/case" -w /case opencfd/openfoam-run:2512 \
      -c "cd /case && postProcess -dict system/streamLines3D -latestTime -fields '(U)'"
  -> <case>/postProcessing/streamLines3D/<time>/*.vtk
Step 3: render with ParaView (render_streamlines.py, run by pvbatch).
"""
import os

CFD = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# run tag, case, velocity field, number of seeds
CASES = [('lam_h1', '388_real', 'U', 250),
         ('lv_lam_h1', '388_real', 'U', 400),
         ('lv_apical_lam_tr_h1', '388_real', 'UMean', 400),
         ('lv_lam_h1', '33_real', 'U', 400)]

DICT = """FoamFile {{ version 2.0; format ascii; class dictionary; object streamLines3D; }}
functions
{{
streamLines3D
{{
    type            streamLine;
    libs            (fieldFunctionObjects);
    setFormat       vtk;
    U               {field};
    fields          ({field});
    direction       forward;
    lifeTime        20000;
    nSubCycle       5;
    cloud           particleTracks;
    seedSampleSet
    {{
        type        patchSeed;
        patches     (lumen_inlet);
        axis        xyz;
        maxPoints   {n};
    }}
}}
}}
"""

if __name__ == '__main__':
    for tag, case, field, n in CASES:
        c = os.path.join(CFD, 'runs', tag, case)
        open(os.path.join(c, 'system', 'streamLines3D'), 'w').write(DICT.format(field=field, n=n))
        print('wrote', os.path.relpath(c, CFD), field, n)
