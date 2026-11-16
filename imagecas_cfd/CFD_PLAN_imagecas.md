# CFD plan: LVOT flow in ImageCAS systolic cases

Version 2026-11-11. Do this work after the LVOT pipeline (`README_lvot_curtain.md`) is complete.
Fix all outputs in section 7 before you run the first case.

## 1. Purpose
Find the effect of LVOT shape (ellipticity, curtain length) on LVOT flow at peak systole.
Keep the flow rate constant, so that geometry is the only variable.

## 2. Case selection
1. Use only cases with `systolic_flag` = True (LV/myocardium volume ratio < 0.6). There are approximately 25.
2. Mark the curtain end manually in these cases (the rule cannot find it in systole).
3. Select approximately 10 cases with a similar LVOT area at 2 mm (±15% of the median).
4. Make sure that these cases have a wide range of ellipticity.
5. Add 2–3 diastolic cases with a similar area for comparison.

## 3. Geometry
1. Merge the LV label (3) and the aorta label (6). Do not include the LA.
2. Make a surface with marching cubes. Make it smooth with a volume-preserving (Taubin) filter.
3. Measure LVOT area and ellipticity on the surface at 2 mm and at the curtain end.
4. If they differ by more than 2% from the label values, decrease the smoothing.
5. Do not model the aortic cusps. The root is an open passage (the valve is open in systole).
6. Cut the aorta 10 mm above the sinotubular junction.
7. Add a straight outlet tube of 10 diameters. Blend it to a circle of equal area.
8. Close the mitral orifice with a wall.
9. Cut the LV with a plane perpendicular to the LV long axis, at mid-cavity. This plane is the inlet (primary model).

## 4. Boundary conditions

### 4.1 Primary model (single inlet)
1. Inlet: the mid-cavity LV plane. Use a flat velocity profile.
2. Make sure that the inlet is more than 20 mm below the valve plane. If it is not, move it toward the apex.
3. Walls: all other surfaces, no-slip. This includes the full LVOT.
4. Outlet: pressure = 0.

### 4.2 Sensitivity model (wall inflow)
1. Inflow: the LV wall more than 20 mm below the valve plane, uniform normal velocity. Use the full LV cavity (no inlet plane).
2. Walls: the LV wall from the valve plane to 20 mm, and the aorta, no-slip.
3. Outlet: pressure = 0.

### 4.3 Inlet position test (one case)
1. Run the primary model with the inlet plane at 15, 20, 25 and 30 mm below the valve plane.
2. Use the smallest distance where the outputs in section 7 change by less than 3%.

## 5. Flow and fluid
1. Flow rate: 400 mL/s (peak systole). Basis: stroke volume approximately 80 mL × peak-flow/stroke-volume ratio approximately 5 per second.
2. Sensitivity: 320 and 480 mL/s (±20%).
3. Blood: Newtonian, density 1060 kg/m³, viscosity 3.5 mPa·s.
4. Flow: steady. Turbulence model k-ω SST (LVOT Reynolds number approximately 6,000–7,000).

## 6. Solver and mesh
1. Solver: OpenFOAM simpleFoam (steady), same settings as the UK Biobank mesh study.
2. Mesh: snappyHexMesh. Base cell 0.5 mm. Refine to 0.2 mm from the valve plane to 5 mm below the curtain end.
3. Add 3–5 boundary layers on the LVOT walls.
4. Do a mesh convergence test on one case with three meshes. Calculate the grid convergence index.
5. If wall shear stress does not converge, report only pressure and velocity.
6. Stop the solver when residuals are below 1e-5 and the outputs change by less than 0.5% over 200 iterations.

## 7. Outputs (fixed before the first run)
1. Pressure drop: from a plane 20 mm below the valve to a plane 2 mm above the valve.
2. Loss coefficient: pressure drop / (½ ρ V²), V = mean velocity at 2 mm.
3. At 2 mm, mid-LVOT and the curtain end: peak velocity and mean velocity.
4. Normalized flow displacement (NFD) at 2 mm, mid-LVOT, the curtain end and 5 mm above the valve plane:
   1. Calculate the "center of velocity": the centroid of the cross-section, weighted by the forward through-plane velocity. Use forward flow only.
   2. Calculate the distance from the area centroid to the center of velocity.
   3. Divide this distance by the equivalent diameter of the cross-section. This value is the NFD.
   4. Record the direction of the displacement. Divide it into two components: along the major axis and along the minor axis of the LVOT ellipse.
   5. Also record the component toward the septum (opposite to the LA direction).
   6. Also give each component divided by the axis length in that direction (major or minor).
5. Wall shear stress: mean and peak on the septal side and on the curtain side, separately (primary model only).
6. Effective flow area: flow rate / peak velocity, and its level.

## 7a. Analysis of NFD and ellipticity
1. Across cases: find the relation between NFD and ellipticity at each level. Include LVOT area as a covariate.
2. In the ellipticity test (section 8): plot NFD against ellipticity (1.1, 1.3, 1.5) for each case. This gives the effect of ellipticity alone.
3. Find the direction of the displacement: along the major axis, along the minor axis, or toward the septum.
4. NFD is the same parameter that 4D flow MRI uses for flow eccentricity in the aorta (Sigovan et al., JMRI 2011;34:1226–1230). Thus you can compare the CFD results with 4D flow data later.

## 8. Ellipticity test at constant area (2–3 cases)
1. In each short-axis slice from the valve plane to the curtain end, stretch the LVOT along its major axis and compress it along its minor axis by the same factor. The area does not change.
2. Blend the change to zero over 3 mm above the valve plane and 3 mm below the curtain end.
3. Make three shapes per case: ellipticity 1.1, 1.3 and 1.5.
4. Use the same flow, mesh rules and boundary conditions for all shapes.

## 9. Plausibility check
1. Compare the simulated LVOT peak velocity with the normal echo range (approximately 0.7–1.1 m/s).
2. If a normal-size LVOT gives a velocity far outside this range, stop. Examine the flow rate and the geometry.

## 10. Where to run
1. Run on the Mac with OpenFOAM in Docker, controlled by Claude Code.
2. Approximately 30–40 runs, 1–3 million cells each, 1–3 hours each on 4 cores.

## 11. Basis in the literature
- Static systolic CT geometry with steady peak systolic flow, validated against Doppler peak velocity (r = 0.97): Kohli et al., Ann Biomed Eng 2021;49:1449–1461.
- Flow rate from LV volume change over all CT phases: Kohli et al., Front Cardiovasc Med 2022;9:869259 (not possible in ImageCAS: one phase per patient).
- End-systolic gated CT geometry, pressure drop from LV apex to neo-LVOT: Pasta et al. (TMVR CFD).
- Steady mean systolic flow with one LV inlet: Kirchner et al., Front Cardiovasc Med 2025;10.3389/fcvm.2025.1665934.

## 12. Limits (write these in the paper)
1. Flow rate is from normal reference data, not from each patient.
2. No wall motion. No aortic valve leaflets.
3. Steady flow at one instant (peak systole).
4. The LVOT end in systolic cases is from manual marks.
