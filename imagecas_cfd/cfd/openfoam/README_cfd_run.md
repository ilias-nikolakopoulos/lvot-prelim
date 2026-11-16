# CFD runs: real vs circular twin (2026-11-13)

Solver: OpenFOAM v2512 (OpenCFD), Docker image `opencfd/openfoam-run:2512`, native arm64.
Setup followed the [OpenFOAM Docker wiki](https://gitlab.com/openfoam/core/openfoam/-/wikis/precompiled/docker).
The interactive `openfoam-docker` script is in `~/bin/`. Batch runs use `run_all.sh`, which calls `docker run` directly.
Docker Desktop must be running.

## Files
- `make_cases.py`: writes `cfd/runs/<tag>/<id>_<twin>/`, one OpenFOAM case per twin.
- `run_all.sh <tag> [glob]`: mesh and solve each case in turn. Cases with `solve.done` are skipped, so it is safe to restart. `MESH_ONLY=1` stops after meshing.
- `postprocess.py --tag <tag>`: writes `cfd/cfd_results_<tag>.csv` and prints the echo-error split.

## Setup
- **Flow:** steady simpleFoam, Newtonian ν = 3.5e-6 m²/s, Q = 400 mL/s (typical peak-systolic LVOT flow), the same in every twin. Re ≈ 6,000.
- **Boundary conditions:**
  - inlet (LV side): plug flow
  - outlet (aorta): p = 0
  - wall: no-slip
- **Models:**
  - `lam_h1`: laminar. This is the main run, as agreed.
  - `kOmegaSST_h1`: k-ω SST. A turbulence-model check on pair 1 (388/142).
- **Mesh** (snappyHexMesh):
  - background 1 mm, aligned with the LVOT frame
  - 0.5 mm at the wall and in the LVOT segment (d = −2 to 16 mm)
  - 3 wall layers
  - case 388_real: 528k cells, checkMesh OK, 100% layer coverage, 45 s to mesh
- **Echo sample volume:** U·n averaged over the slice centroid plus 6 points on a 1.5 mm ring, at 5, 7.5 and 10 mm (probe every 10 iterations).
- **Pressure drop:** centroid probes 1 Deq into the inlet extension and 2 Deq into the aorta.
- **Echo planes:** raw U and p on cutting planes at each echo depth, saved in each time folder, for figures.

## Output (`postprocess.py`)
For each case and depth:
- `area_ratio` = π·AP²/4 ÷ true area. It is 1 for the circle.
- `profile_ratio` = sample-volume velocity ÷ true mean velocity.
- `SV_ratio` = area_ratio × profile_ratio (echo SV ÷ true SV).
- `dp_Pa`
- `drift_pct`: change in the sample-volume velocity over the last 200 iterations. Above about 1% means not converged.

## Timing
About 0.65 s per iteration on 4 cores (388_real). The cap is 2000 iterations, so the worst case is about 22 min per case. It stops earlier if residuals reach p < 1e-5 and U < 1e-6.
Worst case is about 9 h for the 24 laminar cases and 1.5 h for the 4 SST cases.

## Commands (from the `imagecas` folder)
```bash
caffeinate -i cfd/openfoam/run_all.sh lam_h1 > cfd/runs/lam_h1.log 2>&1
caffeinate -i cfd/openfoam/run_all.sh kOmegaSST_h1 > cfd/runs/kOmegaSST_h1.log 2>&1
uv run --python 3.12 --with numpy --with pandas python -I cfd/openfoam/postprocess.py --tag lam_h1
```

## Still to do
- **Convergence:** steady laminar flow at Re ≈ 6,000 may not settle. Check `drift_pct` and the residuals in `log.simpleFoam`. If it oscillates, compare with SST or run a short transient (pimpleFoam).
- **Mesh check:** rerun one pair at a finer mesh (`--h 0.7 --tag lam_h0.7`) and compare SV_ratio.
- **Licence:** ImageCAS is CC BY-NC 4.0, and the `heartchambers_highres` labels are for non-commercial use only. Cite ImageCAS (Zeng 2023), the DTU labels (Hansen 2025) and TotalSegmentator (Wasserthal 2023) plus nnU-Net.

## Checks on pair 1 (388/142), done 2026-11-14
- **Mesh:** `lam_h0.7` (about 1.25M cells, 2.4× the cells) vs `lam_h1` (about 0.53M). The profile ratio changes by ≤ 0.2%, and SV_ratio by ≤ 0.002.
- **Turbulence model:** `kOmegaSST_h1` vs laminar. The profile ratio changes by ≤ 0.1%. SST was active (peak eddy viscosity ≈ 12× ν), but in a segment this short its effect on the LVOT velocity is negligible.
- **Shape effect (real/circle profile ratio)** is the same in all three runs: 1.035–1.038 for 142 and 1.007–1.008 for 388.
- So the laminar 1 mm results are converged with respect to mesh and turbulence model. The remaining uncertainty is the inflow (straight tube with plug inlet, no LV cavity), which these checks don't test.

## LV-cavity inflow variant (`lv_lam_h1`, pair 1), set up 2026-11-14
- **Geometry** (`cfd/build_lv_twins.py` → `geometry/{id}/{twin}_lv.stl`, `lv_info.json`, `qc_lv.png`):
  - the patient's LV cavity (DTU label 3), shared by both twins
  - a 4 mm signed-distance blend (12–16 mm) into the twin's LVOT segment (2–12 mm)
  - the same aortic outlet extension as before
  - cross-section areas at 5/7.5/10 mm match the tube twins within 0.7%
- **Inflow:** uniform inflow through the endocardium deeper than 20 mm and more than 3 mm from the LA label (patch `inlet`), Q = 400 mL/s. This stands in for wall contraction in a rigid model. The rest of the LV, the LVOT and the aorta are no-slip walls.
- **Mesh:** the inlet surface is at 1 mm, the walls and LVOT at 0.5 mm, with 3 layers on the walls. 142_real: 538k cells, checkMesh OK, ~0.4 s/iteration.
- **Pressure probe:** deep in the LV cavity instead of in the inlet tube.
- **Run:**
  ```bash
  caffeinate -i cfd/openfoam/run_all.sh lv_lam_h1 > cfd/runs/lv_lam_h1.log 2>&1
  ```
- **Post-process:**
  ```bash
  uv run --python 3.12 --with numpy --with pandas python -I cfd/openfoam/postprocess.py --tag lv_lam_h1
  uv run --python 3.12 --with numpy --with scipy --with pandas --with matplotlib python -I cfd/openfoam/plane_figures.py --tag lv_lam_h1 --one 142
  ```

## Methodological principles from the CFD literature and how this pilot handles them (2026-11-14)
Context: those studies mostly look at the jet *downstream* of a (stenotic) aortic valve. This pilot looks *upstream*, in the LVOT at the echo PW-Doppler site, in ImageCAS coronary-CT patients (not an AS cohort).

| Principle | What we do | Status |
|---|---|---|
| Inflow boundary condition drives asymmetry; a flat/parabolic inlet underestimates it | Main runs (`lam_h1`) use a straight tube with plug inflow, which removes LV-driven skew by design. `lv_lam_h1` adds the patient's own LV cavity with wall-distributed inflow. `lv_apical_lam_h1` restricts inflow to the deepest half of the LV, to test sensitivity to how inflow is distributed. | Pair 1 only; to run |
| Valve/orifice modelling | No leaflets (open valve). We sample 5–10 mm *below* the valve, where leaflets matter less than for the downstream jet. The outlet is a straight extension, with no sinuses or aorta. | Limitation; not tested |
| Turbulence: RANS is poor for post-stenotic jets; LES needed | Our region is upstream, accelerating flow (favourable pressure gradient). k-ω SST changed the result by ≤ 0.1% at peak. Deceleration and transition are not modelled (steady peak). | Checked at peak; transient not tested |
| Rigid walls, Newtonian blood | Both assumed. Newtonian is reasonable at LVOT shear rates. LVOT shape change through systole is only approximated by the end-systolic vs diastolic comparison. | Limitation |
| Mesh dependence | 2.4× finer mesh: ≤ 0.2% change. | Checked (pair 1) |
| Standard asymmetry metrics | Normalized flow displacement (velocity-weighted centre offset ÷ D_eq) and flow angle, as used in 4D-flow MRI (Sigovan et al. 2011), plus the in-plane/secondary velocity fraction. In `plane_figures.py` → `peak_displacement_<tag>.csv`. | Done |
| Validation against in vivo flow (4D flow MRI / PC-MRI) | None available for ImageCAS. Results are model comparisons (real vs same-area circle under identical assumptions), not validated absolute velocities. | Limitation; state in the talk |
| Small cohort | 6 pairs (pilot); LV-inflow variants in 1 pair. | Pilot only |

Straight-tube baseline (`lam_h1`): flow displacement 0.01–0.02, flow angle 3–8° (follows the centreline curve), secondary velocity 8–17% of mean, the same in real and circle twins.

## LV-inflow results, pair 1 (run 2026-11-14 10:20–13:15)
Tables: `cfd_results_{lv_lam_h1,lv_apical_lam_h1}.csv` and `peak_displacement_{...}.csv`. Figures: `qc/velocity_sections_{lv_lam_h1,lv_apical_lam_h1}_{142,388}.png`.

**Convergence:**
- `lv_lam_h1`: all 4 cases steady. The sample-volume velocity changes by 0.0% over the last 500 iterations. 142_real and 388_* reached 2000 iterations with residuals of 1e-5 to 5e-3, but the probes are flat.
- `lv_apical_lam_h1`: the 142 cases converged. **The 388 cases did not converge:** the sample-volume velocity swings 11–18% over the last 500 iterations. The apical jet is probably unsteady and needs a transient (pimpleFoam) run. The 388 apical values below are a single snapshot and should not be quoted.

SV_ratio (real twin) and shape effect (real ÷ circle profile ratio):

| case, depth | tube SV | LV SV | apical SV | tube shape | LV shape | apical shape |
|---|---|---|---|---|---|---|
| 142, 5 mm | 0.893 | 0.914 | 0.910 | 1.035 | 1.002 | 1.001 |
| 142, 7.5 mm | 0.819 | 0.823 | 0.819 | 1.037 | 0.998 | 0.998 |
| 142, 10 mm | 0.740 | 0.720 | 0.717 | 1.038 | 0.994 | 0.992 |
| 388, 5 mm | 0.901 | 0.954 | (0.867)* | 1.008 | 0.982 | (0.934)* |
| 388, 7.5 mm | 0.856 | 0.896 | (0.810)* | 1.007 | 0.982 | (0.960)* |
| 388, 10 mm | 0.870 | 0.898 | (0.818)* | 1.008 | 0.982 | (0.940)* |

\* not converged.

**Findings:**
- **The area error is unchanged** (area_ratio 0.64–0.81), so echo still underestimates SV by 5–28%.
- **LV inflow increases the profile ratio** in both twins: about 1.13 (142) and 1.17–1.21 (388), vs about 1.07–1.15 with the tube. On the LA side under the aorto-mitral curtain the flow separates and is near zero, so the core speeds up.
- **The shape effect (real vs circle) disappears or slightly reverses** (0.98–1.00). With realistic inflow, the velocity pattern is set by the LV, not by the LVOT cross-section shape: the real and circle twins show the same skew.
- **Asymmetry:** in 388, flow angle 16–21° (vs 9–10° in the tube), secondary flow 40–49% of mean (vs 18–20%), peak 1.26–1.34× mean, fast side away from the LA. 142 stays almost flat (angle 3–4°). Flow displacement is still small (≤0.05) because the fast region is broad.
- **Inflow distribution:** whole-LV vs apical inflow changes 142 by ≤0.5%. 388 is not steady with apical inflow, so it cannot be compared.

## Transient check of the unsteady apical case (`lv_apical_lam_tr_h1`, 388), set up 2026-11-14
- **Solver:** pimpleFoam, laminar, constant Q = 400 mL/s, backward Euler in time, max Courant 0.8 (Δt ≈ 0.085 ms), 0.5 s of physical time.
- **Start:** reuses the `lv_apical_lam_h1` mesh and final steady U and p (`make_cases.py --transient --init-from lv_apical_lam_h1`).
- **Averaging:** `fieldAverage` (UMean, UPrime2Mean) from t = 0.2 s, giving a 0.3 s window (about one systole). Probes are written every time step.
- **Time:** about 0.73 s per step, so about 70–75 min per case and about 2.5 h for both.
- **Post-processing:** `postprocess.py` reports the time-mean v_sv, plus v_sv_std_pct and v_sv_p2p_pct (fluctuation, % of mean velocity). `plane_figures.py --field UMean` draws the time-averaged sections (output named `<tag>_mean`).
- **Caveat:** constant peak flow for 0.5 s is not physiological. This answers whether the apical jet is unsteady at peak flow and what the echo sample volume sees on average. A pulsatile systolic waveform would be the next step.
```bash
caffeinate -i cfd/openfoam/run_all.sh lv_apical_lam_tr_h1 > cfd/runs/lv_apical_lam_tr_h1.log 2>&1
uv run --python 3.12 --with numpy --with pandas python -I cfd/openfoam/postprocess.py --tag lv_apical_lam_tr_h1
uv run --python 3.12 --with numpy --with scipy --with pandas --with matplotlib python -I cfd/openfoam/plane_figures.py --tag lv_apical_lam_tr_h1 --one 388 --field UMean
```

## Echo-style view (removed 2026-11-15)
The pseudo-echo / PW-gate placement figure (`echo_view.py`, `sample_echo_view.py`) is no longer part of the project. Its code and figures were moved to the Trash.

### Transient results (`lv_apical_lam_tr_h1`, run 2026-11-14 15:02–18:06)
Files: `cfd_results_lv_apical_lam_tr_h1.csv`, `peak_displacement_lv_apical_lam_tr_h1_mean.csv`, `qc/velocity_sections_lv_apical_lam_tr_h1_mean_*.png` and `qc/transient_trace_lv_apical_lam_tr_h1.png` (`transient_trace.py`).

**The flow is genuinely unsteady at constant peak flow.** The gate velocity fluctuates with an SD of 2.4–3.1% of the mean and a peak-to-peak range of 10–19%, irregularly, with no single period. This is shear-layer instability of the apical jet at Re ≈ 6,000, not a solver problem.

Time means over 0.2–0.5 s:

| 388, depth | SV real | SV circle | profile real | profile circle | shape (real ÷ circle) |
|---|---|---|---|---|---|
| 5 mm | 0.879 | 1.140 | 1.090 | 1.140 | 0.956 |
| 7.5 mm | 0.814 | 1.111 | 1.064 | 1.111 | 0.958 |
| 10 mm | 0.818 | 1.132 | 1.085 | 1.132 | 0.958 |

- **The steady snapshot was close:** steady-solver SV 0.867/0.810/0.818 vs time mean 0.879/0.814/0.818 (within 1.5%). The shape effect of about 0.94–0.96 also holds.
- **Time-averaged section:** a strong jet on the septal side (1.55–1.6× the mean). The gate sits on the steep slope down from it, and the flow centre is 0.9–1.5 mm off (flow displacement 0.03–0.06, angle 14–20°). This matches the steady result.
- **Echo SV ÷ true SV** stays 0.81–0.88: the area error still dominates.
- **Beat-to-beat noise:** a single PW reading taken at a random instant could vary by about ±3% (1 SD) from flow unsteadiness alone.

## LV inflow extended to all 12 IDs (`lv_lam_h1`), set up 2026-11-15
- `build_lv_twins.py --ids 146,167,253,33,467,576,586,587,663,978`: all 20 twins built, with area checks at 5/7.5/10 mm within 0.7%.
- `make_cases.py --variant lv` for the same IDs. The 142 and 388 cases are untouched (solve.done), so run_all.sh skips them.
- Run: `caffeinate -i cfd/openfoam/run_all.sh lv_lam_h1 > cfd/runs/lv_lam_h1_all.log 2>&1` (about 5–8 h; 663 is the largest mesh).
- Steady only. Any case whose gate probe still swings over the last 500 iterations gets a transient rerun (see the transient section).

### LV-inflow results, all 6 pairs (`lv_lam_h1`, run 2026-11-15 17:34–20:56)
Files: `cfd_results_lv_lam_h1.csv`, `peak_displacement_lv_lam_h1.csv`, `qc/velocity_sections_lv_lam_h1_{d5p0,d7p5,d10p0,33,142,388,576}.png`.

**Convergence:** all 24 cases converged.
- 12 cases met the residual targets in 260–830 iterations. Their probes are flat over the last 100 iterations.
- The other 12 reached the 2,000-iteration cap with flat probes. The least settled is 33_real at 5 mm (±0.8% over the last 500 iterations).
- No case needs a transient rerun.

**Echo SV ÷ true SV (median), with whole-LV inflow vs straight tube:**

| Arm | Twin | 5 mm | 7.5 mm | 10 mm | Tube, real (5 / 7.5 / 10 mm) |
|---|---|---|---|---|---|
| Diastolic | real | 0.81 | 0.70 | 0.68 | 0.81 / 0.71 / 0.71 |
| Diastolic | circle | 1.10 | 1.09 | 1.06 | |
| Systolic | real | 0.92 | 0.89 | 0.86 | 0.84 / 0.83 / 0.82 |
| Systolic | circle | 1.23 | 1.23 | 1.21 | |

**Findings:**
- **The shape effect disappears in every case.** The real ÷ circle gate velocity is 0.98–1.04 (median 0.996–1.002), vs 1.015–1.034 with the tube. With realistic inflow, the LVOT cross-section shape does not change what the gate reads.
- **Diastolic scans: almost no change from the tube.** The gate reads 1.06–1.14× the mean. Echo underestimates SV by 19–32% because of the AP-diameter area.
- **Systolic scans: the gate over-reads more** (1.16–1.47× the mean in the circle twin at 5 mm). This partly offsets the area error.
  - In 2 of 6 pairs (33 and 576), echo SV is 0.99–1.23 at 5–7.5 mm, so echo **overestimates** SV there.
  - Cause: a large separated (near-zero) region on the septal side, just below the valve, in both twins. The core speeds up to 1.4–1.5× the mean.
- **The gate over-read follows LV size.** Spearman ρ = −0.74 between LV label volume (30–152 mL) and the gate ÷ mean velocity at 5 mm. The smallest cavities (33: 30 mL; 576: 46 mL) give the largest over-read.
- **Asymmetry (median):**
  - Diastolic: flow angle 4–5°, secondary flow 14–21% of the mean, flow displacement 0.01–0.02.
  - Systolic: flow angle 9–11°, secondary flow 35–41%, flow displacement 0.01–0.04.

**Caveat:** the systolic scans are end-systolic, so the LV is at its smallest. Peak ejection happens in a larger cavity. Pushing 400 mL/s through an end-systolic cavity probably exaggerates the flow convergence and the septal separation. So the systolic-arm profile values are an upper bound. The diastolic arm (a larger cavity) is probably closer to the cavity size at peak ejection, but it has the open-mitral-leaflet geometry.

## Long-axis flow pattern (2026-11-15)
- `long_axis_flow.py write|plot --id 388 --twin real` samples U (or UMean) on the A5C-like plane with `postProcess -dict system/sampleLongAxis` and draws in-plane streamlines for the tube, whole-LV and apical (transient-mean) runs → `qc/long_axis_flow_388_real.png`.
- The apical run shows a septal jet and a large eddy under the anterior mitral leaflet. The whole-LV run shows smooth convergence with a thin slow septal layer near the valve.
- No run includes the diastolic LV vortex: each is ejection only, starting from rest.

## 3-D streamlines and particle pathlines (2026-11-15)
- **Streamlines:** `streamlines_3d.py` writes `system/streamLines3D` (OpenFOAM `streamLine`, seeded on `lumen_inlet`). Rendered with ParaView 6.2 `pvbatch render_streamlines.py` → `qc/streamlines3d_388_33.png`, showing 388 tube / whole-LV / apical and 33 whole-LV.
- **Pathlines:** `make_particle_case.py` copies `lv_apical_lam_tr_h1/388_real` to `lv_apical_lam_tr_h1_particles/388_real` and adds an `icoUncoupledKinematicCloud` function object.
  - 1,500 tracers (20 µm, ρ = 1060), seeded inside the LV 2–45 mm below the valve at t = 0.5 s, run for 0.15 s (run 23:23–23:44).
  - Positions are written every 12 steps by `vtkCloud`. `render_pathlines.py` links them (predicted position + Hungarian assignment; median error 0.001 mm, 1 rejected link out of 216,604) → `qc/pathlines_lv_apical_lam_tr_h1_particles_388_real.png`.
- **Pathline result:** ejected by 150 ms: 100% (2–10 mm), 89% (10–20), 81% (20–30), 77% (30–45). Median time to the valve: 8, 21, 50 and 70 ms. Deep parcels loop through a helical LV vortex before leaving along the septum.

## Emulated PW spectral read-out (2026-11-16, requested via the lvot-ukb session; Huh et al. 2022)
- `echo_spectrum.py write|analyse` and `echo_spectrum_sample.sh`. The gate is 4 mm long (uniform) with a Gaussian beam (FWHM 3 mm), at the centroid at 5 / 7.5 / 10 mm, sampled every 0.25 mm with OpenFOAM `sets` (`cloud`) for `lam_h1` and `lv_lam_h1` (48 cases, 4 s each).
- The spectrum is the volume-weighted histogram of beam velocity (2 cm/s bins). Read-outs: mean, modal, and max (97.5th percentile). Beam angle 0, ±10, ±20°.
- **Result:** modal ≈ mean (+0–4%), max = mean +3–8%. Real-LVOT echo SV with the max read-out is 0.73–0.96 (median 0.85, LV inflow), so the read-out does not explain Huh's +29% VTI. Angle ±10° changes SV ≤ 4%, ±20° up to −10%.
- Outputs: `cfd/echo_spectrum_{lam_h1,lv_lam_h1}.csv`, `cfd/echo_spectrum_summary.csv`, `qc/echo_spectrum_readout.png`.
