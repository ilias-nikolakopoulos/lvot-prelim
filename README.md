# LVOT size and ellipticity in UK Biobank representative meshes

This is a personal side project, kept separate from `../ukb-mesh/`. It reuses only that project's `.venv` (labels are rebuilt locally) and never modifies it or the raw data.
> **All measurements are end-diastolic.** The released meshes were built from end-diastolic frames only (Ugurlu et al. 2025), so systolic LVOT size and shape cannot be measured here. See `docs/DECISIONS.md` #16.

Data: 1,423 representative biventricular meshes (Ugurlu et al., PLoS One 2025; Zenodo 10.5281/zenodo.15649643, CC0) at `data/meshes`. Each mesh is the average of the participants with the same sex, age (46–83) and integer BMI. The results therefore describe typical hearts, not individuals, and CIs and p-values are descriptive only (participant counts per mesh are unknown).

## Reproduce
- **Docs:** step-by-step pipeline with the code for each step in [`docs/PIPELINE.md`](docs/PIPELINE.md); flowchart in [`docs/flowchart.mmd`](docs/flowchart.mmd) (Mermaid); reasons for each choice in [`docs/DECISIONS.md`](docs/DECISIONS.md).
- **Parameters:** everything in [`config.yaml`](config.yaml).
- **Environment:** Python 3.11 + [`requirements.txt`](requirements.txt) (currently `../ukb-mesh/.venv`).
- **Provenance:** every run appends to `outputs/logs/runs.jsonl` (code/config hash, package versions).
- **Code:** each step is marked `# STEP x.y`. Versioned with git.

```
bash run_all.sh                                 # everything, in order (see the file for each step)
PY=../ukb-mesh/.venv/bin/python
$PY scripts/00_labels.py                        # STEP 0: data/labels.csv from raw folder names
$PY scripts/01_measure.py --subset              # STEP 1-5.1 on the 8-mesh subset -> outputs/tables/lvot_subset.csv
$PY scripts/01_measure.py --all                 # all meshes -> outputs/tables/lvot_all.csv (resumable, ~3 s/mesh)
$PY scripts/02_qc_render.py --subset            # STEP 5.2 QC renders -> outputs/figures/qc/<id>.png
$PY scripts/04_gap_demo.py <mesh_id>            # STEP 5.3 chord-closing illustration
$PY scripts/05_outline_demo.py <mesh_id>        # STEP 5.4 how an outline is obtained (3D + en face)
$PY scripts/03_analysis.py                      # STEP 6 figures + tables
$PY scripts/06_area_methods.py                  # STEP 6.6 outline area vs ellipse / echo circular area
$PY scripts/07_area_methods_demo.py             # STEP 6.7 illustration of the four area estimates
```
Full run (~70 min). Run it yourself in Terminal from this folder:
`caffeinate -i ../ukb-mesh/.venv/bin/python scripts/01_measure.py --all > outputs/logs/measure_all.log 2>&1`

## Methods (`src/lvot.py`)
These are geometry only: the `.fiber.vtu` file is read once and µm are converted to mm. UVCs are not used, so the LV/RV junction UVC artefact does not arise.
- **LV cavity.** The surface of all tets (myocardium plus the valve plugs) has 3 closed components: the outer surface, the RV cavity, and the LV cavity bounded by LV tets plus the mitral (tag 4) and aortic (tag 5) plugs. LV cavity volume is the enclosed volume of the LV component, including the LVOT up to the aortic plug. It equals LV end-diastolic volume, because the meshes are end-diastolic. **LV endocardium** = the cavity triangles that come from LV tets.
- **Rings.** For the aortic (tag 5) and mitral (tag 4) plugs, the ring is the set of LV nodes shared with the plug's tets: a band about 2 mm tall. A plane is fitted to it by SVD, with the normal pointing into the LV. The nodes are projected onto the plane, and the median radius is taken in each of 72 angular bins to give a polygon. From that polygon:
  - area → area-derived diameter; perimeter → perimeter-derived diameter;
  - **moment-equivalent ellipse** (same area and second moments as the polygon) → major/minor diameter, ellipticity (major/minor, and 1 − minor/major);
  - Feret (caliper) max/min diameters as a cross-check. On the subset they agree with the moment ellipse to within ~0.5 mm.

  A direct (Fitzgibbon) ellipse fit was tried first. It exaggerated ellipticity on D-shaped loops, so it was dropped.
- **LVOT sections.** These are planes perpendicular to the aortic-ring normal, at depths of 0–15 mm (1 mm steps) below the ventricular edge of the aortic ring band. Each plane is intersected exactly with the LV endocardium and the segments are chained into loops.
  - *Strict (primary):* the loop must be closed and must encircle the axis. Sections are kept only while they are contiguous from 0–1 mm.
  - **Finding:** the LV endocardium below the aortic ring is a fully closed tube for only **2–5 mm** (pilot of 60 meshes: 2 mm in 8, 3 mm in 41, 4 mm in 11). Below that, the section opens on the **mitral side**: the gap centre lies 0–15° from the direction of the mitral centre. The meshes have no mitral leaflets, so the anterior mitral leaflet, which forms the posterior LVOT wall, is missing. The aorto-mitral curtain is only ~3 mm long.
  - **Primary LVOT measure = the strict section 2 mm below the annulus** (closed in every pilot mesh).
  - *Chord-closed (secondary):* for open sections, the wall points within 0.9 × annulus diameter of the axis are ordered by angle and the gap is closed with a straight chord (an approximation of the missing leaflet). This is kept only if the gap is ≤ 120°. **Secondary measure = 4 mm** (gap 0–85°, always ≤ 120°). The 6 mm section was rejected: its gap exceeded 120° in 5 of 30 women and 0 of 30 men, which would make the missingness sex-biased.
  - Columns are `lvot{0,2,4,6,8,10}_*` (strict) and `lvot{d}cc_*` (chord-closed, plus `gap_deg` and `gap_vs_mitral_deg`).
- **Angles and lengths.**
  - Aorto-mitral angle = the angle between the ring-plane normals.
  - LV long axis = mitral-ring centre → farthest LV endo point (apex); LV length = that distance.
  - `ao_axis_vs_lv_long_angle` = the angle between the aortic-ring normal and the LV long axis, used as an aorto-septal angle proxy.
  - `lvot_centreline_vs_ao_normal` = the tilt of the line through the strict section centroids. It is noisy over only 2–5 mm.
- **QC flags.**
  - `ring_gap`: angular coverage below 95%.
  - `ao_size_implausible`: annulus outside 12–32 mm.
  - `mv_size_implausible`: mitral ring outside 20–55 mm.
  - `no_closed_lvot_section` / `lvot_closed_lt2mm`: the strict sections fail.
  - `surface_regions_k`: the surface does not have exactly 3 components.

  Per-depth reasons are given in `section_reasons`. Exclusions are logged by `03_analysis.py` to `outputs/tables/exclusions.csv`.

## Analysis (`scripts/03_analysis.py`)
- Colours: women = vermillion `#D55E00`, men = blue `#0072B2` (Okabe-Ito, CVD-validated).
- `fig_raw`: all meshes vs age (BMI not held fixed).
- `fig_common`: common support only (same age and BMI present in both sexes; 622 pairs), per-age means.
- `fig_adjusted`: OLS `y ~ age_c + sex + age_c×sex + BMI` (age centred at 65) on all meshes, predicted at BMI 27.
- `fig_sexdiff`: M − F by age, showing the common-support per-age paired mean, the model line ± 95% CI, and paired age-block means.
- Tables: `model_coefs.csv` (the sex difference at 50/65/80, also as %; age slopes by sex; the age×sex p-value; `narrows_with_age` = the interaction has the opposite sign to the sex effect) and `sexdiff_by_age_block.csv` (age blocks 46–54, 55–59, 60–64, 65–69, 70–74, 75–83).
- Indexed size: diameter / LV cavity volume^(1/3) (mm/mL^⅓) and diameter as % of LV length.

## Results
_See [`WRITEUP.md`](WRITEUP.md) (caveats final; results after the full run)._
