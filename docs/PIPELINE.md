# LVOT pipeline — step by step

> Phase: **end-diastole only**. The dataset contains no systolic meshes (`DECISIONS.md` #16).

Each step below lists its purpose, inputs, method, the code that does it, its outputs and the check that it worked. In the code, every step is marked with a comment `# STEP x.y` (search for it). Line links point to the current version and may drift slightly as the code changes. The flowchart source is [`flowchart.mmd`](flowchart.mmd), and the reasons behind each choice are in [`DECISIONS.md`](DECISIONS.md).

**Environment.** Python 3.11 with the exact package versions in [`requirements.txt`](../requirements.txt) (currently the `../ukb-mesh/.venv` interpreter). All parameters are in [`config.yaml`](../config.yaml). Every script run appends a line to `outputs/logs/runs.jsonl` with the time, arguments, code hash, config hash and package versions ([`src/provenance.py`](../src/provenance.py)).

**Reproduce everything:** `bash run_all.sh`. The full measurement step takes ~70 min, so run that one yourself in Terminal with `caffeinate`.

---

## STEP 0 — Label table
| | |
|---|---|
| Purpose | List every mesh with its sex, age and BMI |
| Input | Raw mesh folders under `data_dir` (read-only) |
| Method | Find every `<id>/<id>.fiber.vtu`, skipping hidden and `._` files. Parse `<sex>_age_<a>_bmi_<b>`. Check that ids are unique and that there is one mesh per sex × age × BMI. |
| Code | [`scripts/00_labels.py`](../scripts/00_labels.py), `parse_name()` in [`src/lvot.py:206`](../src/lvot.py:206) |
| Output | `data/labels.csv` (mesh_id, sex, age, bmi, folder) |
| Check | 1,423 meshes (F 761, M 662), ages 46–83, 622 common-support pairs. Identical to `../ukb-mesh/outputs/tables/labels.csv`. |

## STEP 1 — Load one mesh
Code: `measure_mesh()` in [`src/lvot.py:216`](../src/lvot.py:216), called once per mesh by [`scripts/01_measure.py`](../scripts/01_measure.py).

| Sub-step | Method | Code |
|---|---|---|
| 1.1 Read | `pv.read` the `.fiber.vtu` once; check all cells are tets; points µm × 1e-3 → mm. elemTag: 0 RV, 1 LV, 4 mitral plug, 5 aortic plug, 6 tricuspid, 7 pulmonary. | [`STEP 1.1`](../src/lvot.py:219) |
| 1.2 LV myocardium volume | Sum of tag-1 tet volumes (`tet_volumes()`) | [`STEP 1.2`](../src/lvot.py:233) |

## STEP 2 — LV cavity and endocardium (geometry, no UVCs)
| Sub-step | Method | Code |
|---|---|---|
| 2.1 Boundary surface | Boundary surface of **all** tets (myocardium + all 4 plugs), keeping the original cell id so each triangle carries its tet's tag. Split into connected components. | [`STEP 2.1`](../src/lvot.py:238) |
| 2.2 LV cavity | Expect 3 closed components: outer surface, RV cavity, LV cavity. The LV cavity is the component containing aortic and mitral plug faces and no RV faces. Its enclosed volume = LV cavity volume (`lv_cavity_ml`, LV EDV proxy). | [`STEP 2.2`](../src/lvot.py:245) |
| 2.3 LV endocardium | Cavity triangles that come from LV tets (plug faces removed) | [`STEP 2.3`](../src/lvot.py:259) |

Check: `n_regions == 3`, otherwise flag `surface_regions_k`.

## STEP 3 — Valve rings and axes
| Sub-step | Method | Code |
|---|---|---|
| 3.1 Ring nodes | LV nodes shared with the plug's tets: tag 5 = aortic, tag 4 = mitral. This is a ~2 mm band. | [`STEP 3.1-3.4`](../src/lvot.py:264) |
| 3.2 Plane | SVD plane fit to the ring nodes; normal oriented into the LV; nodes projected onto the plane; `band_max_depth` = how far the band reaches towards the ventricle | `ring_measure()` [`src/lvot.py:171`](../src/lvot.py:171) |
| 3.3 Outline | Centre from a direct ellipse fit, then the median radius in each of 72 angular bins (5° each) → polygon. Area → area-derived diameter; perimeter → perimeter-derived diameter. Angular coverage = fraction of bins filled. | `ring_measure()` |
| 3.4 Diameters / shape | **Moment-equivalent ellipse** of the polygon (same area and second moments) → major/minor diameter, ellipticity = major/minor and 1 − minor/major. Feret (caliper) max/min as a cross-check. | `moment_ellipse()` [`src/lvot.py:73`](../src/lvot.py:73) |
| 3.5 Aorto-mitral | Angle between the two ring-plane normals; distance between the ring centres | [`STEP 3.5`](../src/lvot.py:275) |
| 3.6 LV axes | Apex = LV endo point farthest from the mitral ring centre. LV length = that distance. LV long axis = mitral centre → apex. `ao_axis_vs_lv_long_angle` = angle between the aortic-ring normal and the long axis (aorto-septal angle proxy). | [`STEP 3.6`](../src/lvot.py:280) |
| 3.7 Echo-style annulus diameters | `d_plax` = edge-to-edge line through the outline's area centroid towards the mitral centre (the parasternal long-axis direction, septum → aorto-mitral side). `d_perp` = the orthogonal line. `plax_vs_minor_deg` = angle between the PLAX line and the minor axis of the moment ellipse. | [`STEP 3.7`](../src/lvot.py), `echo_diams()`, `chord_through()` |

## STEP 4 — LVOT cross-sections
| Sub-step | Method | Code |
|---|---|---|
| 4.1 Planes | Perpendicular to the aortic-ring normal, at depths 0–15 mm in 1 mm steps. Depth 0 = the ventricular edge of the aortic ring band, where the plug meets the LV endocardium. | [`STEP 4.1`](../src/lvot.py:287) |
| 4.2 Intersection | Exact intersection of the plane with the LV endocardium (`endo.slice`). The 2-point segments are chained into polylines, merging duplicate points; each polyline is closed (every vertex has 2 neighbours) or open. | [`STEP 4.2`](../src/lvot.py:298), `chain_loops()` [`src/lvot.py:129`](../src/lvot.py:129) |
| 4.3 Strict section | Keep a closed loop that encircles the LVOT axis (taking the outermost if there are several). Reject it if its area-derived diameter is > 1.6 × the annulus diameter, which means it has opened into the cavity. The axis point = centre of the previous closed section, starting at the ring centre. Measures: area (shoelace formula on the outline), area/perimeter-derived diameter, moment-ellipse major/minor, Feret, echo-style `d_plax` / `d_perp` (as in 3.7). | [`STEP 4.3`](../src/lvot.py:302) |
| 4.4 Chord-closed section (secondary) | Used when the section is open. Take wall points within 0.9 × annulus diameter of the axis. Find the largest empty angular sector (72 bins) = the gap, and record its size and its angle from the mitral direction. If the gap is ≤ 120°, order the points by angle, close with a straight chord, and measure as in 4.3 (including `cc_area`, Feret, `d_plax`, `d_perp`). | [`STEP 4.4`](../src/lvot.py:319) |
| 4.5 Contiguity | Strict sections count only as a contiguous run starting at 0 or 1 mm | [`STEP 4.5`](../src/lvot.py:354) |
| 4.6 Report depths | Write `lvot{0,2,4,6,8,10}_*` (strict) and `lvot{d}cc_*` (chord-closed, plus `gap_deg` and `gap_vs_mitral_deg`) | [`STEP 4.6`](../src/lvot.py:362) |
| 4.7 Tube length, tilt | `lvot_tube_len` = deepest strict section. `lvot_centreline_vs_ao_normal` = line through the strict section centroids vs the ring normal. | [`STEP 4.7`](../src/lvot.py:380) |

**Analysis choices:** primary LVOT = strict section at **2 mm**; secondary = chord-closed at **4 mm** (`config.yaml → analysis`).

## STEP 5 — Quality control
| Sub-step | Method | Code / output |
|---|---|---|
| 5.1 Flags | `ring_gap` (ring coverage < 0.95), `ao_size_implausible` (12–32 mm), `mv_size_implausible` (20–55 mm), `no_closed_lvot_section`, `lvot_closed_lt2mm`, `surface_regions_k`. A reason is recorded for every depth (`section_reasons`). | [`STEP 5.1`](../src/lvot.py:392) → columns `qc_flags` and `section_reasons` |
| 5.0 Run control | One row per mesh written as soon as it is done; meshes already in the CSV are skipped (resumable); errors go to `<out>_errors.csv`; progress goes to the log | [`scripts/01_measure.py`](../scripts/01_measure.py) → `outputs/tables/lvot_all.csv`, `outputs/logs/*.log` |
| 5.2 QC renders | Per mesh: en-face view (ring nodes, strict and chord-closed sections, ellipses, mitral direction) and side view (LV endo, aortic and mitral rings, depth planes) | [`scripts/02_qc_render.py`](../scripts/02_qc_render.py) → `outputs/figures/qc/<id>.png` |
| 5.3 Gap demo | Step-by-step picture of chord closure at 2/4/6/8 mm | [`scripts/04_gap_demo.py`](../scripts/04_gap_demo.py) → `outputs/figures/qc/gap_demo_<id>.png` |
| 5.4 Outline demo | 3D view of the LVOT wall cut by the 2 mm and 4 mm planes, plus en-face views: closed chain → polygon (2 mm); open chain + chord → polygon (4 mm); area → area-derived diameter; echo PLAX calliper line and its circular area | [`scripts/05_outline_demo.py`](../scripts/05_outline_demo.py) → `outputs/figures/qc/outline_demo_<id>.png` |
| 5.5 Outline, step by step | 9 panels for one mesh: (1) tagged mesh; (2) outer skin = 3 closed surfaces; (3) LV cavity from inside, plug faces dropped, aortic rim plane; (4) cutting planes; (5) plane × triangles → segments; (6) segments → closed chain at 2 mm + checks; (7) open chain at 4 mm; (8) chord closure; (9) polygon → shoelace area and diameters | [`scripts/08_outline_steps.py`](../scripts/08_outline_steps.py) → `outputs/figures/qc/outline_steps_<id>.png` |
| 5.6 Save outlines | Per mesh, saves the 2D outlines (annulus and mitral ring polygons, strict and chord-closed LVOT outlines at the report depths) and the PLAX directions, so any new diameter or area method can be computed without re-reading the meshes. Resumable. | [`scripts/09_save_outlines.py`](../scripts/09_save_outlines.py) → `data/outlines/<id>.npz` |

## STEP 6 — Age × sex analysis
Code: [`scripts/03_analysis.py`](../scripts/03_analysis.py).

| Sub-step | Method | Code / output |
|---|---|---|
| 6.1 Exclusions | Processing errors → whole mesh. Implausible aortic ring → aortic and LVOT measures set to NaN. Implausible mitral ring → mitral measures NaN. A failed strict 2 mm or chord-closed 4 mm section is logged with its reason. | `# STEP 6.1` → `outputs/tables/exclusions.csv` |
| 6.2 Indexing | diameter / LV cavity volume^(1/3); diameter as % of LV length | `# STEP 6.2` |
| 6.3 Common support | F/M pairs with identical age and BMI (622 pairs); within each pair BMI is fixed | `# STEP 6.3` |
| 6.4 Statistics | OLS `y ~ age_c + sex + age_c×sex + BMI` on all meshes (age centred at 65, BMI at 27). Reported: sex difference at 50/65/80 (mm and %), age slope by sex, age×sex p-value, `narrows_with_age`. Paired M − F mean ± 95% CI per age block (46–54 … 75–83). | `# STEP 6.4` → `model_coefs.csv`, `sexdiff_by_age_block.csv` |
| 6.5 Figures | `fig_raw` (all meshes), `fig_common` (common support), `fig_adjusted` (model at BMI 27), `fig_sexdiff` (M − F vs age). Women `#D55E00`, men `#0072B2`. | `# STEP 6.5` → `outputs/figures/analysis/` |
| 6.6 Area methods | For the annulus, LVOT 2 mm and LVOT 4 mm (chord-closed), the measured outline area is compared with: π/4·major·minor (moment ellipse); π/4·Feret max·min; π/4·longest·perpendicular, both through the centre (CT reading, `ct_diams()`, from the saved outlines of STEP 5.6); π/4·D_PLAX·D_perp (two orthogonal views); π/4·D_PLAX² (standard echo, circular LVOT). Reported: % difference by sex; paired M − F difference in echo bias. | [`scripts/06_area_methods.py`](../scripts/06_area_methods.py) → `area_methods_summary.csv`, `area_methods_by_mesh.csv`, `fig_area_methods.png` |
| 6.7 Area methods illustration | One outline (LVOT 2 mm) drawn once per method, showing the diameters used and the assumed shape (area and % error in the legend): moment ellipse, Feret ellipse, minor–major axis ellipse (longest × perpendicular through the centre, as in a CT reading), echocardiographic area (approximation of the echo PLAX view). Method names in bold `#102D40`; assumed shapes in `#B94F00`; font Roboto from `assets/fonts/`. The two-view (PLAX × perpendicular) ellipse is hidden (`SHOW_TWO_VIEW = False`) but kept in the code and in the tables. Below: mean ± SD % difference by method and sex for the aortic annulus, 2 mm and 4 mm from the annulus (bars coloured by level, semi-transparent; women lighter, men darker; title "Estimating the elliptical LVOT area"). | [`scripts/07_area_methods_demo.py`](../scripts/07_area_methods_demo.py) → `fig_area_methods_demo.png` |
| 6.8 Ellipticity by method | Ellipticity = major / minor for the three ellipse methods (moment axes; widest / narrowest calliper width; longest diameter through the centre / the one at 90°) at the annulus, 2 mm and 4 mm, by sex: median, IQR, 5th–95th percentile. Box plots in the same style as 6.7 (Roboto, level colours, women lighter, men darker). | [`scripts/10_ellipticity_methods.py`](../scripts/10_ellipticity_methods.py) → `ellipticity_methods.csv`, `fig_ellipticity_methods.png` |
| 6.9 PLAX plane definitions (deferred, not presented) | Mitral centre = centre of the ellipse fitted to the mitral rim nodes (projected on their plane). Three PLAX plane definitions: A aortic axis + mitral centre (used in 3.7); B aortic axis + LV apex; C aortic centre + mitral centre + apex. For each: angle from A and from the LVOT minor axis, D_PLAX and echo-circle error at 2 mm. Runs on the subset + 30 pilot pairs (68 meshes). | [`scripts/11_plax_definitions.py`](../scripts/11_plax_definitions.py) → `plax_definitions.csv`, `fig_plax_definitions.png` |

## STEP 7 — Write-up
`WRITEUP.md`: methods summary, key numbers, figures, caveats. The results are filled in after the full run. The caveats section already states that **all measurements are end-diastolic** (no systolic meshes exist in this dataset).


## STEP 7b — Highlights deck
`outputs/LVOT_highlights.pptx` (12 slides): key messages, data and method, the short closed tube, size by sex and age (native charts), indexing, ellipticity, echo underestimation, area methods, caveats, reproducibility and next steps. Built from the current tables and figures by [`scripts/12_highlights_deck.js`](../scripts/12_highlights_deck.js) (pptxgenjs; run with `NODE_PATH` pointing to a folder where pptxgenjs is installed; theme colours applied with the pptx skill's `apply_theme.js`). Font Arial. The deferred PLAX-plane and two-view analyses are not shown.
---
### Data flow
```
raw meshes (read-only) ─▶ 00_labels ─▶ data/labels.csv
                                          │
                     ┌────────────────────┴─────────────┐
                01_measure (STEPS 1-5.1, resumable)     02_qc_render / 04_gap_demo / 05_outline_demo / 08_outline_steps (STEP 5.2-5.5)
                     │                                    │
       outputs/tables/lvot_all.csv                 outputs/figures/qc/
                     │
                03_analysis + 06_area_methods (STEP 6) ─▶ tables (model_coefs, sexdiff_by_age_block, exclusions) + figures/analysis
                     │
                WRITEUP.md (STEP 7)
every script ─▶ outputs/logs/runs.jsonl (provenance)
```
