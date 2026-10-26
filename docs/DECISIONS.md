# Decisions and steps taken (chronological log)

One entry per choice: what was tried, what was seen, what was decided, and where the evidence is. Add new entries at the bottom.

### 2026-10-23
1. **Separate project.** `lvot-ukb/` sits beside `ukb-mesh/`, which it never modifies. It reuses only the `ukb-mesh` Python environment. The label table was rebuilt independently (STEP 0) and is identical to `ukb-mesh`'s.
2. **Geometry, not UVCs, for the LV endocardium.** The brief suggested LV nodes with ρ < 0.01, excluding nodes within 2 mm of the RV. Inspection showed instead that the boundary surface of all tets *including the valve plugs* has exactly 3 closed components: outer, RV cavity and LV cavity. Taking the LV cavity component directly gives the exact endocardium and an exact cavity volume. It also avoids the LV/RV junction UVC artefact. *Evidence:* female_age_65_bmi_26 gave regions of 358 / 120 / 114 mL; all 68 meshes checked have 3 regions.
3. **Plug geometry.** The aortic plug is ~3.3 mm thick and its interface band with the LV is ~2 mm tall; its cavity face domes ~1 mm into the LV. So LVOT depth 0 is set at the *ventricular edge of the ring band*, not at the ring centre plane.
4. **Exact plane slicing, not point slabs.** Sections are exact plane ∩ surface curves. VTK's `strip()` split closed loops at duplicate points, which made them look open. It was replaced by our own segment chaining (`chain_loops`).
5. **Finding: the closed LVOT tube is only 2–5 mm long.** At 1 mm steps, the strict closed section exists for 2–5 mm (60-mesh pilot: 2 mm in 8, 3 mm in 41, 4 mm in 11). Deeper sections open with a gap whose centre is 0–15° from the mitral direction. This is the missing anterior mitral leaflet: the meshes have no leaflets, and the aorto-mitral curtain is ~3 mm. The earlier probe's "closed for 5–10 mm" used coarse 30° angular coverage of points, which hid the gap.
6. **Primary LVOT = strict section at 2 mm.** It is closed in every pilot mesh.
7. **Secondary LVOT = chord-closed section at 4 mm.** Open sections are closed with a straight chord across the gap; rejected if the gap is > 120°. At 6 mm the gap was > 120° in 5 of 30 women and 0 of 30 men, which would make the missingness sex-biased, so 6 mm was not used. At 4 mm the gap was 0–85° and never rejected. *Caveat:* the chord approximates the leaflet and probably makes the mitral side slightly smaller and flatter than it is. *Evidence:* `outputs/figures/qc/gap_demo_female_age_65_bmi_26.png`.
8. **Ellipse method.** A direct least-squares (Fitzgibbon) ellipse fit gave ellipticity ~1.45–1.75 on D-shaped loops that looked much less oval. It was replaced by the **moment-equivalent ellipse** (same area and second moments), which is exact for true ellipses and robust otherwise. Feret calipers agree within ~0.5 mm. The direct fit is kept only to find the ring centre for angular binning.
9. **Ring outline.** Ring nodes form a band, so the outline is a 72-bin median-radius polygon around the fitted centre. Area and perimeter come from this polygon.
10. **QC thresholds.** Ring coverage ≥ 0.95; annulus 12–32 mm; mitral 20–55 mm; a loop more than 1.6 × the annulus diameter = opened into the cavity. In the subset plus pilot (68 meshes) the only flag was 1 `ring_gap`.
11. **LV size for indexing.** Diameter / LV cavity volume^(1/3), plus diameter as % of LV length (mitral centre → farthest endo point).
12. **Statistics.**
    - BMI held fixed in two ways: common support (paired F/M with the same age and BMI) and OLS adjustment (`age_c + sex + age_c×sex + BMI`, predicted at BMI 27).
    - "Narrows with age" = the age×sex coefficient has the opposite sign to the sex effect.
    - CIs and p-values are descriptive: meshes are averages with unknown participant counts.
13. **Colours.** Women `#D55E00` (vermillion), men `#0072B2` (blue), from Okabe-Ito; they pass the CVD validator (ΔE 21.9 protan).
14. **Pilot.** 30 random common-support pairs (seed 2026), 60 meshes, written into `lvot_all.csv` so the full run resumes from them. No errors; ~2.9 s/mesh → full run ≈ 70 min.

### 2026-10-24
15. **Reproducibility.**
    - All parameters moved to `config.yaml`; pinned `requirements.txt`; local `data/labels.csv`.
    - Every script logs its provenance (code and config hashes, package versions) to `outputs/logs/runs.jsonl`.
    - `# STEP x.y` markers in the code match `docs/PIPELINE.md` and `docs/flowchart.mmd`.
    - `run_all.sh` runs the pipeline end to end; git is used for versioning.
    - Re-running the subset after the refactor gave identical numbers (max abs diff 0.0).
16. **Cardiac phase: end-diastole only.** The request was to compare LVOT ellipticity in systole vs diastole. The dataset has one mesh per sex × age × BMI, and the `.vtu` holds no displacement or motion arrays (only elemTag, fiber, sheet). Ugurlu et al. (PLoS One 2025) state that the surface meshes were built from the end-diastolic frames; end-systolic frames were used only to validate volumes, and no ES or motion meshes are released (Zenodo 10.5281/zenodo.15649643). **All measurements here are end-diastolic; systolic LVOT shape cannot be measured from this dataset.** The missing anterior mitral leaflet also matters most in systole, when it forms the posterior LVOT wall.
    *Decision (user, 2026-10-24):* keep the analysis end-diastolic only and state this clearly as a caveat in the write-up. No end-systolic data will be sought.
17. **Area vs ellipse vs echo (STEP 6.6).**
    - **How each area is computed:**
      - outline area = shoelace formula on the outline;
      - ellipse estimates = π/4 × two diameters (moment axes, Feret max/min, or PLAX × orthogonal);
      - echo estimate = π/4·D_PLAX², assuming a circular LVOT.
    - **How D_PLAX is defined:** the edge-to-edge line through the outline's area centroid towards the mitral centre. This mimics the parasternal long-axis calliper (septum → aorto-mitral side). New columns: `*_d_plax`, `*_d_perp`, `*_plax_vs_minor_deg`, `lvot{d}cc_area`/Feret. The pilot and subset were regenerated with these columns (60/60 and 8/8 meshes OK).
    - **Pilot (30 common-support pairs, end-diastole):**
      - The moment ellipse is within +1% of the outline area. Feret and two-view ellipses are within about −6% to +6%.
      - The echo circle underestimates the area by about 9% at the annulus, ~30% at LVOT 2 mm (F −32.5%, M −27.7%) and ~32% at 4 mm (F −35.9%, M −28.9%).
      - The PLAX line runs within ~13–19° (median) of the LVOT minor axis, so echo mostly measures the short axis.
      - The echo underestimation is larger in women, by paired M − F +4.8 percentage points (95% CI +3.1 to +6.5) at 2 mm and +6.9 (+5.1 to +8.7) at 4 mm. There is no sex difference at the annulus.
    - **Caveat:** end-diastole only. Clinical LVOT diameters are measured in mid-systole, when the anterior mitral leaflet forms the posterior wall and the shape may differ.
18. **Full run (2026-10-24).**
    - 1,423/1,423 meshes, 0 errors, 67 min; all have 3 surface components.
    - LVOT 2 mm is missing for 10 meshes (closed tube only to 1 mm).
    - 62 aortic rings have outline coverage between 0.875 and 0.95 (polygon unaffected).
    - The 4 mm gap is ≥ 100° in 8 meshes (6 women, 2 men); excluding them does not change ellipticity (F 1.544 → 1.543, M 1.450 → 1.449).
    - Robust-z outliers (MAD-based, inflated by the small spread of averaged meshes): the top three were rendered (female_age_73_bmi_43, male_age_66_bmi_42, female_age_59_bmi_42). They are genuine shapes; none were excluded.
19. **Bug fix in `03_analysis.py` `roll()`.** `pd.Series(series, index=age)` re-indexed by row label instead of re-labelling, so the rolling-mean *lines* in `fig_raw` and `fig_common` were wrong. Dots, tables and models were unaffected. Fixed by passing numpy arrays. Also added `lvot4cc_d_area_per_edv13` to the analysed measures.
20. **CT-style ellipse added (STEP 5.6, 6.6, 6.7).**
    - **Definition:** the longest edge-to-edge line through the outline's area centroid, and the line through the centroid at 90° to it (`ct_diams()`, 0.5° search); area = π/4 × longest × perpendicular. Unlike the Feret widths, both lines go through the centre and are forced to 90°, as in a typical CT reading.
    - **Saved outlines:** the outlines had not been saved, so `09_save_outlines.py` now writes per-mesh 2D outlines to `data/outlines/` (resumable). Saved outlines reproduce the table's area and PLAX diameter exactly.
    - **Pilot (68 meshes):** the CT ellipse is close to the Feret ellipse (annulus +5.3% F / +4.0% M; LVOT 2 mm −3.3 / −2.8%; 4 mm −3.6 / −2.5%). Its longest diameter is within 0.3 mm of the widest calliper width and its perpendicular within 0.2 mm of the narrowest. Forcing 90° changes little because the outlines are close to elliptical.
    - **Full cohort (2026-10-25):** `09_save_outlines.py --all` ran for all 1,423 meshes (67 min, 0 errors); saved outline areas match the table exactly, and the only missing outlines are the 10 meshes without a strict 2 mm section. CT ellipse vs outline area: annulus +5.2% F / +4.2% M; LVOT 2 mm −3.5 / −2.8%; LVOT 4 mm −3.8 / −2.4% (Feret: +5.8 / +5.3, −2.0 / −1.6, −3.3 / −2.1).
21. **Figure scope (user, 2026-10-25).** The area-methods figure shows four methods (moment, Feret, CT, echo circle). The two-view (PLAX × perpendicular) ellipse and the PLAX-plane sensitivity analysis (±10–20° rotation, 2 mm off-centre) are kept for later analyses, not presented now. The two-view code stays in `07_area_methods_demo.py` (`SHOW_TWO_VIEW`) and in `06_area_methods.py` / `area_methods_summary.csv`. Bars are coloured by level (yellow annulus, green 2 mm, purple 4 mm), semi-transparent, women lighter and men darker (hatching removed at the user's request); main title "Estimating the elliptical LVOT area"; the formula boxes were removed and each estimate's area and % error moved into the legend.
22. **Naming in the figure (user, 2026-10-25).** The "CT ellipse" is now called the **minor–major axis ellipse** (longest diameter through the centre × the one at 90°; code and table columns keep the `ct_` / `ellipse_ct` names). Panel 4 is named **Echocardiographic area** (was "echo circle") and described only as "approximation of the echo PLAX view". Method names are in bold `#102D40` in the panel titles and on the bar labels; the assumed ellipses/circle are drawn in `#B94F00`.
23. **Font (user, 2026-10-25).** The area-methods figure uses **Roboto**. Static TTFs (Regular, Medium, Bold, Italic) and the Apache 2.0 licence from Google's official release (`github.com/googlefonts/roboto`, v2.138, `roboto-unhinted.zip`) are stored in `assets/fonts/` and registered by `07_area_methods_demo.py`, so nothing is installed system-wide and the figure renders identically elsewhere. The "mitral" label was moved inside each panel.
24. **PLAX plane definition (2026-10-25, deferred analysis).**
    - **Mitral centre:** the centre of the ellipse fitted to the mitral rim nodes (LV nodes touching the mitral plug), projected onto their best-fit plane.
    - **Uniqueness:** with that definition, the PLAX plane (aortic axis + mitral centre) is unique per mesh. Other reasonable definitions give different planes.
    - **Comparison on 68 meshes (subset + pilot), LVOT 2 mm:**
      - B (aortic axis + LV apex) is a median 5° from A (max 19°).
      - C (aortic centre + mitral centre + apex) is a median 3° from A (max 11°).
      - All three planes lie within 11–14° (median) of the LVOT minor axis.
      - Echo-circle error, median, women / men: A −31.9 / −27.6%, B −33.4 / −28.1%, C −32.8 / −27.7%.
    - **Conclusion:** the ~30% underestimation and its sex difference do not depend on the definition.
    - **Status:** kept for further analyses, not presented now (with #21). `11_plax_definitions.py`, `fig_plax_definitions.png`.
25. **Which plug is the aortic valve? (2026-10-25)**
    - **Source:** the files label regions only by number (elemTag 0, 1, 4–7). Ugurlu et al. 2025 list six regions (LV, RV, aortic, tricuspid, pulmonary and mitral valves) but give no number-to-region map.
    - **Method:** the identities were inferred from anatomy and checked on 40 meshes (`15_plug_identity.py`, `plug_identity.csv`).
    - **Which ventricle each plug touches:** tags 4 and 5 touch only the LV, and tags 6 and 7 only the RV, in 40/40 meshes.
    - **Tag 5 vs tag 4 (the two LV plugs):**
      - Tag 5 is the smaller one: median 19.4 mm vs 35.3 mm (40/40).
      - Tag 5 is oblique to the LV long axis (59°), while tag 4 is nearly perpendicular to it (6°). This is the expected outflow vs inflow valve orientation.
      - Tag 5 is nearer the RV outflow plug, tag 7 (32 vs 51 mm; 40/40).
    - **RV plugs:** tag 7 (22 mm, oblique 27°) matches the pulmonary valve, and tag 6 (34 mm, perpendicular 4°) matches the tricuspid.
    - **LVOT gap:** the open side of the LVOT sections faces tag 4 in every mesh (gap 0–15° from its direction), consistent with aorto-mitral continuity.
    - **Conclusion:** tag 5 = aortic, 4 = mitral, 6 = tricuspid, 7 = pulmonary. This is an inference, not a documented label; the authors' code may confirm it.
26. **"What is a mesh?" slide (2026-10-25).** `16_mesh_explainer.py` → `outputs/figures/mesh_explainer.png` (4 panels: whole heart by region; surface triangles; slice showing tetrahedra filling the wall; one tetrahedron), and `17_mesh_slide.js` → `outputs/slide_what_is_a_mesh.pptx` (full-slide figure with speaker notes).
27. **One comparison set for all sex comparisons (user, 2026-10-25).**
    - **Change:** the ellipticity (STEP 6.8) and area-method (STEP 6.6, 6.7) summaries and figures now use the matched same-age, same-BMI pairs (622; 617 at 2 mm, where both meshes must have a closed section), like the size figures. Previously they used all meshes.
    - **Effect:** values changed by at most 0.01 in ellipticity and 0.2 percentage points in the echo error (e.g. 4 mm ellipticity in women 1.55 → 1.54; echo error at 2 mm −32.7 → −32.5% women, −28.6 → −28.8% men).
    - **Unchanged:** `area_methods_by_mesh.csv` keeps all meshes, with a `paired` flag. The adjusted regression still uses all 1,423 meshes, and `fig_raw` still shows all meshes.
    - **Rebuilt:** the write-up tables and the highlights deck.
28. **AVA with a circular vs elliptical LVOT (2026-10-25).**
    - **Principle:** with the same LVOT and AV VTIs, AVA scales with the LVOT area used, so the correction factor is k = elliptical (measured outline) area / circular area π/4·D_PLAX².
    - **k, median (IQR), matched pairs:**
      - annulus 1.10 (1.08–1.12) in both sexes;
      - LVOT 2 mm: women 1.48 (1.44–1.54), men 1.40 (1.35–1.46);
      - LVOT 4 mm: women 1.56 (1.50–1.63), men 1.42 (1.35–1.51).
      - The minor–major axis ellipse gives k 1.43 (women) and 1.36 (men) at 2 mm.
    - **Equivalent cut-off:** a corrected AVA of 1.0 cm² corresponds to a circular AVA of about 0.67 cm² (women) and 0.71 cm² (men) at 2 mm, and 0.91 cm² at the annulus.
    - **Reclassification at 2 mm:** a circular AVA of 0.8 cm² becomes about 1.19 cm² (women) and 1.12 cm² (men), so about 97–100% of hearts would no longer be "severe". At the annulus almost none would be (1%).
    - **Size-independent approach:** the mesh LVOTs are end-diastolic and smaller than clinical systolic LVOTs, so a fixed-VTI sweep on absolute mesh AVAs was rejected. The simulation applies the unitless k to a range of circular AVAs instead.
    - **Public data:** no public per-patient dataset with LVOT diameter, LVOT VTI and AV VTI was found (searched 2026-10-25). Open AS echo datasets (e.g. TMED-2, EchoNext) give images or severity grades only. MIMIC echo data need credentialed access.
    - **Caveat:** AVA cut-offs were validated with the circular method, so a correction like this needs its own outcome-based thresholds.

