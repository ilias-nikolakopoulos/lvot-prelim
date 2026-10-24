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

