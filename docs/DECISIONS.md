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
