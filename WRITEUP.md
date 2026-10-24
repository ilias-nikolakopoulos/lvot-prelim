# LVOT size and ellipticity by age and sex — UK Biobank representative meshes

> **All results describe the heart at end-diastole.** The released meshes were built from end-diastolic cine frames only (Ugurlu et al., PLoS One 2025). The dataset contains no systolic meshes and no motion data, so how LVOT size or ellipticity changes between diastole and systole, in men or women, **cannot be assessed here**.

_Status: draft. Methods and caveats are final. Results are filled in after the full 1,423-mesh run (`scripts/03_analysis.py`). Pilot numbers (60 meshes) are in `docs/DECISIONS.md`._

## Methods (summary)
Full detail is in `docs/PIPELINE.md`; reasons for each choice are in `docs/DECISIONS.md`.
- **Data:** 1,423 representative end-diastolic biventricular meshes. Each is the average of the UK Biobank participants with the same sex, age (46–83) and integer BMI.
- **LV size:** the LV cavity is the closed surface made of the LV myocardium and the mitral and aortic valve plugs. Its volume = LV end-diastolic volume. LV length = mitral ring centre to apex.
- **Aortic annulus:** LV nodes shared with the aortic plug, fitted with a plane. Reported: area- and perimeter-derived diameters, and moment-ellipse major/minor diameter and ellipticity.
- **LVOT:** exact cross-sections perpendicular to the aortic axis, at depths measured from the ventricular edge of the annulus.
  - Primary: the closed section at 2 mm.
  - Secondary: the 4 mm section, where the open mitral side is bridged by a straight chord.
- **Other measures:** mitral ring diameter, aorto-mitral plane angle, and the angle between the aortic axis and the LV long axis.
- **Statistics:** measure vs age by sex, holding BMI fixed in two ways:
  - common support: men and women with the same age and BMI (622 pairs);
  - regression: y ~ age + sex + age×sex + BMI.

  Sizes are reported in mm and indexed to LV size (÷ EDV^(1/3), and as % of LV length).

## Results
_To be completed after the full run._
- Per-mesh table: `outputs/tables/lvot_all.csv`
- Sex differences by age block: `outputs/tables/sexdiff_by_age_block.csv`
- Model coefficients: `outputs/tables/model_coefs.csv`
- Figures: `outputs/figures/analysis/` (raw, common support, adjusted, sex difference)

## Caveats
1. **End-diastole only.** Every measurement is end-diastolic. Systolic LVOT dimensions and shape, and their diastole-to-systole change, are not available in this dataset. Clinical LVOT measurements for valve sizing and flow calculations are usually taken in **mid-systole**, so the absolute values here are not directly comparable with them. The LVOT is also known to change shape through the cycle.
2. **No valve leaflets.** The meshes have no mitral leaflets, so the anterior mitral leaflet (the posterior LVOT wall) is missing. The LVOT is a fully closed tube for only 2–5 mm below the annulus. Deeper sections are bridged with a straight chord, which approximates the leaflet and probably makes the mitral side slightly smaller and flatter than it really is. This matters most in systole, which is not represented here anyway.
3. **Averaged hearts.** Each mesh averages many participants with unknown counts. Results describe typical hearts, not individual variation. Averaging tends to smooth and slightly shrink shapes. CIs and p-values are descriptive only.
4. **Annulus definition.** The annulus is where the aortic plug meets the LV myocardium. This may sit inside the true hinge points, so annulus sizes may be smaller than CT-derived values. Comparisons by sex and age are more reliable than absolute values.
5. **Cross-sectional data.** Age differences are between age groups, not changes within a person over time. Survivor and cohort effects are possible.
