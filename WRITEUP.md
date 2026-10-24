# LVOT size and ellipticity by age and sex — UK Biobank representative meshes

> **All results describe the heart at end-diastole.** The released meshes were built from end-diastolic cine frames only (Ugurlu et al., PLoS One 2025). The dataset contains no systolic meshes and no motion data, so how LVOT size or ellipticity changes between diastole and systole, in men or women, **cannot be assessed here**.

_Status: complete (2026-10-24). Full run on all 1,423 meshes, 0 processing errors._

## Methods (summary)
Full detail is in `docs/PIPELINE.md`; reasons for each choice are in `docs/DECISIONS.md`.
- **Data:** 1,423 representative end-diastolic biventricular meshes. Each is the average of the UK Biobank participants with the same sex, age (46–83) and integer BMI.
- **LV size:** the LV cavity is the closed surface made of the LV myocardium and the mitral and aortic valve plugs. Its volume = LV end-diastolic volume. LV length = mitral ring centre to apex.
- **Aortic annulus:** LV nodes shared with the aortic plug, fitted with a plane. Reported: area- and perimeter-derived diameters, and moment-ellipse major/minor diameter and ellipticity.
- **LVOT:** exact cross-sections perpendicular to the aortic axis, at depths measured from the ventricular edge of the annulus.
  - Primary: the closed section at 2 mm.
  - Secondary: the 4 mm section, where the open mitral side is bridged by a straight chord.
- **Other measures:** mitral ring diameter, aorto-mitral plane angle, and the angle between the aortic axis and the LV long axis.
- **Area methods:** the outline area (shoelace formula) is compared with ellipse estimates (π/4 × two diameters) and with the standard echo estimate π/4·D_PLAX². D_PLAX is the diameter along the parasternal long-axis direction, from the septum to the aorto-mitral side.
- **Statistics:** measure vs age by sex, holding BMI fixed in two ways:
  - common support: men and women with the same age and BMI (622 pairs);
  - regression: y ~ age + sex + age×sex + BMI.

  Sizes are reported in mm and indexed to LV size (÷ EDV^(1/3), and as % of LV length).

## Results
All 1,423 meshes were measured, with 0 errors; every mesh had the expected 3 surface components.
- **Strict closed LVOT tube:** 3 mm long in most meshes (1 mm: 10; 2 mm: 212; 3 mm: 1,025; 4 mm: 163; 5 mm: 13).
- **Primary LVOT (2 mm):** available for 1,413 meshes. The 10 missing ones (6 women, 4 men) have a closed tube only to 1 mm.
- **Secondary LVOT (4 mm, chord-closed):** available for all 1,423 meshes. The gap is larger in women (median 60° vs 45°) and is ≥ 100° in only 8 meshes; excluding those does not change the results.

Typical values are medians across meshes:

| End-diastole | Women | Men |
|---|---|---|
| LV cavity volume (EDV) | 116 mL | 150 mL |
| LV length | 83 mm | 92 mm |
| Aortic annulus, area-derived diameter | 17.7 mm | 20.5 mm |
| LVOT 2 mm, area-derived diameter (major × minor) | 18.5 mm (22.5 × 15.4) | 21.4 mm (25.5 × 18.3) |
| LVOT 4 mm (chord-closed), area-derived diameter | 20.3 mm | 23.2 mm |
| Annulus ellipticity (major/minor) | 1.14 | 1.14 |
| LVOT 2 mm ellipticity | 1.46 | 1.39 |
| LVOT 4 mm ellipticity | 1.55 | 1.44 |
| Mitral ring diameter | 34.1 mm | 37.3 mm |
| Aorto-mitral plane angle | 60° | 63° |

### Sex difference by age (BMI held fixed: 622 same-age, same-BMI pairs)
Paired men − women, mean (95% CI):

| Measure | 46–54 | 55–59 | 60–64 | 65–69 | 70–74 | 75–83 |
|---|---|---|---|---|---|---|
| Annulus diameter, mm | +2.88 (2.73–3.04) | +2.80 | +2.75 | +2.73 | +2.60 | +2.62 (2.45–2.80) |
| … as % of women's | +16.2% | +15.8% | +15.5% | +15.4% | +14.7% | +14.9% |
| LVOT 2 mm diameter, mm | +2.92 (2.78–3.05) | +2.86 | +2.91 | +2.87 | +2.83 | +2.87 (2.75–2.99) |
| LVOT 4 mm diameter, mm | +2.88 (2.71–3.06) | +2.78 | +2.87 | +2.96 | +2.96 | +2.97 (2.81–3.12) |
| Annulus ellipticity | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | −0.01 |
| LVOT 2 mm ellipticity | −0.06 (−0.07 to −0.05) | −0.07 | −0.06 | −0.05 | −0.06 | −0.06 (−0.07 to −0.04) |
| LVOT 4 mm ellipticity | −0.09 (−0.10 to −0.07) | −0.09 | −0.08 | −0.08 | −0.09 | −0.09 (−0.10 to −0.08) |
| Annulus / EDV^⅓ | +7.2% | +6.5% | +6.1% | +6.1% | +5.1% | +5.4% |
| LVOT 2 mm / EDV^⅓ | +6.5% | +6.1% | +6.2% | +6.2% | +5.7% | +6.2% |

The full table, with CIs for every block, is in `outputs/tables/sexdiff_by_age_block.csv`.

**Size (absolute)**
- Men's annulus and LVOT are about **2.7–2.9 mm (15–16%) larger** than women's at every age.
- **Annulus:** the gap narrows slightly with age, from 2.9 to 2.6 mm. The adjusted age×sex interaction is p = 0.002, because men's annulus shrinks a little with age (−0.016 mm/yr) while women's barely changes (−0.005 mm/yr).
- **LVOT:** the gap does **not** narrow (interaction p = 0.34 at 2 mm, 0.52 at 4 mm). Both sexes lose about 0.7 mm of LVOT diameter between their 50s and late 70s (2 mm: −0.024 mm/yr in women, −0.027 mm/yr in men).

**Size (indexed to LV size)**
- The LV cavity shrinks faster with age than the LVOT (EDV −0.7 to −0.8 mL/yr), so the LVOT and annulus get relatively larger with age, more so in women.
- Indexed to EDV^⅓, men's outflow is still about 5–7% larger, so the sex gap is not fully explained by LV size.
- The indexed annulus gap narrows with age (7.2% → 5.4%, p = 0.004). The indexed LVOT gap stays at about 6%.

**Shape (ellipticity)**
- The **annulus is nearly round** (major/minor about 1.14–1.15), with no sex or age difference.
- The **LVOT just below is clearly oval**: about 1.46 in women and 1.39 in men at 2 mm, and about 1.55 vs 1.44 at 4 mm.
- **Women's LVOT is more elliptical** at every age (−0.06 at 2 mm and −0.09 at 4 mm, men − women). This difference is constant with age (interaction p = 0.78 and 0.87).
- Ellipticity falls only slightly with age (about −0.001 per year in both sexes, ≈ −0.04 over the age range).
- The short axis runs from the aorto-mitral side to the septum; the long axis runs side to side.

**Angles**
- The aorto-mitral angle (59° → 61° in women, 62° → 65° in men) and the angle between the aortic axis and the LV long axis (54° → 59° in women, 57° → 63° in men) both **widen with age**, by about 4–6° from the youngest to the oldest block.
- Both angles are 3–4° larger in men.

### Outline area vs ellipse and echo estimates (`outputs/tables/area_methods_summary.csv`)
% difference from the measured outline area:

| Area estimate | Annulus F / M | LVOT 2 mm F / M | LVOT 4 mm F / M |
|---|---|---|---|
| Ellipse, moment axes π/4·a·b | +0.7 / +0.8% | +0.9 / +0.9% | +0.6 / +0.5% |
| Ellipse, Feret max × min | +5.8 / +5.3% | −2.0 / −1.6% | −3.3 / −2.1% |
| Ellipse, PLAX × orthogonal | +3.1 / +2.7% | −5.1 / −4.7% | −6.6 / −5.2% |
| **Echo circle π/4·D_PLAX²** | **−9.4 / −9.3%** | **−32.7 / −28.6%** | **−36.0 / −29.8%** |

- **Ellipse estimates** are within about 7% of the true area: the LVOT is close to elliptical.
- **Echo circle:** the PLAX diameter lies within about 12–19° (median) of the LVOT minor axis, so a circular assumption based on it **underestimates LVOT area by about 30%**. This flows directly into stroke volume and continuity-equation valve area.
- **By sex:** the underestimation is larger in women. In paired comparisons it is +3.8 percentage points (95% CI 3.5–4.1) at 2 mm and +5.8 (5.5–6.2) at 4 mm. There is no sex difference at the annulus (−0.2, −0.5 to +0.1).
- **By age:** the underestimation shrinks slightly with age (+1.6 to +1.9 percentage points per decade), in line with the slightly rounder LVOT in older hearts.
- **Caveat:** this is at end-diastole; echo LVOT diameters are measured in mid-systole.

### Figures (`outputs/figures/analysis/`)
- `fig_raw.png`: all meshes.
- `fig_common.png`: common support, BMI held fixed.
- `fig_adjusted.png`: regression at BMI 27.
- `fig_sexdiff.png`: men − women vs age.
- `fig_area_methods.png`: area estimates.
- `fig_area_methods_demo.png`: what each area estimate measures, drawn on one outline.
- QC renders are in `outputs/figures/qc/`. These include `outline_steps_*` (every step from mesh to outline), `outline_demo_*` (how an outline is obtained) and `gap_demo_*` (chord closure).

The edge ages (46, 82–83) have few pairs and are noisy in the per-age plots.

## Caveats
1. **End-diastole only.** Every measurement is end-diastolic. Systolic LVOT dimensions and shape, and their diastole-to-systole change, are not available in this dataset. Clinical LVOT measurements for valve sizing and flow calculations are usually taken in **mid-systole**, so the absolute values here are not directly comparable with them. The LVOT is also known to change shape through the cycle.
2. **No valve leaflets.** The meshes have no mitral leaflets, so the anterior mitral leaflet (the posterior LVOT wall) is missing. The LVOT is a fully closed tube for only 2–5 mm below the annulus. Deeper sections are bridged with a straight chord, which approximates the leaflet and probably makes the mitral side slightly smaller and flatter than it really is. This matters most in systole, which is not represented here anyway.
3. **Averaged hearts.** Each mesh averages many participants with unknown counts. Results describe typical hearts, not individual variation. Averaging tends to smooth and slightly shrink shapes. CIs and p-values are descriptive only.
4. **Annulus definition.** The annulus is where the aortic plug meets the LV myocardium. This may sit inside the true hinge points, so annulus sizes may be smaller than CT-derived values. Comparisons by sex and age are more reliable than absolute values.
5. **Robust outliers.** Model residuals flag some meshes, mostly at extreme BMI, where fewer participants go into each average. The top ones were rendered (`outputs/figures/qc/`) and are genuine shapes, not method failures, so none were excluded beyond the 10 meshes without a closed 2 mm section.
6. **Cross-sectional data.** Age differences are between age groups, not changes within a person over time. Survivor and cohort effects are possible.
