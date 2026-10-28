# What makes the aorto-mitral angle wider? A literature review (2026-10-27)

Sources: PubMed, searched 2026-10-27. The question: in our meshes the aorto-mitral angle is the strongest correlate of LVOT ellipticity (DECISIONS #35–36), so what does the literature say determines it?

## 1. First, which way is "wider"?

**Our angle and the clinical angle point opposite ways.**
- `ao_mv_angle` ([`src/lvot.py`](../src/lvot.py), STEP 3.5) is the angle between the two ring-plane normals, both pointing into the LV. If the two valves lay flat in one plane, it would be 0°.
- The clinical aorto-mitral angle (AMA), on 3D echo or CT, is the angle between the aortic annulus plane and the mitral annulus plane, measured through the LV at the aorto-mitral curtain. A flat base would be 180°.
- So **clinical AMA ≈ 180° − `ao_mv_angle`.**

**In our meshes, clinical AMA is:**
- Women: median 120.5° (IQR 119.2–121.9).
- Men: 116.7° (115.4–118.1).
- This is in the range reported clinically: 115–140°, with a normal SD of about 10° ([Varghese 2013](https://doi.org/10.1093/ejcts/ezt234); [Thériault-Lauzier 2016](https://doi.org/10.4244/EIJY15M11_09); [Meimoun 2024](https://doi.org/10.1111/echo.15898)).

**Our finding in clinical terms:** a **wider (more obtuse, flatter) clinical AMA goes with a more elliptical LVOT**, and a narrower (more acute) AMA with a rounder one.
- **Per 1 SD of angle:** −0.39 SD ellipticity per SD of `ao_mv_angle`, i.e. +0.39 SD per SD of clinical AMA.
- **Labels to fix:** the "narrow / low angle" labels in DECISIONS #35–36 and `wall_all_phenotypes.csv` mean low `ao_mv_angle`, which is a *wide* clinical AMA.

**Spread is much smaller than in patients:** the SD across our representative meshes is only 2.6° within each sex, against about 10° between individual patients. The meshes are cohort averages, so most individual variation is averaged out.

## 2. What the literature links to a wider (more obtuse) AMA

No study we found models the AMA's determinants in a healthy population. No study links the AMA to LVOT ellipticity. Most of the evidence is cross-sectional, from valve-disease and procedural cohorts.

| Factor | Direction | Evidence |
|---|---|---|
| **LV dilatation / ventricular remodelling** | wider | Ventricular functional MR had the largest AMA among MR types ([Hirano 2026](https://doi.org/10.1007/s11748-026-02272-0)). Non-ischaemic FMR (larger annulus, longer leaflets) had a wider AMA than ischaemic FMR: 141° vs 135° ([Aruta 2018](https://doi.org/10.1186/s12947-018-0145-8)). In dogs with myxomatous MV disease, AMA flattens with stage: 122° → 134° → 141° ([Mizuno 2021](https://doi.org/10.1016/j.jvc.2021.01.004)). |
| **Larger, flatter mitral annulus** (degenerative MR, fibroelastic deficiency) | wider | FED patients had a larger AMA, annulus and annular height than controls ([Ma 2018](https://doi.org/10.1007/s10554-018-1383-4)). The annulus is larger in MV pathology, AF and LV/LA dilatation ([Papadopoulos 2022](https://doi.org/10.1111/echo.15481)). |
| **Larger LV cavity, not hyperdynamic** | wider (the converse: small, hypercontractile LV ↔ narrow) | A narrow AMA (<120°) clustered with LV end-diastolic diameter <45 mm in SAM risk ([Varghese 2013](https://doi.org/10.1093/ejcts/ezt234)). Small LVESV, high EF and small AMA went together in LVOT narrowing after TMVI ([Hayashi 2021](https://doi.org/10.1111/echo.15078)). Patients with LV obstruction after TAVR were smaller, older and more often women, with a more acute AMA: 107° vs 114° ([Ferri 2025](https://doi.org/10.1016/j.ijcard.2024.132956)). |
| **Removing a stenotic aortic valve** | wider | AMA rose from 133° to 142° immediately after surgical AVR, with less aorto-mitral flexion ([Warraich 2013](https://doi.org/10.1016/j.amjcard.2013.06.013)). Low-gradient AS had a narrower AMA than high-gradient AS: 104–106° vs 120° on CT ([See 2024](https://doi.org/10.1016/j.ijcard.2024.132174)). |
| **Left atrial geometry** | wider after LA plication | Ring annuloplasty narrows the AMA. Adding LA plication prevents this or widens it ([Hirano 2026](https://doi.org/10.1007/s11748-026-02272-0); [Takahashi 2022](https://doi.org/10.1111/jocs.16355); [Takahashi 2024](https://doi.org/10.21037/jtd-24-593)). A larger post-repair AMA correlated with LA reverse remodelling (r = 0.63; [Ma 2018](https://doi.org/10.1007/s10554-018-1383-4)). |
| **Short aorto-mitral curtain** | wider | After Commando reconstruction, curtain length was roughly inversely related to AMA ([Simpson 2023](https://doi.org/10.1093/ejcts/ezad155)). |
| **Cardiac phase / contractility** | ~1–11° change | Aorto-mitral flexion of 8°, rising to 11° with inotropes, in sheep ([Timek 2003](https://doi.org/10.1016/s0003-4975(03)01078-6)). FMR patients: 130° in systole vs 131° in diastole on CT ([Thériault-Lauzier 2016](https://doi.org/10.4244/EIJY15M11_10)). |
| **Body size** | no association (animal data) | Pigs: AMA did not change with body weight ([Lipiski 2020](https://doi.org/10.1038/s41598-020-63044-1)). |

**Conditions with a narrow (acute) AMA:**
- **HCM:** obstructive HCM ~97–100° ([Moreno Garijo 2023](https://doi.org/10.1016/j.atssr.2023.09.016)). A smaller AMA independently predicted LVOT obstruction (OR 0.89 per degree; [Lin 2024](https://doi.org/10.1016/j.echo.2024.05.002)) and non-response to mavacamten: 115° vs 126° ([Ha 2026](https://doi.org/10.1093/ehjci/jeag156)).
- **SAM:** a narrow AMA marks SAM risk after MV repair, during dobutamine stress (126° vs 139°; [Meimoun 2024](https://doi.org/10.1111/echo.15898)) and with MAC ([Friend 2022](https://doi.org/10.1016/j.echo.2022.12.011)).
- **Aortic dilatation:** a case report of reduced AMA with aortic dilatation ([Donia 2023](https://doi.org/10.1016/j.jaccas.2023.102191)).
- **Mechanism:** in an ex vivo porcine model, widening the AMA increased coaptation–septal distance and improved outflow haemodynamics ([Park 2022](https://doi.org/10.1007/s10439-022-03091-z)).

## 3. The related aorto-septal angle (age, sex, hypertension)

The aorto-septal angle has the best population data. In our meshes it correlates with clinical AMA (ρ about +0.5 within sex), so its determinants are a reasonable proxy.

- **Age:** narrows with age ([Kong 2022](https://doi.org/10.1016/j.ijcard.2022.04.046), 140 adults aged 20–90 without structural heart disease; cats: 0.55° per year, [Crosara 2023](https://doi.org/10.1016/j.jvc.2023.03.002)).
- **Sex:** women have a wider angle before 70, but narrow faster, so the difference is gone beyond 70 ([Kong 2022](https://doi.org/10.1016/j.ijcard.2022.04.046)).
- **Hypertension and aortic diameter:** the angle narrows with hypertension duration and ascending-aorta diameter, independent of age (n = 1,294; [Balaban 2022](https://doi.org/10.1536/ihj.21-471)).
- **Mechanism proposed:** age- and pressure-related aortic elongation pushes the aortic root towards the LV.
- **Linked to septal changes:** a narrower angle goes with a basal septal bulge (cats; [Novo Matos 2023](https://doi.org/10.1016/j.jvc.2023.09.006)), and with events in haemodialysis patients ([Nakayama 2024](https://doi.org/10.1371/journal.pone.0298637)).

## 4. What our meshes show (Spearman ρ with clinical AMA, within sex, F / M)

Using `lvot_all.csv` and `wall_all.csv`, without z-scoring.

- **Wider AMA goes with:**
  - larger LV EDV (+0.46 / +0.49)
  - greater aortic–mitral centre distance (+0.40 / +0.31)
  - longer LV (+0.37 / +0.37)
  - larger LV myocardial volume (+0.37 / +0.34)
  - larger mitral annulus (+0.35 / +0.33)
  - higher BMI (+0.32 / +0.23)
  - wider aorto-septal angle (+0.48 / +0.50)
- **Narrower AMA goes with:**
  - older age (−0.34 / −0.44); the median falls about 2° from the 40s to the 70s in each sex
  - a larger aortic-axis vs LV-long-axis angle (−0.41 / −0.55)
- **Little or no link:** aortic annulus diameter (−0.11 / −0.08) and septal thickness (≈0).
- **Sex:** women have a wider AMA than men (120.5° vs 116.7°).

**Fit with the literature:** a larger, longer, younger LV with a bigger mitral annulus has a flatter base (wider AMA). Ageing, hypertension-related aortic elongation and a smaller, hyperdynamic LV fold the base (narrower AMA).

**Fit with the LVOT finding:** a wider AMA (as in women, younger and larger-LV meshes) goes with a more elliptical LVOT. A narrower one (older age, HCM- or SAM-type geometry) goes with a rounder LVOT.

## 5. Gaps and what it means for the HCM idea

- **No population reference values:** no study gives determinants of AMA in a general population, and none relates AMA to LVOT shape. UK Biobank CMR could do both. A 3-chamber cine gives the AMA; LVOT ellipticity needs a short-axis or 3D LVOT dataset.
- **Opposite direction in HCM:** obstructive HCM has a *narrow* AMA. Our data say a narrow AMA goes with a *rounder* LVOT. If LVOT ellipticity were an early HCM marker, it would most likely be *lower* ellipticity together with a narrow AMA, not higher. This is a hypothesis to test in individual scans; the averaged meshes can't settle it.
- **Phase:** the AMA changes by a few degrees through the cycle. Our meshes are end-diastolic, while many clinical values are systolic.

## Search record

- **Queries:**
  - `aorto-mitral angle[Title/Abstract]` (26 hits)
  - `aortomitral angle OR "aortic-mitral angle" OR "mitral-aortic angle" OR aortomitral angulation` (65 hits; first 40 screened)
  - aorto-septal angle with age or hypertension (22 + 25 hits; screened by title)
  - aortic root angulation or unfolding with age, hypertension or sex, on CT or MR (70 hits; abstracts of 10 were requested but not returned, so they were not reviewed)
- **Abstracts read:** 51, all in full.
- **Limits:**
  - Not systematic: one database (PubMed), English, abstracts only.
  - First authors and years were checked against the PubMed records.
