# Rationale: why the LVOT ends at the end of the aorto-mitral curtain

This file gives the reasons for the rule in [README_lvot_curtain.md](README_lvot_curtain.md). The README tells you how to apply the rule. This file tells you why the rule is as it is.

**How to read this file.** Each reason has one of three labels:
- **Recorded:** the reason is written in the project files (README, script, notes, tables).
- **Checked:** the statement was tested on the data, and the numbers are given.
- **Not recorded:** the rule was developed on 2026-11-10 and 2026-11-11, and the discussion from those two days is not in the project files. The value is known, but the reason is not written down. Where a reason is probable, it is given after the word "Probable". Treat these as suggestions, not as decisions. The author must confirm or correct them.

---

## 1. The anatomy

- The LVOT is the channel from the left ventricle to the aortic valve.
- **Front and side walls:** muscle (the ventricular septum and the adjacent LV wall).
- **Back wall:** not muscle. It has two parts, one after the other:
  1. **The aorto-mitral curtain.** This is a sheet of fibrous tissue between the aortic valve and the mitral valve. It is fixed and does not move during the heartbeat. The left atrium (LA) is immediately behind it.
  2. **The anterior mitral leaflet.** It hangs from the lower edge of the curtain. The line where it attaches is the hinge. The leaflet moves: it swings into the LV in diastole and closes in systole.
- **The end of the curtain is the hinge.** Above the hinge, the back wall of the LVOT is a fixed structure. Below the hinge, the "wall" is a moving leaflet, and in diastole there is no back wall at all, only the open mitral orifice.
- **So the hinge is the last level where the LVOT is a tube with a wall on all sides.** This is why the rule uses it as the end of the LVOT.
- References for the anatomy (standard texts, added with this file; they were not used to set the parameters):
  - Piazza N, et al. Anatomy of the aortic valvar complex and its implications for transcatheter implantation of the aortic valve. *Circ Cardiovasc Interv* 2008;1:74–81. PMID 20031657.
  - Ho SY. Structure and anatomy of the aortic root. *Eur J Echocardiogr* 2009;10:i3–i10. PMID 19131496.

**Link to the UK Biobank part of the project (recorded).** The UK Biobank meshes have no leaflets. In those meshes, the LVOT is a closed tube for only 2–5 mm below the annulus, and then it opens on the mitral side (`lvot-ukb/docs/DECISIONS.md` #5). The CT rule looks for the same anatomical limit, but on images that show the leaflet.

---

## 2. Why the wall-to-LA gap is the signal

**What the rule measures.** In a long-axis plane, the back wall of the LVOT is a thin dark line between two bright, contrast-filled spaces: the LVOT in front and the LA behind. The rule follows this dark line downward. At each depth, it measures the distance from the dark line to the front edge of the LA label. This distance is the gap.

**Why the gap shows the hinge (recorded in the README, section 2.3 and the NOTE in 2.8).**
- **Above the hinge,** the dark line is the curtain. The LA is directly behind it. The gap is approximately zero.
- **Below the hinge in diastole,** the dark line is the open anterior mitral leaflet. The leaflet points into the LV, away from the LA. The gap opens and continues to increase with depth.
- **The depth where the gap starts to open is the hinge,** and thus the end of the curtain.

**Why it works in diastole but not in systole (recorded).**
- In systole, the mitral valve is closed. The closed leaflet stays against the LA. The gap does not open, so the rule finds no end.
- The rule does not try to give a result for these scans. It records "not found", and systolic scans are flagged (LV/myocardium volume ratio < 0.6) for a manual mark.
- **Checked:** 14 of 998 scans have the systolic flag. The rule found no end in 9 of them. It gave an end in 5, at 10–12 mm. These 5 values are not reliable (`cfd/README_case_selection.md`).

**Why the CT image is used, and not the labels alone (see section 3, alternative 2).** The labels do not contain the curtain or the leaflets. The dark line exists only in the CT image.

---

## 3. Alternatives that were tried

Three methods are in the project files. They were made in this sequence.

### 3.1 Fixed depths (prototype, `lvot_cas.py`, 2026-11-10): superseded
- **What it did:** short-axis slices at fixed depths of 0–12 mm below the valve plane. Run on cases 1–200.
- **What it showed (checked):** see `tracings_53.png`. At 0 mm the slice still cuts the aortic cusps. At 2 mm the outline is clean. From approximately 4–6 mm, a part of the outline is on the LA side, with no tissue wall.
- **Why a fixed depth is a weak definition (checked on the full run):**
  - The curtain end is at a different depth in each patient: 2.8–12.3 mm, middle half 6.4–8.4 mm.
  - Ellipticity increases with depth: 1.27 at 2 mm and 1.49 at the curtain end.
  - Thus one fixed depth measures a different part of the LVOT in each patient.
- **What was kept:** the 2 mm level. It is in all cases, and it is the same level as the primary measure in the UK Biobank meshes.
- **Not recorded:** the discussion that led from fixed depths to the anatomical end.

### 3.2 Label contact (border scan, 2026-11-10): superseded
- **What it did:** at each depth, it calculated the fraction of the LVOT outline that touches the LA label, the myocardium label, or unlabelled tissue. Files: `lvot_border_scan_1-200.csv`, `lvot_border_by_depth.png`.
- **What it showed (checked, 198 cases):** the outline first touches the LA label (> 5% of the outline) at approximately 3–5 mm in most cases.
- **Recorded:** this scan "led to the curtain rule" (`lvot_ellipticity_prelim_data.md`, step 1).
- **Not recorded:** why the label-contact depth was not used as the end of the LVOT.
- **Probable:** the labels have no curtain and no leaflet. Where the wall is thin, the LV label and the LA label touch each other across it. Thus label contact shows where the wall is thin, not where the hinge is. In agreement with this, label contact starts at 3–5 mm, but the CT-based rule puts the hinge at a median of 7.4 mm.

### 3.3 "Rule B" (2026-11-11): accepted and locked
- **What it is:** the rule in the README. The test table has the name `curtain_end_ruleB_29cases.csv`. The script `lvot_curtain.py` gives the same result as this table (checked for case 1).
- **"Rule A": not recorded.** The name "rule B" tells us that there was a rule A. There is no description of rule A in the project files, and no result table for it. The author must add what rule A was and why it was rejected.

### 3.4 Other landmarks
- **Not recorded.** The project files do not show that other landmarks (for example, the leaflet tip, or the level of the mitral annulus plane) were tested.

---

## 4. Why each parameter has its value

The values are in `lvot_curtain.py` under "parameters (locked)". The column "What it does" comes from the code. For most values, the reason was not written down on 2026-11-11.

| Parameter | Value | What it does | Why this value |
|---|---|---|---|
| Start depth | 2 mm below the valve plane | The trace and the search for the end start here. | **Not recorded.** Probable: the valve plane comes from the label boundary, and at 0 mm the slice still cuts the cusps (`tracings_53.png`). At 2 mm the LVOT is clear of the valve. It is also the fixed measurement level. |
| Depth step | 0.1 mm | The long-axis planes are resampled on a 0.1 mm grid. One row of the grid is one step. | **Not recorded.** Probable: a fine grid lets the trace move smoothly. It does not add information beyond the CT voxel size. |
| Number of planes | 5 planes, 2 mm apart (−4 to +4 mm) | Five parallel long-axis planes across the LVOT. Each gives its own end. | **Not recorded.** Checked support: one plane alone is not sufficient. Planes of the same case differ by a median of 1.9 mm (cases with 3 or more ends), and many planes give no end. |
| Search window | 6 mm in front of the LA edge, 1.5 mm behind it | The wall is searched for only in this band. | **Not recorded.** Probable: 6 mm in front lets the trace follow the leaflet when it moves away from the LA. 1.5 mm behind allows for an LA label edge that is 0.5–1 mm out of position (README, section 6). |
| Lateral jump | 0.3 mm for each 0.1 mm step | A point is rejected if it is farther than this from the last accepted point. | **Not recorded.** Probable: the wall is a continuous line. This stops the trace from a jump to a different dark structure. |
| Contrast | 30 HU or more | The dark point must be 30 HU or more darker than the tissue within ±1.5 mm. | **Not recorded.** Probable: it makes sure that the point is a real wall between two bright spaces, and not noise in the blood pool. |
| Stop | 1 mm of rejected points | The trace stops after 10 rejected steps in a row. | **Not recorded.** Probable: after 1 mm without a wall, the trace has lost the wall, and points found later are not reliable. |
| Gap smoothing | 1 mm | The gap is averaged over 10 steps. | **Not recorded.** Probable: it removes step-to-step noise before the thresholds are applied. |
| Gap threshold | > 0.2 mm | The first depth with a larger gap is a possible end. | **Not recorded.** Note: 0.2 mm is 2 steps of the 0.1 mm grid, the smallest gap that is more than one grid step. |
| Confirmation | Over the next 3 mm: no decrease, and a rise of more than 0.3 mm | The possible end is accepted only if the gap continues to open. In the code, "no decrease" is tested each 0.5 mm with a tolerance of 0.05 mm. | **Not recorded.** Probable: an open leaflet moves steadily away from the LA. A small local gap that closes again is label or image noise, not the hinge. |
| Combining planes | Median of 3 or more planes. With only 2 planes: the mean, if they agree within 1.5 mm. Otherwise "not found". | Gives one value for each case. | **Recorded in part:** the README says not to calculate the end by other methods when the rule fails (section 2.8). Probable: the median is not moved by one bad plane. Two planes cannot show which one is wrong, so they must agree. One plane has no check at all. |

**No sensitivity test.** The project files contain no test of these values. It is not known how much the result changes if a value changes.

---

## 5. How the rule was checked

- **29 test cases (2026-11-11, recorded).** The rule found the end in 23 of 29 (79%). Median 7.2 mm, range 3.0–10.8 mm. Table: `curtain_end_ruleB_29cases.csv`.
  - **Not recorded:** how the 29 cases were selected.
- **One hand mark (recorded).** Case 119: hand mark 5.3 mm. The rule gave "not found", because only 2 planes gave an end and they did not agree (3.5 and 7.2 mm in the test). The hand mark is between the two values.
- **Full run (2026-11-11, checked).** 998 cases, 0 errors.
  - End found in 762 cases (76%): median 7.4 mm, middle half 6.4–8.4 mm, range 2.8–12.3 mm.
  - Number of planes behind each result: 5 planes in 288 cases, 4 in 207, 3 in 191, and 2 planes that agree in 76.
  - The script gave the same result as the test table for case 1.
- **QC images (recorded).** `lvot_out/qc/{id}.png` shows the five planes, the trace and the end for each case. The README tells you to examine a sample.

---

## 6. Known weaknesses

- **236 of 998 cases have no result (24%).**
  - 96 cases: no plane gave an end.
  - 104 cases: only 1 plane gave an end.
  - 36 cases: 2 planes gave an end, but they differed by more than 1.5 mm.
- **Systole does not explain these.** Only 9 of the 236 have the systolic flag. The cause for the other 227 was not studied.
- **The cases with a result can be a selected group.** If the rule fails more frequently for some anatomy or some image quality, the 762 cases do not represent all 998. This was not tested.
- **Only one hand mark.** There is no comparison with manual marks on a set of cases, and no test of agreement between observers.
- **No systematic QC review.** The extreme values were not all examined: 16 cases are above 11 mm and 3 are below 3 mm. Case 69 (12.0 mm, 2 planes) looks doubtful in its QC image: the trace stays next to the LA edge and no clear opening is visible.
- **The result depends on the LA label edge (recorded).** A label that is 0.5–1 mm out of position can move the end by 1–2 mm when the gap opens slowly.
- **Results from 2 planes are the weakest.** 76 results come from only 2 planes.
- **Diastole only.** Clinical LVOT measurements are made in systole. This rule cannot find the end in systole.
- **Few QC images were examined for this file.** The statements about QC images here come from cases 2, 69 and 119 only.

---

## 7. What the author must add

These items are not in the project files. Only the author can supply them:
1. What rule A was, and why rule B replaced it.
2. Why the label-contact depth (section 3.2) was rejected. Confirm or correct the "Probable" reason.
3. The reason for each value in section 4 that has the label "Not recorded".
4. How the 29 test cases were selected, and how the hand mark of case 119 was made.
5. Any reference that was used to set the rule.
