# LVOT measurement on ImageCAS: end of the LVOT = end of the aorto-mitral curtain

Rule locked on 2026-11-11. Script: `lvot_curtain.py`.

## 1. Data

- CT: ImageCAS (Kaggle `xiaoweixumedicalai/imagecas`, CC BY-NC 4.0). Cite Zeng et al., CMIG 2023.
- Labels: DTU whole-heart labels (Hansen et al., STACOM 2025, arXiv 2510.06090). The chamber labels come from TotalSegmentator `heartchambers_highres`. Do not train other segmentation models on these chamber labels (TotalSegmentator licence).
- Label values: 1 myocardium, 2 LA, 3 LV, 6 aorta.
- Do not use cases 90 and 141 (invalid labels). The script skips them.

## 2. Rule for the end of the aorto-mitral curtain

### 2.1 Set the reference
1. Find the boundary between the LV label and the aorta label.
2. Use this boundary as the aortic valve plane (0 mm).
3. Make the LVOT axis perpendicular to this plane.

### 2.2 Make the long-axis planes
1. Make five long-axis planes through the LVOT.
2. Each plane must contain the LVOT axis and the direction to the left atrium (LA).
3. Put the planes 2 mm apart, at −4, −2, 0, +2 and +4 mm.

### 2.3 Trace the posterior LVOT wall in each plane
1. Start 2 mm below the valve plane. Go down in steps of 0.1 mm.
2. At each step, find the front edge of the LA label.
3. Search from 6 mm in front of this edge to 1.5 mm behind it.
4. Find the darkest point in this search area. This point is the wall.
5. Above the hinge, the wall is the curtain. Below the hinge, the wall is the anterior mitral leaflet.

### 2.4 Apply the safeguards
1. Reject a point if it is more than 0.3 mm to the side of the last accepted point.
2. Reject a point if it is not 30 HU or more darker than the tissue within 1.5 mm.
3. Stop the trace after 1 mm of rejected points.

### 2.5 Measure the gap
1. At each step, measure the distance from the wall point to the LA edge. This distance is the gap.
2. Use only accepted points.
3. Smooth the gap over 1 mm.

### 2.6 Find the curtain end in one plane
1. Start at 2 mm below the valve plane.
2. Find the first step where the gap is more than 0.2 mm.
3. Make sure that the gap does not decrease in the next 3 mm.
4. Make sure that the gap increases by more than 0.3 mm in the next 3 mm.
5. If these conditions are correct, this step is the curtain end.

### 2.7 Find the curtain end for the case
1. If 3 or more planes give a curtain end, use the median.
2. If only 2 planes give a curtain end, and they agree within 1.5 mm, use the mean.
3. In all other conditions, record "not found".
4. The curtain end is the end of the LVOT.

### 2.8 Cases with no result
1. Do not calculate the curtain end by other methods.
2. Report the number of "not found" cases.
3. Mark systolic scans manually (column `systolic_flag` = True, LV/myocardium volume ratio < 0.6).

NOTE: In systole, the closed leaflet stays against the LA. The gap does not open, so the rule cannot find the curtain end.

## 3. LVOT measurements
The script cuts short-axis planes perpendicular to the LVOT axis. It measures the LV label cross-section at each level:
- `2mm`: 2 mm below the valve plane (all cases).
- `mid`: half of the curtain end depth (only if the curtain end is found).
- `curtain_end`: at the curtain end (only if found).

For each level: area (mm²), ellipticity (√ of the ratio of the second-moment eigenvalues), equivalent diameter, major and minor axis.

## 4. How to run (Claude Code on the Mac)
1. Open Terminal. Go to this folder: `cd "data/meshes/imagecas"`.
2. Start Claude Code in this folder. Tell it: "Read README_lvot_curtain.md and run the LVOT pipeline on all downloaded chunks."
3. Or run it yourself. One command does all five chunks:

```bash
caffeinate -i uv run --python 3.12 --with nibabel --with scipy --with numpy --with matplotlib \
  python -I lvot_curtain.py \
  --img-dir ImageCAS_1-200/1-200 ImageCAS_201-400/201-400 ImageCAS_401-600/401-600 ImageCAS_601-800/601-800 ImageCAS_801-1000/801-1000 \
  --lab-dir labels_DTU/ImageCAS-STACOM2025-02-10-2025/segmentations \
  --out lvot_out --workers 4
```

4. To run one chunk only, give only that folder to `--img-dir`. Keep the same `--out`.
5. The script skips cases that are already done. You can stop it and start it again.
6. Use 4 workers or fewer. Each worker uses approximately 2 GB of memory.

## 5. Outputs (in `--out`)
- `cases/{id}.json`: one result per case.
- `qc/{id}.png`: five long-axis planes per case, with the trace and the curtain end. Examine a sample of these.
- `lvot_curtain_results.csv`: all cases in one table. The script makes it again at the end of each run.

## 6. Known limits
- The curtain end depends on the LA label edge. A label that is 0.5–1 mm out of position can move the curtain end by 1–2 mm when the gap increases slowly.
- In the 29 test cases, the rule found the curtain end in 23 (79%). Median 7.2 mm, range 3.0–10.8 mm.
- Reference hand mark: case 119 = 5.3 mm (rule: not found).
