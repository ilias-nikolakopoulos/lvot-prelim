# CFD twin geometries (2026-11-12)

Script: `build_twins.py`. Summary: `summarize_twins.py` → `twin_geometry_summary.csv`.
Run from the `imagecas` folder:

```bash
uv run --python 3.12 --with nibabel --with scipy --with numpy --with scikit-image --with matplotlib python -I cfd/build_twins.py --ids 388,142,576,167,33,978,586,587,253,467,146,663
uv run --python 3.12 --with pandas --with numpy --with trimesh python -I cfd/summarize_twins.py
```

About 8 s per case.

## What the twins are
- Slices are cut perpendicular to the LVOT axis (the frame from `lvot_curtain.py`), from 2 mm to 12 mm below the valve plane, 0.5 mm apart.
  - The 0 mm slice cuts the sinuses and coronary ostia, so the segment starts at 2 mm.
  - Below about 12 mm the slices open into the LV cavity.
- **real**: the LVOT lumen in each slice, smoothed 1 mm along the axis.
- **circle**: in each slice, a circle with the same area and the same centroid as the real lumen.
- Both twins have the same extensions: an inlet 2 equivalent diameters long (LV side) and an outlet 5 long (aorta side). Each extension repeats the end slice.
- Both twins have the same volume (within 0.04%) and both are watertight.
- STL units are mm. Each STL has the solids `inlet`, `outlet` and `wall`, ready for snappyHexMesh. Scale them by 0.001 for SI units.
- Valve leaflets are not modelled.

## Lumen contour and the mitral-leaflet fix
The DTU LV label runs through the open anterior mitral leaflet into the inflow tract (in diastole, from about 6 mm down).
The contour therefore comes from 360 rays cast from the lumen centre. Each ray stops at the first of:
- the label edge, or
- the first dark run inside the label (CT < myo + 0.65·(blood − myo) for ≥ 0.4 mm).

A dark run within 0.6 mm of the label edge is treated as wall blur, so the ray keeps the label edge.
A 25° circular median filter removes rays that slip through the gap at the leaflet tip.
The same method is applied to every case.

QC: `qc/twins_contours_systolic.png`, `qc/twins_contours_diastolic.png` (cyan = lumen, orange = circle, magenta = label).
- 2–10 mm: the contours follow the leaflet well.
- 10–12 mm, diastolic: small notches remain near the leaflet tip (978, 142), where the lumen is ambiguous.
  **So 5 and 7.5 mm are the primary echo depths, and 10 mm is secondary.**

## Geometric (area) part of the echo error
Echo area = π·AP²/4, where AP is the width through the centroid toward the LA. It is compared with the true slice area.

| Depth | Systolic median | Diastolic median | Pairs where systolic error is smaller |
|---|---|---|---|
| 5 mm | −24% | −27% | 2/6 |
| 7.5 mm | −26% | −34% | 5/6 |
| 10 mm | −27% | −37% | 5/6 |

The CFD will add the profile part: centreline PW-Doppler velocity vs true mean velocity, real vs circular twin.
