# CFD pilot: case selection (2026-11-12)

Question: elliptical vs circular LVOT of the same area. Primary endpoint: the echo stroke-volume error
(π(AP/2)² × centreline PW-Doppler velocity vs true flow), split into an area part and a profile part.
Design: patient twins only (real LVOT vs circular twin: same centreline, same area in every slice).
Steady flow at peak systole, same Q in both twins, valve leaflets not modelled (open valve).
The audience is an internal progress talk.

## Cases: `case_selection.csv`
Six systolic-flagged scans (LV/myocardium < 0.6), each paired with a diastolic scan
(ratio 0.9–1.3, curtain end found in 4–5 planes) matched on ellipticity and area at 2 mm.

| Pair | Systolic | Diastolic | Ellipticity 2 mm | Area 2 mm (mm²) |
|---|---|---|---|---|
| 1 | 388 | 142 | 1.18 | ~495 |
| 2 | 576 | 167 | 1.23 | ~305 |
| 3 | 33 | 978 | 1.24 | ~338 |
| 4 | 586 | 587 | 1.27 | ~295 |
| 5 | 253 | 467 | 1.33 | ~373 |
| 6 | 146 | 663 | 1.37 | ~452 |

Systolic reserves: 13, 166, 370, 825.

## Excluded after QC (`qc/`)
- 364: stair-step (band) artifact through the LVOT at +5 to +12 mm.
- 592: blurry, low contrast; posterior-wall trace fails.
- 539: small LV–aorta interface and an irregular aortic-root label.
- 808: irregular long-axis frame; LA and root shape near the valve look wrong.
- Diastolic 479 (stair-step artifact) and 29 (label leak at +12 mm) were replaced by 167 and 978.

## Caveats for the talk
- The systolic-flagged scans are **end-systolic**: at −3 mm the aortic root shows leaflet tissue in the
  aorta label, so the valve is closed. They are closer to ejection than diastolic scans, but not peak ejection.
- Curtain ends found in systolic scans (10–12 mm) are not reliable: the closed mitral leaflet lies against the LA.
  The CFD domain uses a fixed segment and does not depend on the curtain end.
- Systolic ellipticity at 2 mm covers only 1.16–1.37, so this pilot is not a wide dose-response curve.
