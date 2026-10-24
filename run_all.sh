#!/usr/bin/env bash
# Reproduce the LVOT pipeline end to end (see docs/PIPELINE.md). Run from this folder.
# The external SSD must be mounted. STEP 1-5.1 on all meshes takes ~70 min and is resumable:
# run it under caffeinate in Terminal, e.g.
#   caffeinate -i bash run_all.sh > outputs/logs/run_all.log 2>&1
set -euo pipefail
cd "$(dirname "$0")"
PY="${PY:-../ukb-mesh/.venv/bin/python}"   # Python 3.11 with requirements.txt installed

$PY scripts/00_labels.py                      # STEP 0   labels -> data/labels.csv
$PY scripts/01_measure.py --subset            # STEP 1-5.1 on the 8-mesh dev subset
$PY scripts/02_qc_render.py --subset          # STEP 5.2 QC renders
$PY scripts/04_gap_demo.py female_age_65_bmi_26   # STEP 5.3 gap-bridging demo
$PY scripts/05_outline_demo.py female_age_65_bmi_26   # STEP 5.4 outline (plane cut -> polygon) demo
$PY scripts/08_outline_steps.py female_age_65_bmi_26  # STEP 5.5 every step from mesh to outline
$PY scripts/01_measure.py --all               # STEP 1-5.1 all meshes (resumable, ~70 min)
$PY scripts/09_save_outlines.py --all         # STEP 5.6 save 2D outlines (resumable, ~65 min)
$PY scripts/03_analysis.py                    # STEP 6   tables + figures
$PY scripts/06_area_methods.py                # STEP 6.6 outline area vs ellipse / echo circle
$PY scripts/07_area_methods_demo.py           # STEP 6.7 illustration of the five area estimates
echo "done; provenance in outputs/logs/runs.jsonl"
