#!/usr/bin/env bash
# Reproduce the LVOT report (STEPS 0-9) end to end. See REPRODUCE.md (what each step makes) and docs/PIPELINE.md (methods).
# Run from this folder. The raw meshes must be at data_dir in config.yaml. About 4 h in total; the long steps
# (01 --all, 09 --all, 22 --all) are resumable. On macOS run it under caffeinate in Terminal:
#   caffeinate -i bash run_all.sh > outputs/logs/run_all.log 2>&1
# STEPS 10-12 (strain models, scripts 25-28) are not part of the report and are not run here.
set -euo pipefail
cd "$(dirname "$0")"
PY="${PY:-.venv/bin/python}"   # Python 3.11 with requirements.txt installed

# --- needs the raw meshes ---
$PY scripts/00_labels.py                              # STEP 0    labels -> data/labels.csv
$PY scripts/01_measure.py --subset                    # STEP 1-5.1 on the 8-mesh dev subset
$PY scripts/02_qc_render.py --subset                  # STEP 5.2  QC renders
$PY scripts/02_qc_render.py female_age_73_bmi_43 male_age_66_bmi_42 female_age_59_bmi_42   # STEP 5.2 outliers checked
$PY scripts/04_gap_demo.py female_age_65_bmi_26       # STEP 5.3  chord-closure demo
$PY scripts/05_outline_demo.py female_age_65_bmi_26   # STEP 5.4  outline demo
$PY scripts/08_outline_steps.py female_age_65_bmi_26  # STEP 5.5  every step from mesh to outline
$PY scripts/01_measure.py --all                       # STEP 1-5.1 all meshes (~70 min, resumable)
$PY scripts/09_save_outlines.py --all                 # STEP 5.6  2D outlines (~65 min, resumable)
$PY scripts/15_plug_identity.py                       # STEP 3.0  which plug is which (check)
$PY scripts/11_plax_definitions.py                    # STEP 6.9  PLAX plane definitions (68 meshes, ~4 min)
$PY scripts/16_mesh_explainer.py                      # STEP 7d   "what is a mesh" figure
$PY scripts/22_wall_geometry.py --n 40 --render 8     # STEP 9    wall dev subset + QC pictures
$PY scripts/22_wall_geometry.py --all                 # STEP 9    wall thickness, all meshes (~1.5 h, resumable)
$PY scripts/24_thickness_method_check.py              # STEP 9.9  thickness method check

# --- tables only (06 also needs data/outlines/, 07 re-reads one mesh) ---
$PY scripts/03_analysis.py                            # STEP 6.1-6.5 tables + figures
$PY scripts/06_area_methods.py                        # STEP 6.6  outline area vs ellipse / echo circle
$PY scripts/07_area_methods_demo.py                   # STEP 6.7  area methods illustration
$PY scripts/10_ellipticity_methods.py                 # STEP 6.8  ellipticity by method
$PY scripts/13_size_by_age_figure.py                  # STEP 6.10 diameter by age
$PY scripts/13_size_by_age_figure.py --area           # STEP 6.10 area by age
$PY scripts/18_ava_simulation.py                      # STEP 6.11 AVA circular vs elliptical LVOT
$PY scripts/14_project_flowchart.py                   # STEP 7c   study flowchart
$PY scripts/19_ellipticity_depth_slide.py             # STEP 7e   ellipticity slide figures
$PY scripts/21_shape_vs_ellipticity.py                # STEP 8.1  shape vs ellipticity (tables)
$PY scripts/23_wall_phenotypes.py                     # STEP 9.6-9.8 feature map + phenotypes
# PowerPoint files (optional, need Node.js + pptxgenjs): see REPRODUCE.md, STEP 7
echo "done; provenance in outputs/logs/runs.jsonl"
