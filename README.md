# LVOT size and ellipticity in UK Biobank representative meshes

Data: 1,423 representative biventricular meshes (Ugurlu et al., PLoS One 2025; Zenodo 10.5281/zenodo.15649643, CC0) at `data/meshes`. Each mesh is the average of the participants with the same sex, age (46–83) and integer BMI. The results therefore describe typical hearts, not individuals, and CIs and p-values are descriptive only (participant counts per mesh are unknown).

## Start here
- **Reproduce every result and figure (STEPS 0–9), step by step:** [`lvot_prelim.md`](lvot_prelim.md).
- **Figures and the code for each:** [`FIGURES.md`](FIGURES.md).
- **Results and caveats:** [`WRITEUP.md`](WRITEUP.md).
- **Setup on a new machine:**
  ```bash
  python3.11 -m venv .venv && .venv/bin/pip install -r requirements.txt
  ```
  Then set `data_dir` in `config.yaml` to the unzipped meshes from Zenodo, and use `PY=.venv/bin/python`.
- **Separate validation from CT images:** code, tables and QC figures in [`imagecas_cfd/`](imagecas_cfd/README.md); refresh with `bash scripts/sync_imagecas_cfd.sh`.
- **Without the meshes:** the tables in `outputs/tables/` are included, so the analysis figures can be redrawn without the raw meshes.

## Reproduce
- **Docs:** step-by-step pipeline with the code for each step in [`docs/PIPELINE.md`](docs/PIPELINE.md); flowchart in [`docs/flowchart.mmd`](docs/flowchart.mmd) (Mermaid); reasons for each choice in [`docs/DECISIONS.md`](docs/DECISIONS.md).
- **Parameters:** everything in [`config.yaml`](config.yaml).
- **Environment:** Python 3.11 + [`requirements.txt`](requirements.txt) (this machine: `~/.venvs/lvot-ukb`).
- **Provenance:** every run appends to `outputs/logs/runs.jsonl` (code/config hash, package versions).
- **Code:** each step is marked `# STEP x.y`. Versioned with git.

```
bash run_all.sh                                 # everything, in order (see the file for each step)
PY=../ukb-mesh/.venv/bin/python
$PY scripts/00_labels.py                        # STEP 0: data/labels.csv from raw folder names
$PY scripts/01_measure.py --subset              # STEP 1-5.1 on the 8-mesh subset -> outputs/tables/lvot_subset.csv
$PY scripts/01_measure.py --all                 # all meshes -> outputs/tables/lvot_all.csv (resumable, ~3 s/mesh)
$PY scripts/02_qc_render.py --subset            # STEP 5.2 QC renders -> outputs/figures/qc/<id>.png
$PY scripts/04_gap_demo.py <mesh_id>            # STEP 5.3 chord-closing illustration
$PY scripts/05_outline_demo.py <mesh_id>        # STEP 5.4 how an outline is obtained (3D + en face)
$PY scripts/08_outline_steps.py <mesh_id>       # STEP 5.5 every step from mesh to outline (9 panels)
$PY scripts/03_analysis.py                      # STEP 6 figures + tables
$PY scripts/09_save_outlines.py --all           # STEP 5.6 save 2D outlines per mesh (resumable, ~65 min)
$PY scripts/06_area_methods.py                  # STEP 6.6 outline area vs ellipse / echo circular area
$PY scripts/07_area_methods_demo.py             # STEP 6.7 illustration of the five area estimates
$PY scripts/10_ellipticity_methods.py           # STEP 6.8 ellipticity by method: median and IQR
$PY scripts/13_size_by_age_figure.py [--area]   # STEP 6.10 diameter (or area) by age, median and IQR
$PY scripts/18_ava_simulation.py                # STEP 6.11 AVA circular vs elliptical LVOT, same VTIs
$PY scripts/14_project_flowchart.py             # STEP 7c study flowchart (PRISMA style)
NODE_PATH=<dir with pptxgenjs> node scripts/12_highlights_deck.js   # STEP 7b highlights deck -> outputs/LVOT_highlights.pptx
```
Full run (~70 min). Run it yourself in Terminal from this folder:
`caffeinate -i ../ukb-mesh/.venv/bin/python scripts/01_measure.py --all > outputs/logs/measure_all.log 2>&1`

## Results
_See [`WRITEUP.md`](WRITEUP.md)._

## Licence
Code: MIT ([`LICENSE`](LICENSE)). Source data keep their own terms: UK Biobank representative meshes CC0 (Ugurlu et al. 2025), Chan et al. 2024 data CC BY 4.0, ImageCAS under its own licence.
