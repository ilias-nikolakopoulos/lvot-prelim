# LVOT size and ellipticity in UK Biobank representative meshes

Data: 1,423 representative biventricular meshes (Ugurlu et al., PLoS One 2025; Zenodo 10.5281/zenodo.15649643, CC0). Each mesh is the average of the participants with the same sex, age (46–83) and integer BMI. The results therefore describe typical hearts, not individuals, and CIs and p-values are descriptive only (participant counts per mesh are unknown).

## Start here
- **Reproduce every result and figure (STEPS 0–9), step by step:** [`lvot_prelim.md`](lvot_prelim.md).
- **Figures and the code for each:** [`FIGURES.md`](FIGURES.md).
- **Results and caveats:** [`WRITEUP.md`](WRITEUP.md).
- **Setup on a new machine:**
  ```bash
  python3.11 -m venv .venv && .venv/bin/pip install -r requirements.txt
  ```
  Then set `data_dir` in `config.yaml` to the unzipped meshes from Zenodo, and use `PY=.venv/bin/python`.
- **Separate validation from CT images (still in progress):** code, tables and QC figures in [`imagecas_cfd/`](imagecas_cfd/README.md); refresh with `bash scripts/sync_imagecas_cfd.sh`.

## Reproduce
- **Docs:** step-by-step pipeline with the code for each step in [`docs/PIPELINE.md`](docs/PIPELINE.md); flowchart in [`docs/flowchart.mmd`](docs/flowchart.mmd) (Mermaid); reasons for each choice in [`docs/DECISIONS.md`](docs/DECISIONS.md).
- **Parameters:** everything in [`config.yaml`](config.yaml).
- **Environment:** Python 3.11 + [`requirements.txt`](requirements.txt).
- **Provenance:** every run appends to `outputs/logs/runs.jsonl` (code/config hash, package versions).
- **Code:** each step is marked `# STEP x.y`.

**Without the meshes (tables only, seconds each).** The per-mesh tables in `outputs/tables/` are in the repository, so these redraw every analysis result and figure:
```bash
PY=.venv/bin/python
$PY scripts/03_analysis.py                    # STEP 6.1-6.5 tables + figures
$PY scripts/10_ellipticity_methods.py         # STEP 6.8 ellipticity by method
$PY scripts/13_size_by_age_figure.py [--area] # STEP 6.10 diameter (or area) by age
$PY scripts/18_ava_simulation.py              # STEP 6.11 AVA, circular vs elliptical LVOT
$PY scripts/14_project_flowchart.py           # STEP 7c study flowchart
$PY scripts/19_ellipticity_depth_slide.py     # STEP 7e ellipticity slide figures
$PY scripts/21_shape_vs_ellipticity.py        # STEP 8.1 shape vs ellipticity
$PY scripts/23_wall_phenotypes.py             # STEP 9.6-9.8 wall feature map + phenotypes
```

**With the meshes (full pipeline, about 4 h).** Download and unzip the meshes from Zenodo ([10.5281/zenodo.15649643](https://doi.org/10.5281/zenodo.15649643)), set `data_dir` in `config.yaml` to that folder, then run everything in order:
```bash
bash run_all.sh
```

## Results
_See [`WRITEUP.md`](WRITEUP.md)._

## Licence
Code: MIT ([`LICENSE`](LICENSE)). Source data keep their own terms: UK Biobank representative meshes CC0 (Ugurlu et al. 2025), ImageCAS under its own licence.
