# ImageCAS CT and CFD arm

- **Start here:** [`lvot_ellipticity_prelim_data.md`](lvot_ellipticity_prelim_data.md) is the step-by-step guide. The plan is in [`CFD_PLAN_imagecas.md`](CFD_PLAN_imagecas.md).
- **How the end of the LVOT is defined (aorto-mitral curtain):** rule in [`README_lvot_curtain.md`](README_lvot_curtain.md), reasons in [`RATIONALE_aortomitral_curtain.md`](RATIONALE_aortomitral_curtain.md).
- **CT LVOT geometry:** `lvot_cas.py`, `lvot_curtain.py` ([`README_lvot_curtain.md`](README_lvot_curtain.md)), results in `lvot_out/` and the top-level CSVs.
- **CFD:** `cfd/` holds the case selection, the geometry builders (`build_twins.py`, `build_lv_twins.py`), the OpenFOAM case generator and post-processing (`cfd/openfoam/`, see [`README_cfd_run.md`](cfd/openfoam/README_cfd_run.md)), the result tables (`cfd_results_*.csv`, `echo_spectrum_*.csv`, `peak_displacement_*.csv`) and QC figures (`cfd/qc/`).
- **Not included:** the CT images and labels (download from Kaggle with `get_imagecas_chunk.sh`; ImageCAS licence applies), the OpenFOAM run folders (15 GB), the 3D meshes (stl/vtp/vtu) and the per-case QC renders in `lvot_out/qc/`.
- **Paths:** the scripts use the external-drive paths. To run them elsewhere, change the paths at the top of each script.
- **Refresh:** `bash scripts/sync_imagecas_cfd.sh` from the repo root, then commit.
