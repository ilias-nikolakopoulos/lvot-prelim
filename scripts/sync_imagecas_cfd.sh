#!/usr/bin/env bash
# Copy the ImageCAS / CFD code, docs, result tables and QC figures from the external drive into imagecas_cfd/.
# Read-only on the source. Leaves out the CT images and labels, OpenFOAM run folders, meshes (stl/vtp/vtu) and
# per-case QC renders, per-case JSON caches and the unused annulus hand-mark images. Rerun to refresh the snapshot, then commit.
set -euo pipefail
cd "$(dirname "$0")/.."
SRC="${SRC:-data/meshes/imagecas}"   # data/meshes: local link to the mesh folder (see config.yaml)
rsync -a --delete --prune-empty-dirs --exclude '/README.md' --exclude '/RATIONALE_aortomitral_curtain.md' \
  --exclude '__pycache__/' --exclude '*.pyc' --exclude '*.bak' --exclude '.DS_Store' \
  --exclude '/ImageCAS_*/' --exclude '/labels_DTU/' --exclude '/cfd/runs/' --exclude '/lvot_out/qc/' --exclude '/annulus_out/cases/' --exclude '/sinus_2mm_out/cases/' --exclude '/annulus_handmarks/' \
  --exclude '*.tar' --exclude '*.zip' --exclude '*.log' \
  --include '*/' \
  --include '*.py' --include '*.sh' --include '*.md' --include '*.csv' --include '*.json' --include '*.png' --include '*.xlsx' \
  --exclude '*' \
  "$SRC/" imagecas_cfd/
du -sh imagecas_cfd
