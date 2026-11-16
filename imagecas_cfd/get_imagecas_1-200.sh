#!/bin/bash
# Download ImageCAS cases 1-200 (~18 GB) from Kaggle and unpack. CC BY-NC 4.0 - cite Zeng et al., CMIG 2023.
# Needs a Kaggle API token in ~/.kaggle/kaggle.json (Kaggle > Settings > API > Create New Token).
# Run on your Mac (not inside Claude):   bash get_imagecas_1-200.sh
set -e
cd "$(dirname "$0")"
DS=xiaoweixumedicalai/imagecas
for f in 1-200.change2zip 1-200.z01 1-200.z02 1-200.z03 1-200.z04 imageCAS_data_split.xlsx; do
  if [ -f "$f" ]; then echo "have $f"; continue; fi
  echo "downloading $f ..."
  uvx --python 3.12 --with six kaggle datasets download "$DS" -f "$f" -p . --quiet
  [ -f "$f.zip" ] && unzip -o -q "$f.zip" && rm "$f.zip"   # Kaggle sometimes wraps single files in a zip
done
mv -n 1-200.change2zip 1-200.zip 2>/dev/null || true
echo "joining split archive ..."
zip -s 0 1-200.zip --out 1-200_joined.zip
echo "unzipping ..."
unzip -n -q 1-200_joined.zip -d ImageCAS_1-200
echo "done: $(find ImageCAS_1-200 -name '*.nii.gz' | wc -l) nii.gz files"
