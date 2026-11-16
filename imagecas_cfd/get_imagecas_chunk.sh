#!/bin/bash
# Download one ImageCAS chunk (~18 GB), e.g. 201-400 from Kaggle and unpack. CC BY-NC 4.0 - cite Zeng et al., CMIG 2023.
# Needs a Kaggle API token in ~/.kaggle/kaggle.json (Kaggle > Settings > API > Create New Token).
# Run on your Mac (not inside Claude):   bash get_imagecas_$CH.sh
set -e
CH=${1:?usage: bash get_imagecas_chunk.sh 201-400}
cd "$(dirname "$0")"
DS=xiaoweixumedicalai/imagecas
for f in $CH.change2zip $CH.z01 $CH.z02 $CH.z03 $CH.z04; do
  if [ -f "$f" ]; then echo "have $f"; continue; fi
  echo "downloading $f ..."
  uvx --python 3.12 --with six kaggle datasets download "$DS" -f "$f" -p . --quiet
  [ -f "$f.zip" ] && unzip -o -q "$f.zip" && rm "$f.zip"   # Kaggle sometimes wraps single files in a zip
done
mv -n $CH.change2zip $CH.zip 2>/dev/null || true
echo "joining split archive ..."
zip -s 0 $CH.zip --out ${CH}_joined.zip
echo "unzipping ..."
unzip -n -q ${CH}_joined.zip -d ImageCAS_$CH
echo "done: $(find ImageCAS_$CH -name '*.nii.gz' | wc -l) nii.gz files"
rm -f $CH.z0? $CH.zip ${CH}_joined.zip  # free space once unpacked
