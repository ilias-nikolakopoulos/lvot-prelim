#!/bin/bash
# Sample U in the emulated PW gate volumes (system/sampleEchoGate, written by echo_spectrum.py write) for every finished
# case of the given run tags, one container per case (~4 s each). Skips cases already sampled.
# Usage (from the imagecas folder):  bash cfd/openfoam/echo_spectrum_sample.sh lam_h1 lv_lam_h1
export PATH=/Applications/Docker.app/Contents/Resources/bin:$PATH
cd "$(dirname "$0")/../runs" || exit 1
for tag in "$@"; do
  for c in "$tag"/*_*; do
    [ -f "$c/system/sampleEchoGate" ] || continue
    if ls "$c"/postProcessing/echoGate/*/g7p5_U.xy >/dev/null 2>&1; then echo "skip $c"; continue; fi
    docker run --rm --user "$(id -u):$(id -g)" -v "$PWD/$c:/case" -w /case opencfd/openfoam-run:2512 \
      -c 'cd /case && postProcess -dict system/sampleEchoGate -latestTime -fields "(U)" > log.sampleEchoGate 2>&1' \
      && echo "done $c" || echo "FAILED $c"
  done
done
