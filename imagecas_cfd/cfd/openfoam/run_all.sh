#!/bin/bash
# Run OpenFOAM twin cases one after another in the opencfd/openfoam-run:2512 container.
# Usage (from anywhere):  cfd/openfoam/run_all.sh <tag> [case-glob]      e.g.  run_all.sh lam_h1  or  run_all.sh lam_h1 '388_*'
# MESH_ONLY=1 run_all.sh lam_h1 388_real   -> mesh only. Finished cases (solve.done) are skipped, so it can be restarted.
export PATH=/Applications/Docker.app/Contents/Resources/bin:$PATH
IMAGE=opencfd/openfoam-run:2512
TAG=${1:?tag}; PAT=${2:-*_*}
cd "$(dirname "$0")/../runs/$TAG" || exit 1
for c in $PAT; do
    [ -d "$c" ] || continue
    if [ -f "$c/solve.done" ]; then echo "skip $c (done)"; continue; fi
    echo "$(date '+%F %T') start $c"
    docker run --rm --user "$(id -u):$(id -g)" -e MESH_ONLY="$MESH_ONLY" \
        -v "$PWD/$c:/case" -w /case "$IMAGE" -c "cd /case && ./Allrun" || echo "$(date '+%F %T') FAILED $c (see $c/log.*)"
done
echo "$(date '+%F %T') all done"
