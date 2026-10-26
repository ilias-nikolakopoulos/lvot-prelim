# CLAUDE.md — LVOT side project
See README.md for methods, docs/PIPELINE.md for the step-by-step pipeline (keep it, docs/flowchart.mmd and docs/DECISIONS.md updated whenever a step or choice changes; mark code with `# STEP x.y`). Parameters live in config.yaml, not in code. New scripts call src.provenance.log_run(). Rules: never modify ../ukb-mesh/ or raw data (data/meshes, must be mounted).
Python: ~/.venvs/lvot-ukb/bin/python (outside iCloud; `uv venv --python 3.11 ~/.venvs/lvot-ukb && uv pip install -r requirements.txt`). ../ukb-mesh/.venv gets iCloud-evicted and hangs on import. One mesh in memory at a time (MacBook Air).
Long runs (>10 min): give the user the `caffeinate -i ... > outputs/logs/<name>.log 2>&1` command; do not run it yourself. 01_measure.py is resumable (skips mesh_ids already in the output CSV).
Phase: end-diastole only (no systolic meshes exist; user decision 2026-10-24) — keep the caveat at the top of WRITEUP.md.
Status (2026-10-24): full run done (1,423 meshes, 0 errors); analysis, area-method comparison and WRITEUP.md complete.
