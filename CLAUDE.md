# CLAUDE.md — LVOT side project
See README.md for methods, docs/PIPELINE.md for the step-by-step pipeline (keep it, docs/flowchart.mmd and docs/DECISIONS.md updated whenever a step or choice changes; mark code with `# STEP x.y`). Parameters live in config.yaml, not in code. New scripts call src.provenance.log_run(). Rules: never modify ../ukb-mesh/ or raw data (data/meshes, must be mounted).
Python: ../ukb-mesh/.venv/bin/python. One mesh in memory at a time (MacBook Air).
Long runs (>10 min): give the user the `caffeinate -i ... > outputs/logs/<name>.log 2>&1` command; do not run it yourself. 01_measure.py is resumable (skips mesh_ids already in the output CSV).
Status (2026-10-23): method built and checked on the 8-mesh subset plus a 60-mesh pilot (30 random common-support pairs, already in outputs/tables/lvot_all.csv). Full run pending, then run 03_analysis.py and write WRITEUP.md.
