"""STEPS 1-5.1 — Measure aortic annulus, LVOT sections and mitral ring for each mesh. Resumable.

Writes one row per mesh to outputs/tables/<out>.csv as it goes (skips meshes already there);
failures go to outputs/tables/<out>_errors.csv. One mesh in memory at a time.
Usage:
  python scripts/01_measure.py --subset            # 8-mesh subset from config.yaml
  python scripts/01_measure.py --all               # all 1,423 meshes -> lvot_all.csv
  python scripts/01_measure.py --ids <id> [...]
"""
import argparse
import csv
import sys
import time
import traceback
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src.config import load_config  # noqa: E402
from src.lvot import measure_mesh  # noqa: E402
from src.provenance import log_run  # noqa: E402

ap = argparse.ArgumentParser()
g = ap.add_mutually_exclusive_group(required=True)
g.add_argument("--subset", action="store_true")
g.add_argument("--all", action="store_true")
g.add_argument("--ids", nargs="+")
g.add_argument("--sample-pairs", type=int, help="N random common-support F/M pairs (seed 2026), written to lvot_all")
ap.add_argument("--out", default=None)
a = ap.parse_args()

cfg = load_config()
labels = pd.read_csv(ROOT / cfg["labels"]).set_index("mesh_id")
if a.subset:
    ids, out = cfg["subset"], a.out or "lvot_subset"
elif a.all:
    ids, out = list(labels.index), a.out or "lvot_all"
elif a.sample_pairs:
    lab = labels.reset_index()
    pr = lab[lab.sex == "F"].merge(lab[lab.sex == "M"], on=["age", "bmi"], suffixes=("_F", "_M"))
    pr = pr.sample(a.sample_pairs, random_state=cfg["seed"])
    ids, out = list(pr.mesh_id_F) + list(pr.mesh_id_M), a.out or "lvot_all"
else:
    ids, out = a.ids, a.out or "lvot_ids"

out_csv = ROOT / "outputs/tables" / f"{out}.csv"
err_csv = ROOT / "outputs/tables" / f"{out}_errors.csv"
done = set(pd.read_csv(out_csv)["mesh_id"]) if out_csv.exists() else set()
todo = [i for i in ids if i not in done]
print(f"{len(ids)} meshes, {len(done)} already done, {len(todo)} to do -> {out_csv.name}", flush=True)

fields = None
if out_csv.exists():
    fields = list(pd.read_csv(out_csv, nrows=0).columns)
t0 = time.time()
n_ok = n_err = 0
for k, mid in enumerate(todo, 1):
    t = time.time()
    try:
        row = measure_mesh(labels.loc[mid, "folder"])
        row = {c: (round(v, 4) if isinstance(v, float) else v) for c, v in row.items()}
        new = not out_csv.exists()
        if fields is None:
            fields = list(row)
        with open(out_csv, "a", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
            if new:
                w.writeheader()
            w.writerow(row)
        n_ok += 1
        print(f"[{k}/{len(todo)}] {mid} {time.time()-t:.1f}s ao {row['ao_d_area']:.1f} "
              f"mv {row['mv_d_area']:.1f} lvot_len {row['lvot_tube_len']} LVcav {row['lv_cavity_ml']:.0f} "
              f"flags={row['qc_flags'] or '-'}", flush=True)
    except Exception as e:
        new = not err_csv.exists()
        with open(err_csv, "a", newline="") as f:
            w = csv.writer(f)
            if new:
                w.writerow(["mesh_id", "error"])
            w.writerow([mid, repr(e)])
        n_err += 1
        print(f"[{k}/{len(todo)}] {mid} ERROR {e!r}", flush=True)
        traceback.print_exc()
print(f"finished in {(time.time()-t0)/60:.1f} min", flush=True)
log_run(__file__, output=str(out_csv.relative_to(ROOT)), n_requested=len(ids), n_already_done=len(done),
        n_ok=n_ok, n_errors=n_err, minutes=round((time.time() - t0) / 60, 1))
