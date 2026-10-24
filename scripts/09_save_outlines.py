"""STEP 5.6 — Save each mesh's 2D outlines so new diameter/area methods need no re-read of the meshes.

Per mesh -> data/outlines/<mesh_id>.npz (resumable: skips meshes already saved):
  ao_poly, mv_poly          ring outlines (72-bin polygons) in their own ring planes, mm
  ao_u                      unit in-plane PLAX direction for the annulus (towards the mitral centre)
  u_plax                    unit in-plane PLAX direction for the LVOT section planes
  lvot{d}_xy                strict closed outline at depth d (absent if not a closed tube)
  lvot{d}cc_xy              outline used for lvot{d}cc_* (strict if closed, else chord-closed; absent if rejected)
for d in the report depths. All outlines are ordered polygons (no repeated end point).
Usage: python scripts/09_save_outlines.py --all | --subset | --ids <id> ... | --from-table lvot_all
"""
import argparse
import sys
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src.config import load_config  # noqa: E402
from src.lvot import REPORT_DEPTHS, measure_mesh  # noqa: E402
from src.provenance import log_run  # noqa: E402

ap = argparse.ArgumentParser()
g = ap.add_mutually_exclusive_group(required=True)
g.add_argument("--all", action="store_true")
g.add_argument("--subset", action="store_true")
g.add_argument("--ids", nargs="+")
g.add_argument("--from-table", help="mesh ids from outputs/tables/<name>.csv")
a = ap.parse_args()
cfg = load_config()
labels = pd.read_csv(ROOT / cfg["labels"]).set_index("mesh_id")
if a.all:
    ids = list(labels.index)
elif a.subset:
    ids = cfg["subset"]
elif a.from_table:
    ids = list(pd.read_csv(ROOT / "outputs/tables" / f"{a.from_table}.csv").mesh_id)
else:
    ids = a.ids
odir = ROOT / "data/outlines"
odir.mkdir(parents=True, exist_ok=True)
todo = [i for i in ids if not (odir / f"{i}.npz").exists()]
print(f"{len(ids)} meshes, {len(ids) - len(todo)} already saved, {len(todo)} to do -> {odir.relative_to(ROOT)}", flush=True)

t0 = time.time()
n_ok = n_err = 0
for k, mid in enumerate(todo, 1):
    t = time.time()
    try:
        r = measure_mesh(labels.loc[mid, "folder"], keep_geometry=True)
        G = r["_geom"]
        ao, mv = G["rings"]["ao"], G["rings"]["mv"]
        n, e1, e2, org = ao["normal"], ao["e1"], ao["e2"], G["origin"]
        md = mv["centre"] - org
        u_plax = np.array([md @ e1, md @ e2]); u_plax /= np.linalg.norm(u_plax)
        dm = mv["centre"] - ao["centre"]
        ao_u = np.array([dm @ ao["e1"], dm @ ao["e2"]]); ao_u /= np.linalg.norm(ao_u)
        arr = {"ao_poly": ao["poly"], "mv_poly": mv["poly"], "ao_u": ao_u, "u_plax": u_plax}
        for s in G["sections"]:
            if s["depth"] not in REPORT_DEPTHS:
                continue
            dd = f"{s['depth']:g}".replace(".", "p")
            if s.get("kept"):
                arr[f"lvot{dd}_xy"] = s["xy"]
                arr[f"lvot{dd}cc_xy"] = s["xy"]
            elif s.get("cc"):
                arr[f"lvot{dd}cc_xy"] = s["cc_xy"]
        np.savez_compressed(odir / f"{mid}.npz", **arr)
        n_ok += 1
        print(f"[{k}/{len(todo)}] {mid} {time.time() - t:.1f}s {len(arr)} arrays", flush=True)
    except Exception as e:
        n_err += 1
        print(f"[{k}/{len(todo)}] {mid} ERROR {e!r}", flush=True)
        traceback.print_exc()
print(f"finished in {(time.time() - t0) / 60:.1f} min", flush=True)
log_run(__file__, n_requested=len(ids), n_ok=n_ok, n_errors=n_err, output_dir=str(odir.relative_to(ROOT)),
        minutes=round((time.time() - t0) / 60, 1))
