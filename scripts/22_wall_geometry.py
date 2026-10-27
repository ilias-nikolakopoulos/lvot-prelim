"""STEP 9 — Regional LV wall thickness, subaortic septum and aorto-septal angle per mesh. Resumable.

Uses src.lvot.measure_mesh (geometry) + src.wall.wall_geometry. One row per mesh to
outputs/tables/<out>.csv (skips meshes already there); errors to <out>_errors.csv.
--render N saves QC pictures for the first N meshes processed: outputs/figures/qc/wall/<mesh_id>.png.
Usage:
  python scripts/22_wall_geometry.py --n 40 --render 8   # dev subset (config subset + matched pairs) -> wall_subset
  python scripts/22_wall_geometry.py --all               # all 1,423 meshes -> wall_all (~1-1.5 h; caffeinate)
"""
import argparse
import csv
import sys
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src.config import load_config  # noqa: E402
from src.lvot import measure_mesh  # noqa: E402
from src.provenance import log_run  # noqa: E402
from src.wall import SECTORS, wall_geometry  # noqa: E402

ap = argparse.ArgumentParser()
g = ap.add_mutually_exclusive_group(required=True)
g.add_argument("--n", type=int, help="dev subset: config subset + random matched pairs (seed), first n")
g.add_argument("--all", action="store_true")
g.add_argument("--ids", nargs="+")
ap.add_argument("--render", type=int, default=0)
ap.add_argument("--out", default=None)
a = ap.parse_args()

cfg = load_config()
labels = pd.read_csv(ROOT / cfg["labels"])
if a.all:
    ids, out = list(labels.mesh_id), a.out or "wall_all"
elif a.ids:
    ids, out = a.ids, a.out or "wall_ids"
else:  # same selection as scripts/15_plug_identity.py
    pr = labels[labels.sex == "F"].merge(labels[labels.sex == "M"], on=["age", "bmi"], suffixes=("_F", "_M")).sample(30, random_state=cfg["seed"])
    ids = list(dict.fromkeys(cfg["subset"] + [m for p in zip(pr.mesh_id_F, pr.mesh_id_M) for m in p]))[: a.n]
    out = a.out or "wall_subset"
labels = labels.set_index("mesh_id")
out_csv, err_csv = (ROOT / "outputs/tables" / f"{out}{s}.csv" for s in ("", "_errors"))
done = set(pd.read_csv(out_csv)["mesh_id"]) if out_csv.exists() else set()
todo = [i for i in ids if i not in done]
print(f"{len(ids)} meshes, {len(done)} already done, {len(todo)} to do -> {out_csv.name}", flush=True)
QC = ROOT / "outputs/figures/qc/wall"


def render(mid, row, w):
    """STEP 9.5 — QC: endo points coloured by thickness, seen from the septum and from the side; PLAX strip in black."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    QC.mkdir(parents=True, exist_ok=True)
    X, t, k = w["X"], w["thick"], w["keep"]
    B = np.vstack([w["e_t"], w["la"], w["e_s"]])  # view coords: x = anterior, y = towards apex, z = septal
    fig, axs = plt.subplots(1, 3, figsize=(15, 5.5))
    for ax, (i, j), ttl in zip(axs[:2], ((0, 1), (2, 1)), ("viewed from the septum (x = anterior)", "viewed from anterior (x = septal)")):
        Y = (X - w["mvc"]) @ B.T
        o = np.argsort(Y[:, 3 - i - j])  # far side first
        sc = ax.scatter(Y[o, i], -Y[o, j], c=np.where(k, t, np.nan)[o], s=1.5, cmap="viridis", vmin=4, vmax=14)
        A = (w["aosa_pts"] - w["mvc"]) @ B.T
        ax.scatter(A[:, i], -A[:, j], s=3, c="k")
        ax.set_aspect("equal"); ax.set_title(ttl, fontsize=10); ax.axis("off")
    fig.colorbar(sc, ax=axs[:2], shrink=.6, label="wall thickness (mm)")
    ax = axs[2]
    M = np.array([[row[f"wt_{lv}_{s}"] for s in SECTORS] for lv in cfg["wall"]["levels"]])
    im = ax.imshow(M, cmap="viridis", vmin=4, vmax=14)
    for (r_, c_), v in np.ndenumerate(M):
        ax.text(c_, r_, f"{v:.1f}", ha="center", va="center", color="w")
    ax.set_xticks(range(4), SECTORS); ax.set_yticks(range(3), list(cfg["wall"]["levels"]))
    ax.set_title(f"subaortic septum {row['wt_subaortic_septum']:.1f} mm · AoSA {row['aosa_deg']:.0f}°", fontsize=10)
    fig.suptitle(mid)
    fig.savefig(QC / f"{mid}.png", dpi=110, bbox_inches="tight")
    plt.close(fig)


fields = list(pd.read_csv(out_csv, nrows=0).columns) if out_csv.exists() else None
t0, n_ok, n_err = time.time(), 0, 0
for kk, mid in enumerate(todo, 1):
    t = time.time()
    try:
        m = measure_mesh(labels.loc[mid, "folder"], keep_geometry=True)
        w = wall_geometry(m["_geom"])
        geo = w.pop("_wall")
        row = {c: m[c] for c in ("mesh_id", "sex", "age", "bmi")} | {c: (round(v, 4) if isinstance(v, float) else v) for c, v in w.items()}
        fields = fields or list(row)
        new = not out_csv.exists()
        with open(out_csv, "a", newline="") as f:
            wr = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
            if new:
                wr.writeheader()
            wr.writerow(row)
        if n_ok < a.render:
            render(mid, row, geo)
        n_ok += 1
        print(f"[{kk}/{len(todo)}] {mid} {time.time()-t:.1f}s subAo sept {row['wt_subaortic_septum']:.1f} "
              f"mid sept {row['wt_mid_septal']:.1f} mid lat {row['wt_mid_lateral']:.1f} AoSA {row['aosa_deg']:.0f}", flush=True)
    except Exception as e:
        new = not err_csv.exists()
        with open(err_csv, "a", newline="") as f:
            wr = csv.writer(f)
            if new:
                wr.writerow(["mesh_id", "error"])
            wr.writerow([mid, repr(e)])
        n_err += 1
        print(f"[{kk}/{len(todo)}] {mid} ERROR {e!r}", flush=True)
        traceback.print_exc()
print(f"finished in {(time.time()-t0)/60:.1f} min", flush=True)
log_run(__file__, output=str(out_csv.relative_to(ROOT)), n_requested=len(ids), n_ok=n_ok, n_errors=n_err,
        minutes=round((time.time() - t0) / 60, 1))
