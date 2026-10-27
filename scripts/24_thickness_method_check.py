"""STEP 9.9 (check) — Wall thickness: nearest-surface distance (used) vs perpendicular ray (CMR-like).

For each LV endocardial point, the ray starts at the point and runs along the endocardial normal into the
myocardium; thickness = distance to the first hit on the epicardium or RV cavity surface (max 40 mm).
Segment medians by both methods and their ratio, plus the mid septal / lateral ratio by both methods.
Output: outputs/tables/thickness_method_check.csv
Usage: python scripts/24_thickness_method_check.py   (the 8-mesh config subset; needs raw meshes)
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyvista as pv
import vtk
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src.config import load_config  # noqa: E402
from src.lvot import measure_mesh  # noqa: E402
from src.provenance import log_run  # noqa: E402
from src.wall import SECTORS, wall_geometry  # noqa: E402

cfg = load_config()
lab = pd.read_csv(ROOT / cfg["labels"]).set_index("mesh_id")
rows = []
for mid in cfg["subset"]:
    m = measure_mesh(lab.loc[mid, "folder"], keep_geometry=True)
    G = m["_geom"]
    w = wall_geometry(G)["_wall"]
    # outward cavity normals (into the myocardium), mapped to the endo points
    cav = G["cavity"].compute_normals(point_normals=True, cell_normals=False, auto_orient_normals=True,
                                      consistent_normals=True)
    nrm = np.asarray(cav.point_data["Normals"])[cKDTree(np.asarray(cav.points)).query(w["X"])[1]]
    S = pv.merge([w["epi"], w["rvc"]]).extract_surface()
    loc = vtk.vtkCellLocator()  # first hit along the ray (vtkOBBTree returned no hits here)
    loc.SetDataSet(S)
    loc.BuildLocator()
    t, x, pc, sid, cid, cell = vtk.reference(0.), [0.] * 3, [0.] * 3, vtk.reference(0), vtk.reference(0), vtk.vtkGenericCell()
    ray = np.full(len(w["X"]), np.nan)
    for i, (p, n) in enumerate(zip(w["X"], nrm)):
        if w["keep"][i] and loc.IntersectWithLine(p + .05 * n, p + 40 * n, 1e-6, t, x, pc, sid, cid, cell):
            ray[i] = np.linalg.norm(np.array(x) - p)
    near = w["thick"]
    r = {"mesh_id": mid, "ray_miss_frac": float(np.isnan(ray[w["keep"]]).mean())}
    for lv, (lo, hi) in cfg["wall"]["levels"].items():
        for sc in SECTORS:
            s = w["keep"] & (w["h"] >= lo) & (w["h"] < hi) & (w["sector"] == sc)
            r[f"near_{lv}_{sc}"] = np.median(near[s])
            r[f"ray_{lv}_{sc}"] = np.nanmedian(ray[s])
    for k in ("near", "ray"):
        r[f"{k}_sep_lat_mid"] = r[f"{k}_mid_septal"] / r[f"{k}_mid_lateral"]
    rows.append(r)
    print(f"{mid}: mid septal {r['near_mid_septal']:.1f} / {r['ray_mid_septal']:.1f}  mid lateral "
          f"{r['near_mid_lateral']:.1f} / {r['ray_mid_lateral']:.1f}  sep/lat {r['near_sep_lat_mid']:.2f} / "
          f"{r['ray_sep_lat_mid']:.2f}  (nearest / ray)  ray misses {r['ray_miss_frac']:.1%}", flush=True)
D = pd.DataFrame(rows)
D.to_csv(ROOT / "outputs/tables/thickness_method_check.csv", index=False)
print("\nmedian ray / nearest by segment:")
for lv in cfg["wall"]["levels"]:
    print(f"  {lv:7s}" + "  ".join(f"{sc} {np.median(D[f'ray_{lv}_{sc}'] / D[f'near_{lv}_{sc}']):.2f}" for sc in SECTORS))
print(f"mid septal/lateral: nearest {D.near_sep_lat_mid.median():.2f}, ray {D.ray_sep_lat_mid.median():.2f}")
log_run(__file__, output="outputs/tables/thickness_method_check.csv", n=len(D))
