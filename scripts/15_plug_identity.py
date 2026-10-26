"""STEP 3.0 (check) — Is elemTag 5 the aortic valve? Anatomical tests of the four valve plugs (tags 4-7).

The files label regions only by number; the paper (Ugurlu et al. 2025) lists six regions (LV, RV, aortic, mitral,
tricuspid, pulmonary) but not which number is which. For each plug this script records:
  ventricle      which myocardium its contact nodes belong to (LV tag 1 / RV tag 0)
  d_area         area-derived diameter of the contact ring (mm)
  centre         ring centre; pairwise centre distances between plugs
  angle_long     angle between the ring-plane normal and the LV long axis (deg)
Expected anatomy: aortic + mitral on the LV, tricuspid + pulmonary on the RV; aortic smaller than mitral;
aortic and pulmonary roots are neighbours; mitral and tricuspid are neighbours; pulmonary farthest from mitral.
Output: outputs/tables/plug_identity.csv
Usage: python scripts/15_plug_identity.py [--n 40]   (subset + first n/2 pilot pairs; seed from config)
"""
import argparse
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
import pyvista as pv

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src.config import load_config  # noqa: E402
from src.lvot import ring_measure  # noqa: E402
from src.provenance import log_run  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--n", type=int, default=40)
a = ap.parse_args()
cfg = load_config()
lab = pd.read_csv(ROOT / cfg["labels"])
pr = lab[lab.sex == "F"].merge(lab[lab.sex == "M"], on=["age", "bmi"], suffixes=("_F", "_M")).sample(30, random_state=cfg["seed"])
ids = list(dict.fromkeys(cfg["subset"] + [m for p in zip(pr.mesh_id_F, pr.mesh_id_M) for m in p]))[: a.n]
lab = lab.set_index("mesh_id")
NAME = {4: "mitral?", 5: "aortic?", 6: "tricuspid?", 7: "pulmonary?"}

rows = []
for k, mid in enumerate(ids, 1):
    g = pv.read(Path(lab.loc[mid, "folder"]) / f"{mid}.fiber.vtu")
    tags = np.asarray(g.cell_data["elemTag"])
    conn = g.cells.reshape(-1, 5)[:, 1:]
    P = np.asarray(g.points) * 1e-3
    lv, rv = set(np.unique(conn[tags == 1]).tolist()), set(np.unique(conn[tags == 0]).tolist())
    lvc = P[list(lv)].mean(0)
    lvn = np.array(sorted(lv))
    lvp = P[lvn]
    row = {"mesh_id": mid, "sex": lab.loc[mid, "sex"]}
    C, N = {}, {}
    for t in (4, 5, 6, 7):
        pts = np.unique(conn[tags == t])
        in_lv = np.isin(pts, lvn).sum(); in_rv = len(set(pts.tolist()) & rv)
        ring = pts[np.isin(pts, lvn)] if in_lv >= in_rv else np.array(sorted(set(pts.tolist()) & rv))
        R = ring_measure(P[ring], lvc)
        C[t], N[t] = R["centre"], R["normal"]
        row[f"t{t}_ventricle"] = "LV" if in_lv >= in_rv else "RV"
        row[f"t{t}_d_area"] = R["d_area"]
    # LV long axis: tag-4 centre to the LV point farthest from it
    apex = lvp[np.argmax(np.linalg.norm(lvp - C[4], axis=1))]
    la = (apex - C[4]) / np.linalg.norm(apex - C[4])
    for t in (4, 5, 6, 7):
        row[f"t{t}_angle_long"] = np.degrees(np.arccos(abs(N[t] @ la)))
    for t1, t2 in combinations((4, 5, 6, 7), 2):
        row[f"dist_{t1}_{t2}"] = np.linalg.norm(C[t1] - C[t2])
    rows.append(row)
    print(f"[{k}/{len(ids)}] {mid}", flush=True)
D = pd.DataFrame(rows)
D.round(2).to_csv(ROOT / "outputs/tables/plug_identity.csv", index=False)

pd.set_option("display.width", 200)
print("\nventricle each plug touches:")
for t in (4, 5, 6, 7):
    print(f"  tag {t} {NAME[t]:11s}", D[f"t{t}_ventricle"].value_counts().to_dict())
print("\nmedian (min-max):")
for c in [c for c in D.columns if c.startswith(("t4_d", "t5_d", "t6_d", "t7_d", "t4_a", "t5_a", "t6_a", "t7_a", "dist"))]:
    print(f"  {c:16s} {D[c].median():6.1f} ({D[c].min():.1f}-{D[c].max():.1f})")
chk = {
    "aortic (5) smaller than mitral (4)": (D.t5_d_area < D.t4_d_area).mean(),
    "5 nearer to pulmonary (7) than 4 is": (D.dist_5_7 < D.dist_4_7).mean(),
    "pulmonary (7) farthest plug from mitral (4)": (D.dist_4_7 > D[["dist_4_5", "dist_4_6"]].max(axis=1)).mean(),
    "mitral (4) nearer tricuspid (6) than pulmonary (7)": (D.dist_4_6 < D.dist_4_7).mean(),
}
print("\nchecks (fraction of meshes passing):")
for k_, v in chk.items():
    print(f"  {k_:50s} {v:.0%}")
log_run(__file__, n=len(D), checks={k_: round(float(v), 3) for k_, v in chk.items()})
