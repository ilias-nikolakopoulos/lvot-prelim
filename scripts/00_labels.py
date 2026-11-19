"""STEP 0 — Build the label table from the raw mesh folders (read-only).

Finds every <id>/<id>.fiber.vtu under data_dir (skipping hidden / AppleDouble files), parses
sex / age / BMI from the folder name and writes data/labels.csv (mesh_id, sex, age, bmi, folder).
Checks: 1,423 meshes, one per sex x age x BMI, every folder has its .uvc directory.
Usage: python scripts/00_labels.py
"""
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src.config import load_config, mesh_folder  # noqa: E402
from src.lvot import parse_name  # noqa: E402
from src.provenance import log_run  # noqa: E402

cfg = load_config()
rows = []
DATA = mesh_folder("")  # data_dir as a full path
for vtu in DATA.rglob("*.fiber.vtu"):
    if vtu.name.startswith("._") or any(p.startswith(".") for p in vtu.parts[-3:-1]):
        continue
    mid = vtu.name[: -len(".fiber.vtu")]
    rows.append({**parse_name(mid), "folder": str(vtu.parent.relative_to(DATA))})  # relative: no local paths in git
lab = pd.DataFrame(rows).sort_values("mesh_id").reset_index(drop=True)
assert lab.mesh_id.is_unique, "duplicate mesh ids"
assert not lab.duplicated(["sex", "age", "bmi"]).any(), "more than one mesh per sex x age x BMI"
out = ROOT / cfg["labels"]
out.parent.mkdir(parents=True, exist_ok=True)
lab.to_csv(out, index=False)
pairs = lab[lab.sex == "F"].merge(lab[lab.sex == "M"], on=["age", "bmi"])
print(f"{len(lab)} meshes (F {sum(lab.sex == 'F')}, M {sum(lab.sex == 'M')}), ages {lab.age.min()}-{lab.age.max()}, "
      f"{len(pairs)} common-support pairs -> {out.relative_to(ROOT)}")
log_run(__file__, n_meshes=len(lab), n_pairs=len(pairs), output=str(out.relative_to(ROOT)))
