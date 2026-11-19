"""Load config.yaml (project root)."""
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = ROOT / "config.yaml"


def load_config() -> dict:
    with open(CONFIG_PATH) as f:
        return yaml.safe_load(f)


def mesh_folder(folder) -> Path:
    """Mesh folder from data/labels.csv (relative to data_dir) as a full path on this computer."""
    p = Path(folder)
    if p.is_absolute():
        return p
    d = Path(load_config()["data_dir"])
    return (d if d.is_absolute() else ROOT / d) / p
