"""Run provenance: every script appends one JSON line to outputs/logs/runs.jsonl.

Records time, script, arguments, SHA-256 of the code (src/*.py + the script) and config.yaml,
Python and key package versions, plus script-specific results (counts, output files).
"""
import hashlib
import json
import platform
import sys
import time
from importlib import metadata
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LOG = ROOT / "outputs/logs/runs.jsonl"
PKGS = ("numpy", "scipy", "pandas", "pyvista", "vtk", "matplotlib", "pyyaml")


def _sha(paths):
    h = hashlib.sha256()
    for p in sorted(paths):
        h.update(Path(p).read_bytes())
    return h.hexdigest()[:16]


def log_run(script: str, **extra) -> dict:
    script = Path(script).resolve()
    rec = {
        "time": time.strftime("%Y-%m-%d %H:%M:%S"),
        "script": script.relative_to(ROOT).as_posix(),
        "argv": sys.argv[1:],
        "code_sha": _sha(list((ROOT / "src").glob("*.py")) + [script]),
        "config_sha": _sha([ROOT / "config.yaml"]),
        "python": platform.python_version(),
        "packages": {p: _ver(p) for p in PKGS},
        **extra,
    }
    LOG.parent.mkdir(parents=True, exist_ok=True)
    with open(LOG, "a") as f:
        f.write(json.dumps(rec, default=str) + "\n")
    return rec


def _ver(p):
    try:
        return metadata.version(p)
    except metadata.PackageNotFoundError:
        return None
