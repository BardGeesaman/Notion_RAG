"""DVC manager utilities for dataset versioning.

This module provides thin wrappers around the `dvc` CLI for:
- initializing a DVC workspace
- adding/versioning files
- restoring files from cache
- parsing .dvc metadata
"""

from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Any, Dict, Optional

import yaml


def _repo_root() -> Path:
    # amprenta_rag/ingestion/dvc_manager.py -> amprenta_rag -> repo root
    return Path(__file__).resolve().parents[2]


def check_dvc_installed() -> Dict[str, Any]:
    """Check if DVC is installed.

    Returns: {"installed": bool, "version": str}
    """
    try:
        proc = subprocess.run(
            ["dvc", "--version"],
            capture_output=True,
            text=True,
            check=False,
        )
        if proc.returncode == 0:
            return {"installed": True, "version": (proc.stdout or "").strip()}
        return {"installed": False, "version": None}
    except FileNotFoundError:
        return {"installed": False, "version": None}


def init_dvc_workspace() -> None:
    """Initialize DVC in the repo if not already initialized."""
    root = _repo_root()
    if (root / ".dvc").exists() or (root / ".dvc").is_dir():
        return

    info = check_dvc_installed()
    if not info.get("installed"):
        return

    # Best-effort init; do not fail ingestion if DVC init fails.
    subprocess.run(["dvc", "init", "-q"], cwd=str(root), capture_output=True, text=True, check=False)


def _dvc_file_for(path: Path) -> Path:
    # DVC creates "<file>.dvc" for file targets.
    return Path(str(path) + ".dvc")


def get_version_info(dvc_file: str | Path) -> Dict[str, Any]:
    """Parse a .dvc YAML file into a dict."""
    p = Path(dvc_file)
    raw = p.read_text(encoding="utf-8", errors="replace")
    data = yaml.safe_load(raw) or {}
    if isinstance(data, dict):
        data["_raw_yaml"] = raw
    return data if isinstance(data, dict) else {"_raw_yaml": raw}


def version_file(file_path: str | Path) -> Optional[str]:
    """Run `dvc add` on a file and return its MD5 hash (from the .dvc file)."""
    init_dvc_workspace()
    info = check_dvc_installed()
    if not info.get("installed"):
        return None

    root = _repo_root()
    fp = Path(file_path)
    # Run from repo root so .dvc stays consistent.
    proc = subprocess.run(
        ["dvc", "add", str(fp)],
        cwd=str(root),
        capture_output=True,
        text=True,
        check=False,
    )
    if proc.returncode != 0:
        return None

    dvc_file = _dvc_file_for(fp)
    if not dvc_file.exists():
        return None

    meta = get_version_info(dvc_file)
    outs = meta.get("outs") if isinstance(meta, dict) else None
    if isinstance(outs, list) and outs:
        md5 = outs[0].get("md5") if isinstance(outs[0], dict) else None
        return str(md5) if md5 else None
    return None


def fetch_version(dvc_file: str | Path) -> bool:
    """Restore a version from DVC cache using `dvc checkout <dvc_file>`."""
    info = check_dvc_installed()
    if not info.get("installed"):
        return False

    root = _repo_root()
    p = Path(dvc_file)
    proc = subprocess.run(
        ["dvc", "checkout", str(p)],
        cwd=str(root),
        capture_output=True,
        text=True,
        check=False,
    )
    return proc.returncode == 0


def _write_dvc_file(dvc_file: Path, raw_yaml: str) -> None:
    dvc_file.parent.mkdir(parents=True, exist_ok=True)
    dvc_file.write_text(raw_yaml, encoding="utf-8")


def restore_from_metadata(dvc_file: str | Path, raw_yaml: str) -> bool:
    """Restore a specific version by writing the stored .dvc YAML then running checkout."""
    try:
        _write_dvc_file(Path(dvc_file), raw_yaml)
    except Exception:
        return False
    return fetch_version(dvc_file)


__all__ = [
    "check_dvc_installed",
    "init_dvc_workspace",
    "version_file",
    "fetch_version",
    "get_version_info",
    "restore_from_metadata",
]


