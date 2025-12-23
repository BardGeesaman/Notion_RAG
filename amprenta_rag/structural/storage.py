"""Structure file storage utilities."""

from __future__ import annotations

import hashlib
import os
from pathlib import Path
from uuid import UUID


def _base_dir() -> Path:
    return Path(os.environ.get("STRUCTURE_STORE_BASE", "data/structures"))


def _ext_for(file_type: str) -> str:
    ft = (file_type or "").lower().strip()
    if ft in ("pdb", "cif", "mmcif"):
        return "pdb" if ft == "pdb" else "cif"
    if ft in ("prepared", "prepared_pdb"):
        return "pdb"
    return ft or "bin"


def _md5(content: bytes) -> str:
    h = hashlib.md5()  # noqa: S324 - used for integrity checksum, not security
    h.update(content)
    return h.hexdigest()


def save_structure_file(structure_id: UUID, content: bytes, file_type: str) -> str:
    """Save structure file content to disk and return the file path.

    Storage layout:
      data/structures/{structure_id}/{file_type}.{ext}
    """
    if not structure_id:
        raise ValueError("structure_id required")
    if content is None:
        raise ValueError("content required")
    ft = (file_type or "").strip().lower() or "pdb"

    base = _base_dir()
    out_dir = base / str(structure_id)
    out_dir.mkdir(parents=True, exist_ok=True)

    ext = _ext_for(ft)
    out_path = out_dir / f"{ft}.{ext}"
    out_path.write_bytes(content)
    return str(out_path)


__all__ = ["save_structure_file"]


