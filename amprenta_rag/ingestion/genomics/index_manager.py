"""Genomics reference index manager.

Provides CRUD-style helpers to register and manage reference transcriptome indices
for Salmon/Kallisto runs.
"""

from __future__ import annotations

import os
import shutil
from pathlib import Path
from typing import List, Optional

from amprenta_rag.database.models import GenomicsIndex
from amprenta_rag.database.session import db_session


INDEX_ROOT = Path("/data/genomics/indices")


def _safe_slug(s: str) -> str:
    return (
        s.strip()
        .replace("/", "_")
        .replace("\\", "_")
        .replace(" ", "_")
    )


def _index_storage_dir(tool: str, organism: str, version: str) -> Path:
    return INDEX_ROOT / _safe_slug(tool) / f"{_safe_slug(organism)}_{_safe_slug(version)}"


def _dir_size_bytes(p: Path) -> int:
    total = 0
    for root, _dirs, files in os.walk(p):
        for fn in files:
            try:
                total += (Path(root) / fn).stat().st_size
            except Exception:
                continue
    return total


def register_index(organism: str, tool: str, version: str, file_path: str, user_id: str) -> GenomicsIndex:
    """Register a genomics index and place it under the standard storage path."""
    src = Path(file_path)
    dest_dir = _index_storage_dir(tool=tool, organism=organism, version=version)
    dest_dir.parent.mkdir(parents=True, exist_ok=True)

    stored_path: Path
    size_bytes: Optional[int] = None

    # Best-effort copy into managed storage. If it fails, keep original path.
    try:
        if src.is_dir():
            if dest_dir.exists():
                shutil.rmtree(dest_dir)
            shutil.copytree(src, dest_dir)
            stored_path = dest_dir
            size_bytes = _dir_size_bytes(dest_dir)
        else:
            dest_dir.mkdir(parents=True, exist_ok=True)
            stored_path = dest_dir / src.name
            shutil.copy2(src, stored_path)
            try:
                size_bytes = stored_path.stat().st_size
            except Exception:
                size_bytes = None
    except Exception:
        stored_path = src
        try:
            size_bytes = _dir_size_bytes(src) if src.is_dir() else src.stat().st_size
        except Exception:
            size_bytes = None

    with db_session() as db:
        obj = GenomicsIndex(
            organism=organism,
            tool=tool,
            version=version,
            file_path=str(stored_path),
            file_size_bytes=size_bytes,
            uploaded_by=user_id,
            metadata_={"original_path": file_path},
        )
        db.add(obj)
        db.flush()
        db.refresh(obj)
        return obj


def list_indices(tool: str | None = None, organism: str | None = None) -> List[GenomicsIndex]:
    """List available genomics indices (optionally filtered)."""
    with db_session() as db:
        q = db.query(GenomicsIndex)
        if tool:
            q = q.filter(GenomicsIndex.tool == tool)
        if organism:
            q = q.filter(GenomicsIndex.organism == organism)
        return q.order_by(GenomicsIndex.uploaded_at.desc()).all()


def get_index(index_id: str) -> GenomicsIndex:
    """Get a genomics index by id; raises if not found."""
    with db_session() as db:
        obj = db.query(GenomicsIndex).filter(GenomicsIndex.id == index_id).first()
        if obj is None:
            raise ValueError(f"Index not found: {index_id}")
        return obj


def delete_index(index_id: str) -> bool:
    """Delete a genomics index record (best-effort deletes the underlying files)."""
    with db_session() as db:
        obj = db.query(GenomicsIndex).filter(GenomicsIndex.id == index_id).first()
        if obj is None:
            return False

        # Best-effort filesystem cleanup.
        try:
            p = Path(obj.file_path)
            if p.is_dir():
                shutil.rmtree(p)
            elif p.exists():
                p.unlink()
        except Exception:
            pass

        db.delete(obj)
        return True


__all__ = [
    "register_index",
    "list_indices",
    "get_index",
    "delete_index",
    "INDEX_ROOT",
]


