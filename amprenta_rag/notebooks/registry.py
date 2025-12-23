"""Notebook template registry utilities."""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, List, Optional


def _default_registry_path() -> Path:
    root = Path(__file__).resolve().parents[2]  # .../RAG
    return root / "deploy" / "jupyterhub" / "templates" / "registry.json"


@lru_cache(maxsize=1)
def load_registry(registry_path: Optional[str] = None) -> List[Dict[str, Any]]:
    """Load the notebook template registry JSON."""
    p = Path(registry_path) if registry_path else _default_registry_path()
    if not p.exists():
        raise FileNotFoundError(str(p))
    data = json.loads(p.read_text(encoding="utf-8"))
    if not isinstance(data, list):
        raise ValueError("registry.json must be a list of template objects")
    return data


def search_templates(
    query: Optional[str] = None,
    tags: Optional[List[str]] = None,
    omics_type: Optional[str] = None,
    *,
    registry_path: Optional[str] = None,
) -> List[Dict[str, Any]]:
    """Search templates by free-text query and filters."""
    items = load_registry(registry_path=registry_path)
    q = (query or "").strip().lower()
    tagset = {t.strip().lower() for t in (tags or []) if t and str(t).strip()}
    om = (omics_type or "").strip().lower() if omics_type else None

    out: List[Dict[str, Any]] = []
    for it in items:
        if not isinstance(it, dict):
            continue
        if om and str(it.get("omics_type") or "").strip().lower() != om:
            continue
        if tagset:
            it_tags = {str(t).strip().lower() for t in (it.get("tags") or []) if t}
            if not tagset.issubset(it_tags):
                continue
        if q:
            blob = " ".join(
                [
                    str(it.get("id") or ""),
                    str(it.get("title") or ""),
                    str(it.get("description") or ""),
                    " ".join([str(t) for t in (it.get("tags") or [])]),
                ]
            ).lower()
            if q not in blob:
                continue
        out.append(it)
    return out


def get_template(template_id: str, *, registry_path: Optional[str] = None) -> Dict[str, Any]:
    """Get a single template by id."""
    tid = (template_id or "").strip()
    if not tid:
        raise KeyError("template id is empty")
    for it in load_registry(registry_path=registry_path):
        if isinstance(it, dict) and str(it.get("id") or "") == tid:
            return it
    raise KeyError(f"template not found: {tid}")


def resolve_notebook_path(template: Dict[str, Any], *, registry_path: Optional[str] = None) -> Path:
    """Resolve the notebook_path from a registry entry to an absolute path."""
    p = Path(registry_path) if registry_path else _default_registry_path()
    base = p.parent
    nb_rel = str(template.get("notebook_path") or "").strip()
    if not nb_rel:
        raise ValueError("template.notebook_path missing")
    nb_path = (base / nb_rel).resolve()
    return nb_path


__all__ = ["load_registry", "search_templates", "get_template", "resolve_notebook_path"]


