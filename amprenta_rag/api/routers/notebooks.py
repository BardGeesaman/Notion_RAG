"""Notebook template registry API endpoints."""

from __future__ import annotations

import os
from typing import Any, Dict, List, Optional

from fastapi import APIRouter, HTTPException, Query
from fastapi.responses import FileResponse

from amprenta_rag.notebooks.registry import get_template, load_registry, resolve_notebook_path, search_templates


router = APIRouter(prefix="/notebooks", tags=["Notebooks"])


@router.get("/templates", response_model=List[Dict[str, Any]])
def list_templates(
    query: Optional[str] = None,
    tags: Optional[List[str]] = Query(None, description="Filter by tags (repeatable)"),
    omics_type: Optional[str] = None,
) -> List[Dict[str, Any]]:
    """List notebook templates with optional search/filter."""
    return search_templates(query=query, tags=tags, omics_type=omics_type)


@router.get("/templates/{template_id}", response_model=Dict[str, Any])
def get_template_detail(template_id: str) -> Dict[str, Any]:
    try:
        return get_template(template_id)
    except KeyError as e:
        raise HTTPException(status_code=404, detail=str(e))


@router.get("/templates/{template_id}/download")
def download_template(template_id: str) -> FileResponse:
    """Return the raw .ipynb file for a template."""
    try:
        tpl = get_template(template_id)
    except KeyError as e:
        raise HTTPException(status_code=404, detail=str(e))

    nb_path = resolve_notebook_path(tpl)
    if not nb_path.exists():
        raise HTTPException(status_code=404, detail="Notebook file not found on disk")
    if nb_path.suffix != ".ipynb":
        raise HTTPException(status_code=400, detail="Invalid notebook file type")

    return FileResponse(
        path=str(nb_path),
        media_type="application/x-ipynb+json",
        filename=nb_path.name,
    )


@router.get("", response_model=List[Dict[str, Any]])
def list_notebooks() -> List[Dict[str, Any]]:
    """List available notebooks with Voila/Jupyter URLs.

    Only includes notebooks that exist on disk under the templates directory.
    """
    items = load_registry()
    out: List[Dict[str, Any]] = []
    jhub = os.environ.get("JUPYTERHUB_URL", "http://localhost:8000").rstrip("/")
    voila = os.environ.get("VOILA_URL", "").rstrip("/") or jhub

    for tpl in items:
        if not isinstance(tpl, dict):
            continue
        try:
            nb_path = resolve_notebook_path(tpl)
        except Exception:
            continue
        if not nb_path.exists():
            continue
        nb_rel = str(tpl.get("notebook_path") or "").strip()
        if not nb_rel:
            continue
        tpl_id = str(tpl.get("id") or "")
        out.append(
            {
                "path": nb_rel,
                "title": tpl.get("title") or tpl_id or nb_rel,
                "description": tpl.get("description") or "",
                "tags": tpl.get("tags") or [],
                "preview_image": tpl.get("preview_image"),
                "voila_url": f"{voila}/voila/render/{nb_rel}",
                "jupyter_url": f"{jhub}/notebooks/{nb_rel}",
                "source_url": f"/api/notebooks/templates/{tpl_id}/download" if tpl_id else None,
            }
        )
    return out


@router.get("/{path:path}/metadata", response_model=Dict[str, Any])
def get_notebook_metadata(path: str) -> Dict[str, Any]:
    """Lookup metadata for a notebook by its relative path within templates."""
    nb_rel = (path or "").strip()
    if not nb_rel:
        raise HTTPException(status_code=400, detail="path is required")

    jhub = os.environ.get("JUPYTERHUB_URL", "http://localhost:8000").rstrip("/")
    voila = os.environ.get("VOILA_URL", "").rstrip("/") or jhub

    for tpl in load_registry():
        if not isinstance(tpl, dict):
            continue
        if str(tpl.get("notebook_path") or "").strip() != nb_rel:
            continue
        nb_path = resolve_notebook_path(tpl)
        if not nb_path.exists():
            raise HTTPException(status_code=404, detail="Notebook file not found on disk")
        tpl_id = str(tpl.get("id") or "")
        return {
            "path": nb_rel,
            "title": tpl.get("title") or tpl_id or nb_rel,
            "description": tpl.get("description") or "",
            "tags": tpl.get("tags") or [],
            "preview_image": tpl.get("preview_image"),
            "voila_url": f"{voila}/voila/render/{nb_rel}",
            "jupyter_url": f"{jhub}/notebooks/{nb_rel}",
            "source_url": f"/api/notebooks/templates/{tpl_id}/download" if tpl_id else None,
        }
    raise HTTPException(status_code=404, detail="Notebook not found in registry")


__all__ = ["router"]


