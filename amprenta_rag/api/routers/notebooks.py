"""Notebook template registry API endpoints."""

from __future__ import annotations

from typing import Any, Dict, List, Optional

from fastapi import APIRouter, HTTPException, Query
from fastapi.responses import FileResponse

from amprenta_rag.notebooks.registry import get_template, resolve_notebook_path, search_templates


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


__all__ = ["router"]


