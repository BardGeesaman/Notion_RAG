"""Phenotype (HPO) API endpoints."""

from __future__ import annotations

from typing import Any, Dict, List

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel

from amprenta_rag.phenotype.phenotype_mapper import get_mapper

router = APIRouter(prefix="/phenotypes", tags=["Phenotypes"])


class ExpandQueryRequest(BaseModel):
    query: str


class ExpandQueryResponse(BaseModel):
    hpo_ids: List[str]
    genes: List[str]
    gene_count: int


@router.get("/{hpo_id}/genes", response_model=List[str])
def genes_for_hpo(hpo_id: str) -> List[str]:
    try:
        mapper = get_mapper()
        return mapper.get_genes_for_hpo(hpo_id)
    except Exception as e:  # noqa: BLE001
        raise HTTPException(status_code=500, detail=f"Phenotype lookup failed: {e}")


@router.post("/expand-query", response_model=ExpandQueryResponse)
def expand_query(payload: ExpandQueryRequest) -> ExpandQueryResponse:
    try:
        mapper = get_mapper()
        out: Dict[str, Any] = mapper.expand_query(payload.query)
        return ExpandQueryResponse(**out)
    except Exception as e:  # noqa: BLE001
        raise HTTPException(status_code=500, detail=f"Phenotype expand failed: {e}")


__all__ = ["router"]


