"""Evidence graph API (edge queries)."""

from __future__ import annotations

from typing import Any, Dict, List, Optional
from uuid import UUID

from fastapi import APIRouter, HTTPException, Query
from pydantic import BaseModel

from amprenta_rag.graph.edge_builder import EdgeBuilder, Direction


router = APIRouter(prefix="/graph", tags=["Graph"])


class NeighborResponse(BaseModel):
    source_entity_type: str
    source_entity_id: UUID
    target_entity_type: str
    target_entity_id: UUID
    relationship_type: str
    confidence: Optional[float] = None
    evidence_source: Optional[str] = None
    provenance: Optional[Dict[str, Any]] = None


@router.get("/neighbors", response_model=List[NeighborResponse])
def get_neighbors(
    entity_type: str = Query(...),
    entity_id: UUID = Query(...),
    direction: Direction = Query("both"),
    relationship_type: Optional[List[str]] = Query(None),
    min_confidence: float = Query(0.0, ge=0.0, le=1.0),
    limit: int = Query(200, ge=1, le=5000),
) -> List[NeighborResponse]:
    try:
        builder = EdgeBuilder()
        neighbors = builder.get_neighbors(
            entity_type=entity_type,
            entity_id=entity_id,
            direction=direction,
            relationship_types=relationship_type,
            min_confidence=min_confidence,
            limit=limit,
        )
        return [
            NeighborResponse(
                source_entity_type=n.source_entity_type,
                source_entity_id=n.source_entity_id,
                target_entity_type=n.target_entity_type,
                target_entity_id=n.target_entity_id,
                relationship_type=n.relationship_type,
                confidence=n.confidence,
                evidence_source=n.evidence_source,
                provenance=n.provenance,
            )
            for n in neighbors
        ]
    except Exception as e:  # noqa: BLE001
        raise HTTPException(status_code=500, detail=f"Graph neighbor query failed: {e}")


__all__ = ["router"]


