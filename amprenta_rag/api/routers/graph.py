"""Evidence graph API (edge queries + traversal)."""

from __future__ import annotations

from typing import Any, Dict, List, Optional
from uuid import UUID

from fastapi import APIRouter, HTTPException, Query
from pydantic import BaseModel

from amprenta_rag.graph.cytoscape import to_cytoscape_json
from amprenta_rag.graph.edge_builder import EdgeBuilder, Direction
from amprenta_rag.graph.analytics import compute_graph_analytics
from amprenta_rag.graph.traversal import k_hop_subgraph, shortest_path
from amprenta_rag.connectivity.string_client import (
    get_interactions,
    get_interaction_partners,
    interactions_to_cytoscape,
    DEFAULT_SPECIES,
)


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


class TraverseRequest(BaseModel):
    entity_type: str
    entity_id: UUID
    depth: int = 2
    relationships: Optional[List[str]] = None


class GraphSubgraphResponse(BaseModel):
    nodes: List[Dict[str, Any]]
    edges: List[Dict[str, Any]]
    cytoscape: Dict[str, Any]
    truncated: bool = False


class PathRequest(BaseModel):
    source_type: str
    source_id: UUID
    target_type: str
    target_id: UUID
    relationships: Optional[List[str]] = None


class GraphPathResponse(BaseModel):
    found: bool
    nodes: List[Dict[str, Any]] = []
    edges: List[Dict[str, Any]] = []
    cytoscape: Optional[Dict[str, Any]] = None


class GraphAnalyticsRequest(BaseModel):
    nodes: List[Dict[str, Any]]
    edges: List[Dict[str, Any]]
    metrics: List[str]


class GraphAnalyticsResponse(BaseModel):
    degree_centrality: Optional[Dict[str, float]] = None
    communities: Optional[Dict[str, int]] = None
    community_sizes: Optional[Dict[str, int]] = None
    modularity: Optional[float] = None


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


@router.post("/traverse", response_model=GraphSubgraphResponse)
def traverse(payload: TraverseRequest) -> GraphSubgraphResponse:
    try:
        out = k_hop_subgraph(
            payload.entity_type,
            payload.entity_id,
            depth=payload.depth,
            relationships=payload.relationships,
            max_nodes=500,
        )
        nodes = [{"entity_type": t, "entity_id": str(i)} for (t, i) in out["nodes"]]
        edges = out["edges"]
        cy = to_cytoscape_json(out["nodes"], edges)
        return GraphSubgraphResponse(nodes=nodes, edges=edges, cytoscape=cy, truncated=bool(out.get("truncated")))
    except Exception as e:  # noqa: BLE001
        raise HTTPException(status_code=500, detail=f"Graph traversal failed: {e}")


@router.post("/path", response_model=GraphPathResponse)
def path(payload: PathRequest) -> GraphPathResponse:
    try:
        out = shortest_path(
            payload.source_type,
            payload.source_id,
            payload.target_type,
            payload.target_id,
            relationships=payload.relationships,
            timeout_s=5.0,
            max_nodes=500,
        )
        if out is None:
            return GraphPathResponse(found=False)
        nodes = [{"entity_type": t, "entity_id": str(i)} for (t, i) in out["nodes"]]
        edges = out["edges"]
        cy = to_cytoscape_json(out["nodes"], edges)
        return GraphPathResponse(found=True, nodes=nodes, edges=edges, cytoscape=cy)
    except Exception as e:  # noqa: BLE001
        raise HTTPException(status_code=500, detail=f"Graph path failed: {e}")


@router.post("/analytics", response_model=GraphAnalyticsResponse)
def analytics(payload: GraphAnalyticsRequest) -> GraphAnalyticsResponse:
    try:
        out = compute_graph_analytics(payload.nodes or [], payload.edges or [], payload.metrics or [])
        comms = out.get("communities") if isinstance(out, dict) else None
        sizes: Optional[Dict[str, int]] = None
        if isinstance(comms, dict) and comms:
            tmp: Dict[int, int] = {}
            for _, cid in comms.items():
                try:
                    tmp[int(cid)] = tmp.get(int(cid), 0) + 1
                except Exception:
                    continue
            sizes = {str(k): v for k, v in tmp.items()}
        return GraphAnalyticsResponse(
            degree_centrality=out.get("degree_centrality"),
            communities=out.get("communities"),
            modularity=out.get("modularity"),
            community_sizes=sizes,
        )
    except Exception as e:  # noqa: BLE001
        raise HTTPException(status_code=500, detail=f"Graph analytics failed: {e}")


# --- Protein-Protein Interaction Schemas ---

class PPIRequest(BaseModel):
    """Request for protein-protein interaction network."""
    proteins: List[str]
    species: int = DEFAULT_SPECIES
    min_score: int = 400


class PPIResponse(BaseModel):
    """Response for protein-protein interaction network."""
    nodes: List[Dict[str, Any]]
    edges: List[Dict[str, Any]]
    interaction_count: int


# --- Protein-Protein Interaction Endpoints ---

@router.get("/ppi/{gene_symbol}")
def get_ppi_partners(
    gene_symbol: str,
    species: int = Query(DEFAULT_SPECIES),
    limit: int = Query(50, ge=1, le=200),
    min_score: int = Query(400, ge=0, le=1000),
) -> PPIResponse:
    """Get PPI network for a single gene."""
    interactions = get_interaction_partners(gene_symbol, species, limit, min_score)
    nodes, edges = interactions_to_cytoscape(interactions)
    return PPIResponse(nodes=nodes, edges=edges, interaction_count=len(interactions))


@router.post("/ppi/network")
def get_ppi_network(request: PPIRequest) -> PPIResponse:
    """Get PPI network for a list of genes."""
    interactions = get_interactions(request.proteins, request.species, request.min_score)
    nodes, edges = interactions_to_cytoscape(interactions)
    return PPIResponse(nodes=nodes, edges=edges, interaction_count=len(interactions))


__all__ = ["router"]


