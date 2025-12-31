"""Pathway maps API (KEGG-only MVP).

Provides KGML-based pathway structure, expression overlays, and convenience endpoints
for searching KEGG pathways and returning enriched pathways for a dataset.
"""

from __future__ import annotations

import asyncio
import threading
import time
from typing import Dict, List, Optional, Set
from uuid import UUID

import requests
from fastapi import APIRouter

from amprenta_rag.analysis.pathway.enrichment import perform_pathway_enrichment
from amprenta_rag.analysis.pathway.overlay import compute_overlay
from amprenta_rag.analysis.pathway.structure import edge_style, get_pathway_structure
from amprenta_rag.api import schemas
from amprenta_rag.database.models import Feature, dataset_feature_assoc
from amprenta_rag.database.session import db_session


router = APIRouter(prefix="/pathway-maps", tags=["pathway-maps"])

KEGG_BASE_URL = "https://rest.kegg.jp"
KEGG_SEARCH_DELAY_SECONDS = 0.5
_KEGG_LOCK = threading.Lock()
_LAST_KEGG_SEARCH_TS: float = 0.0


def _rate_limit_kegg_search() -> None:
    """Thread-safe rate limiting for KEGG API calls."""
    global _LAST_KEGG_SEARCH_TS
    with _KEGG_LOCK:  # Thread-safe access to global state
        now = time.monotonic()
        wait = KEGG_SEARCH_DELAY_SECONDS - (now - _LAST_KEGG_SEARCH_TS)
        if wait > 0:
            time.sleep(wait)
        _LAST_KEGG_SEARCH_TS = time.monotonic()


# Sync helper functions for external API calls
def _sync_get_pathway_structure(pathway_id: str):
    """Sync helper for pathway structure fetch via KEGG API."""
    return get_pathway_structure(pathway_id)


def _sync_search_kegg_pathways(query: str):
    """Sync helper for KEGG pathway search API."""
    q = (query or "").strip()
    if not q:
        return []
    
    _rate_limit_kegg_search()
    url = f"{KEGG_BASE_URL}/find/pathway/{q}"
    
    try:
        resp = requests.get(url, timeout=15)
        if resp.status_code != 200:
            return []
        
        results = []
        for line in resp.text.strip().split("\n"):
            if not line:
                continue
            parts = line.split("\t", 1)
            if len(parts) < 2:
                continue
            pid_raw = parts[0].strip()
            pid = pid_raw.split(":", 1)[1] if ":" in pid_raw else pid_raw
            name = parts[1].strip()
            results.append({
                "pathway_id": pid,
                "name": name,
                "organism": _organism_from_pathway_id(pid),
                "gene_count": 0,
            })
            if len(results) >= 50:
                break
        return results
    except Exception:
        return []


def _sync_perform_pathway_enrichment(features: Set[str], types: Set[str]):
    """Sync helper for pathway enrichment via KEGG API."""
    return perform_pathway_enrichment(
        input_features=features,
        input_feature_types=types or {"gene"},
        pathway_sources=["KEGG"],
        p_value_threshold=0.05,
    )


def _organism_from_pathway_id(pathway_id: str) -> str:
    pid = (pathway_id or "").strip()
    return pid[:3] if len(pid) >= 3 else ""


def _load_dataset_features(dataset_id: UUID) -> tuple[Set[str], Set[str], List[Feature]]:
    with db_session() as db:
        feats: List[Feature] = (
            db.query(Feature)
            .join(dataset_feature_assoc, Feature.id == dataset_feature_assoc.c.feature_id)
            .filter(dataset_feature_assoc.c.dataset_id == dataset_id)
            .all()
        )
    features: Set[str] = set()
    types: Set[str] = set()
    for f in feats:
        name = getattr(f, "name", None)
        ftype = getattr(f, "feature_type", None) or getattr(f, "type", None)
        if name:
            features.add(str(name))
        if ftype:
            t = str(ftype)
            if t in ("gene", "rna", "transcript"):
                types.add("gene")
            elif t in ("protein", "proteomics"):
                types.add("protein")
            elif t in ("metabolite", "metabolomics"):
                types.add("metabolite")
            elif t in ("lipid", "lipidomics"):
                types.add("lipid")
    return features, types, feats


@router.get("/structure/{pathway_id}", response_model=schemas.PathwayStructureResponse)
async def get_structure(pathway_id: str) -> schemas.PathwayStructureResponse:
    st = await asyncio.to_thread(_sync_get_pathway_structure, pathway_id)
    nodes = [
        schemas.PathwayNodeSchema(
            id=n.id,
            name=n.name,
            type=n.type,
            x=float(n.x),
            y=float(n.y),
            kegg_ids=list(n.kegg_ids or []),
        )
        for n in (st.nodes or [])
    ]
    edges = []
    for e in st.edges or []:
        style = edge_style(e)
        edges.append(
            schemas.PathwayEdgeSchema(
                source=e.source,
                target=e.target,
                type=e.type,
                subtype=e.subtype,
                style=style.get("style", "solid"),
                color=style.get("color", "#9E9E9E"),
            )
        )
    return schemas.PathwayStructureResponse(
        pathway_id=st.pathway_id,
        name=st.name,
        nodes=nodes,
        edges=edges,
        organism=st.organism or _organism_from_pathway_id(st.pathway_id),
    )


@router.post("/overlay", response_model=schemas.PathwayOverlayResponse)
def build_overlay(request: schemas.PathwayOverlayRequest) -> schemas.PathwayOverlayResponse:
    st = get_pathway_structure(request.pathway_id)
    overlays = compute_overlay(
        structure=st,
        expression_data=request.expression_data or {},
        colormap=request.colormap,
        vmin=float(request.vmin),
        vmax=float(request.vmax),
    )
    out = [
        schemas.PathwayOverlayNodeSchema(
            node_id=o.node_id,
            gene_symbol=o.gene_symbol,
            value=float(o.value),
            color=o.color,
            label=o.label,
        )
        for o in overlays
    ]
    return schemas.PathwayOverlayResponse(pathway_id=str(request.pathway_id), overlays=out)


@router.get("/search", response_model=List[schemas.PathwaySearchResult])
async def search_pathways(query: str) -> List[schemas.PathwaySearchResult]:
    # Search pathways using async thread pool
    search_results = await asyncio.to_thread(_sync_search_kegg_pathways, query)
    
    # Convert to response objects
    results = [
        schemas.PathwaySearchResult(**result_data)
        for result_data in search_results
    ]
    
    return results


@router.get("/enriched/{dataset_id}", response_model=List[schemas.PathwaySearchResult])
async def enriched_pathways(dataset_id: UUID) -> List[schemas.PathwaySearchResult]:
    features, types, _ = _load_dataset_features(dataset_id)
    if not features:
        return []
    
    # Perform pathway enrichment using async thread pool
    results = await asyncio.to_thread(
        _sync_perform_pathway_enrichment,
        features,
        types
    )
    
    out: List[schemas.PathwaySearchResult] = []
    for r in results[:20]:
        pid = r.pathway.pathway_id
        out.append(
            schemas.PathwaySearchResult(
                pathway_id=pid,
                name=r.pathway.name,
                organism=_organism_from_pathway_id(pid),
                gene_count=int(r.input_features),
            )
        )
    return out


def _extract_log2fc(feature: Feature) -> Optional[float]:
    # The spec calls this Feature.metadata["log2fc"] or similar; not all deployments
    # have a metadata column on Feature, so we probe common attribute names.
    md = getattr(feature, "metadata", None)
    if isinstance(md, dict):
        for k in ("log2fc", "log2FC", "log2_fold_change", "log2FoldChange", "fold_change", "fc"):
            if k in md:
                try:
                    return float(md[k])
                except Exception:
                    continue
    ext = getattr(feature, "external_ids", None)
    if isinstance(ext, dict):
        for k in ("log2fc", "log2FC", "log2_fold_change", "log2FoldChange", "fold_change", "fc"):
            if k in ext:
                try:
                    return float(ext[k])
                except Exception:
                    continue
    return None


@router.get("/expression/{dataset_id}", response_model=schemas.PathwayExpressionResponse)
def dataset_expression(dataset_id: UUID) -> schemas.PathwayExpressionResponse:
    # Graceful: missing dataset or missing values -> empty mapping.
    _, _, feats = _load_dataset_features(dataset_id)
    gene_expression: Dict[str, float] = {}
    for f in feats:
        name = getattr(f, "name", None)
        if not name:
            continue
        val = _extract_log2fc(f)
        if val is None:
            continue
        gene_expression[str(name)] = float(val)
    return schemas.PathwayExpressionResponse(dataset_id=dataset_id, gene_expression=gene_expression)


__all__ = ["router"]


