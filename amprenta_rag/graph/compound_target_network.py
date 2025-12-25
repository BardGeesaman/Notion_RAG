"""Compound-target network query helpers (GraphEdge -> Cytoscape-ready JSON)."""

from __future__ import annotations

import math
from typing import Any, Dict, Iterable, List, Optional
from uuid import UUID

from amprenta_rag.database.models import Compound, Feature, GraphEdge
from amprenta_rag.database.session import db_session
from amprenta_rag.utils.uuid_utils import ensure_uuid


def calculate_pic50(ic50_nm: float) -> float:
    """
    pIC50 = -log10(IC50_M) = 9 - log10(IC50_nM)
    """
    x = float(ic50_nm)
    if x <= 0:
        raise ValueError("ic50_nm must be > 0")
    return 9.0 - math.log10(x)


def edge_width_from_pic50(pic50: float) -> float:
    """
    Map pIC50 in [5, 9] to width in [1, 10] (clamped).
    """
    p = float(pic50)
    p = max(5.0, min(9.0, p))
    # 5 -> 1, 9 -> 10
    return 1.0 + (p - 5.0) * (9.0 / 4.0)


def _node_id(entity_type: str, entity_id: UUID) -> str:
    return f"{entity_type}:{entity_id}"


def _is_activator(activity_type: str | None) -> bool:
    t = (activity_type or "").lower()
    return any(k in t for k in ("ec50", "agonist", "activation", "activator"))


class CompoundTargetNetworkService:
    """
    Query GraphEdge for compound-target (compound->feature) activity edges and convert to Cytoscape JSON.
    """

    def __init__(self, session_factory=db_session):
        self._session_factory = session_factory

    def get_compound_target_network(
        self,
        compound_ids: Optional[Iterable[str | UUID]] = None,
        target_ids: Optional[Iterable[str | UUID]] = None,
        filters: Optional[Dict[str, Any]] = None,
        *,
        max_nodes: int = 500,
        limit_edges: int = 5000,
    ) -> Dict[str, Any]:
        """
        Returns:
            {"nodes": [...], "edges": [...], "meta": {...}}
        """
        filt = filters or {}
        ic50_range = filt.get("ic50_range") if isinstance(filt.get("ic50_range"), dict) else {}
        min_nm = ic50_range.get("min_nm")
        max_nm = ic50_range.get("max_nm")
        act_type = filt.get("activity_type")

        comp_u = [ensure_uuid(x) for x in (compound_ids or []) if ensure_uuid(x) is not None]
        tgt_u = [ensure_uuid(x) for x in (target_ids or []) if ensure_uuid(x) is not None]

        edges_out: List[Dict[str, Any]] = []
        node_compound_ids: List[UUID] = []
        node_target_ids: List[UUID] = []

        with self._session_factory() as db:
            q = (
                db.query(GraphEdge)
                .filter(GraphEdge.source_entity_type == "compound")
                .filter(GraphEdge.target_entity_type == "feature")
                .filter(GraphEdge.relationship_type == "activity_against")
            )
            if comp_u:
                q = q.filter(GraphEdge.source_entity_id.in_(comp_u))
            if tgt_u:
                q = q.filter(GraphEdge.target_entity_id.in_(tgt_u))

            q = q.order_by(GraphEdge.updated_at.desc()).limit(int(limit_edges))
            rows: List[GraphEdge] = list(q.all() or [])

            for e in rows:
                prov = e.provenance or {}
                best_nm = prov.get("best_ic50_nm")
                a_type = prov.get("activity_type")

                # Filter by activity_type (exact match) if requested.
                if act_type and str(a_type or "") != str(act_type):
                    continue

                # Filter by IC50 range if requested.
                if best_nm is not None:
                    try:
                        v = float(best_nm)
                    except Exception:
                        v = None
                    if v is not None:
                        if min_nm is not None and v < float(min_nm):
                            continue
                        if max_nm is not None and v > float(max_nm):
                            continue

                # Enforce max nodes during accumulation.
                if e.source_entity_id not in node_compound_ids:
                    if (len(node_compound_ids) + len(node_target_ids) + 1) > int(max_nodes):
                        break
                    node_compound_ids.append(e.source_entity_id)
                if e.target_entity_id not in node_target_ids:
                    if (len(node_compound_ids) + len(node_target_ids) + 1) > int(max_nodes):
                        break
                    node_target_ids.append(e.target_entity_id)

                # Styling metadata
                pic50 = None
                width = 1.0
                if best_nm is not None:
                    try:
                        pic50 = calculate_pic50(float(best_nm))
                        width = edge_width_from_pic50(pic50)
                    except Exception:
                        pic50 = None
                        width = 1.0

                edge_color = "#2E8B57" if _is_activator(str(a_type) if a_type is not None else None) else "#C73E1D"

                edges_out.append(
                    {
                        "data": {
                            "id": str(e.id),
                            "source": _node_id("compound", e.source_entity_id),
                            "target": _node_id("feature", e.target_entity_id),
                            "relationship_type": e.relationship_type,
                            "confidence": e.confidence,
                            "evidence_source": e.evidence_source,
                            "provenance": prov,
                            "best_ic50_nm": best_nm,
                            "activity_type": a_type,
                            "assays_count": prov.get("assays_count"),
                            "assay_ids": prov.get("assay_ids"),
                            "pic50": pic50,
                            "pic50_width": float(width),
                            "edge_color": edge_color,
                        }
                    }
                )

            # Labels
            comp_rows = (
                db.query(Compound.id, Compound.compound_id)
                .filter(Compound.id.in_(node_compound_ids))
                .all()
            )
            comp_label = {cid: (cstr or str(cid)[:8]) for cid, cstr in comp_rows}

            feat_rows = db.query(Feature.id, Feature.name).filter(Feature.id.in_(node_target_ids)).all()
            feat_label = {fid: (nm or str(fid)[:8]) for fid, nm in feat_rows}

        nodes_out: List[Dict[str, Any]] = []
        for cid in node_compound_ids:
            nodes_out.append(
                {
                    "data": {
                        "id": _node_id("compound", cid),
                        "entity_type": "compound",
                        "entity_id": str(cid),
                        "compound_id": comp_label.get(cid),
                        "label": comp_label.get(cid),
                        "shape": "ellipse",
                        "node_type": "compound",
                    }
                }
            )
        for tid in node_target_ids:
            nodes_out.append(
                {
                    "data": {
                        "id": _node_id("feature", tid),
                        "entity_type": "feature",
                        "entity_id": str(tid),
                        "target_name": feat_label.get(tid),
                        "label": feat_label.get(tid),
                        "shape": "diamond",
                        "node_type": "target",
                    }
                }
            )

        return {
            "nodes": nodes_out,
            "edges": edges_out,
            "meta": {"node_count": len(nodes_out), "edge_count": len(edges_out), "max_nodes": int(max_nodes)},
        }

    def expand_from_compound(self, compound_id: str | UUID) -> List[UUID]:
        cid = ensure_uuid(compound_id)
        if cid is None:
            raise ValueError("compound_id required")
        with self._session_factory() as db:
            rows = (
                db.query(GraphEdge.target_entity_id)
                .filter(GraphEdge.source_entity_type == "compound")
                .filter(GraphEdge.source_entity_id == cid)
                .filter(GraphEdge.relationship_type == "activity_against")
                .all()
            )
        return [r[0] for r in rows if r and r[0]]

    def expand_from_target(self, target_id: str | UUID) -> List[UUID]:
        tid = ensure_uuid(target_id)
        if tid is None:
            raise ValueError("target_id required")
        with self._session_factory() as db:
            rows = (
                db.query(GraphEdge.source_entity_id)
                .filter(GraphEdge.target_entity_type == "feature")
                .filter(GraphEdge.target_entity_id == tid)
                .filter(GraphEdge.relationship_type == "activity_against")
                .all()
            )
        return [r[0] for r in rows if r and r[0]]


__all__ = ["CompoundTargetNetworkService", "calculate_pic50", "edge_width_from_pic50"]


