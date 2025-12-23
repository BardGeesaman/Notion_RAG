"""EdgeBuilder service for the Evidence Graph.

MVP: Postgres-backed edge CRUD + neighbor queries.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, Dict, Iterable, List, Literal, Optional
from uuid import UUID

import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import insert as pg_insert

from amprenta_rag.database.models import GraphEdge
from amprenta_rag.database.session import db_session


Direction = Literal["out", "in", "both"]


@dataclass(frozen=True)
class Neighbor:
    source_entity_type: str
    source_entity_id: UUID
    target_entity_type: str
    target_entity_id: UUID
    relationship_type: str
    confidence: Optional[float]
    evidence_source: Optional[str]
    provenance: Optional[Dict[str, Any]]

    def asdict(self) -> Dict[str, Any]:
        d = asdict(self)
        d["source_entity_id"] = str(d["source_entity_id"])
        d["target_entity_id"] = str(d["target_entity_id"])
        return d


class EdgeBuilder:
    """Create edges and query neighbors."""

    def __init__(self, session_factory=db_session):
        self._session_factory = session_factory

    def create_edge(
        self,
        *,
        source_entity_type: str,
        source_entity_id: UUID,
        target_entity_type: str,
        target_entity_id: UUID,
        relationship_type: str,
        confidence: Optional[float] = None,
        evidence_source: Optional[str] = None,
        provenance: Optional[Dict[str, Any]] = None,
    ) -> GraphEdge:
        """Create or update a graph edge (idempotent on unique key)."""
        with self._session_factory() as db:
            # Prefer a Postgres upsert so we always bump updated_at (and remain race-safe).
            try:
                bind = db.get_bind()  # type: ignore[attr-defined]
                dialect_name = getattr(getattr(bind, "dialect", None), "name", None)
            except Exception:
                dialect_name = None

            if dialect_name == "postgresql" and hasattr(db, "execute"):
                now = sa.func.now()
                insert_stmt = pg_insert(GraphEdge).values(
                    source_entity_type=source_entity_type,
                    source_entity_id=source_entity_id,
                    target_entity_type=target_entity_type,
                    target_entity_id=target_entity_id,
                    relationship_type=relationship_type,
                    confidence=confidence,
                    evidence_source=evidence_source,
                    provenance=provenance,
                    created_at=now,
                    updated_at=now,
                )

                excluded = insert_stmt.excluded

                # Keep max confidence when both present. If excluded.confidence is NULL, keep existing.
                new_conf = sa.case(
                    (excluded.confidence.is_(None), GraphEdge.confidence),
                    else_=sa.func.greatest(sa.func.coalesce(GraphEdge.confidence, 0.0), excluded.confidence),
                )

                # If new provenance provided, replace; otherwise keep existing.
                new_prov = sa.case(
                    (excluded.provenance.is_(None), GraphEdge.provenance),
                    else_=excluded.provenance,
                )

                upsert = insert_stmt.on_conflict_do_update(
                    constraint="uq_graph_edge",
                    set_={
                        "confidence": new_conf,
                        "provenance": new_prov,
                        "updated_at": now,  # required by P2: bump updated_at on conflict
                    },
                ).returning(GraphEdge)

                row = db.execute(upsert).first()
                db.commit()
                # row may be a GraphEdge (ORM returning) or a Row with GraphEdge at index 0
                if row is None:
                    # Fallback to ORM lookup
                    return (
                        db.query(GraphEdge)
                        .filter_by(
                            source_entity_type=source_entity_type,
                            source_entity_id=source_entity_id,
                            target_entity_type=target_entity_type,
                            target_entity_id=target_entity_id,
                            relationship_type=relationship_type,
                            evidence_source=evidence_source,
                        )
                        .first()
                    )
                if isinstance(row, GraphEdge):
                    return row
                if hasattr(row, "_mapping") and len(row._mapping) == 1:
                    return list(row._mapping.values())[0]
                return row[0]

            existing: Optional[GraphEdge] = (
                db.query(GraphEdge)
                .filter_by(
                    source_entity_type=source_entity_type,
                    source_entity_id=source_entity_id,
                    target_entity_type=target_entity_type,
                    target_entity_id=target_entity_id,
                    relationship_type=relationship_type,
                    evidence_source=evidence_source,
                )
                .first()
            )

            if existing is not None:
                # Keep the max confidence if both present.
                if confidence is not None:
                    if existing.confidence is None or confidence > existing.confidence:
                        existing.confidence = confidence
                # Shallow merge provenance.
                if provenance:
                    merged = dict(existing.provenance or {})
                    for k, v in provenance.items():
                        merged[k] = v
                    existing.provenance = merged
                try:
                    existing.updated_at = sa.func.now()  # type: ignore[assignment]
                except Exception:
                    pass
                db.commit()
                return existing

            edge = GraphEdge(
                source_entity_type=source_entity_type,
                source_entity_id=source_entity_id,
                target_entity_type=target_entity_type,
                target_entity_id=target_entity_id,
                relationship_type=relationship_type,
                confidence=confidence,
                evidence_source=evidence_source,
                provenance=provenance,
            )
            db.add(edge)
            db.commit()
            try:
                db.refresh(edge)
            except Exception:
                pass
            return edge

    def get_neighbors(
        self,
        *,
        entity_type: str,
        entity_id: UUID,
        direction: Direction = "both",
        relationship_types: Optional[Iterable[str]] = None,
        min_confidence: float = 0.0,
        limit: int = 200,
    ) -> List[Neighbor]:
        """Get neighbors (incoming/outgoing/both) for an entity."""
        rel_set = set(relationship_types) if relationship_types else None

        edges: List[GraphEdge] = []
        with self._session_factory() as db:
            if direction in ("out", "both"):
                q = db.query(GraphEdge).filter_by(source_entity_type=entity_type, source_entity_id=entity_id)
                edges.extend(q.all())
            if direction in ("in", "both"):
                q = db.query(GraphEdge).filter_by(target_entity_type=entity_type, target_entity_id=entity_id)
                edges.extend(q.all())

        neighbors: List[Neighbor] = []
        for e in edges:
            if rel_set is not None and e.relationship_type not in rel_set:
                continue
            conf = float(e.confidence) if e.confidence is not None else None
            if conf is not None and conf < float(min_confidence):
                continue
            neighbors.append(
                Neighbor(
                    source_entity_type=e.source_entity_type,
                    source_entity_id=e.source_entity_id,
                    target_entity_type=e.target_entity_type,
                    target_entity_id=e.target_entity_id,
                    relationship_type=e.relationship_type,
                    confidence=e.confidence,
                    evidence_source=e.evidence_source,
                    provenance=e.provenance,
                )
            )

        # Stable-ish ordering: higher confidence first, then relationship type.
        neighbors.sort(key=lambda n: (-(n.confidence or 0.0), n.relationship_type))
        return neighbors[: int(limit)]


__all__ = ["EdgeBuilder", "Neighbor", "Direction"]


