"""Artifact registry utilities for generated reports."""

from __future__ import annotations

from typing import List, Optional
from uuid import UUID

from amprenta_rag.database.models import ReportArtifact
from amprenta_rag.database.session import db_session


def save_artifact(
    entity_type: str,
    entity_id: UUID,
    format: str,
    file_path: str,
    params_hash: str,
    user_id: Optional[UUID] = None,
) -> ReportArtifact:
    """Persist a report artifact record."""
    with db_session() as db:
        artifact = ReportArtifact(
            entity_type=entity_type,
            entity_id=entity_id,  # type: ignore[arg-type]
            format=format,
            file_path=file_path,
            params_hash=params_hash,
            created_by_id=user_id,  # type: ignore[arg-type]
        )
        db.add(artifact)
        db.commit()
        db.refresh(artifact)
        return artifact


def get_artifact(artifact_id: UUID) -> Optional[ReportArtifact]:
    """Fetch a report artifact by ID."""
    with db_session() as db:
        return db.query(ReportArtifact).filter(ReportArtifact.id == artifact_id).first()


def list_artifacts(
    entity_type: str,
    entity_id: UUID,
    limit: int = 50,
) -> List[ReportArtifact]:
    """List artifacts for an entity, newest first."""
    with db_session() as db:
        return (
            db.query(ReportArtifact)
            .filter(
                ReportArtifact.entity_type == entity_type,
                ReportArtifact.entity_id == entity_id,
            )
            .order_by(ReportArtifact.created_at.desc())
            .limit(limit)
            .all()
        )


def get_cached_artifact(
    entity_type: str,
    entity_id: UUID,
    params_hash: str,
) -> Optional[ReportArtifact]:
    """Return most recent artifact matching cache key."""
    with db_session() as db:
        return (
            db.query(ReportArtifact)
            .filter(
                ReportArtifact.entity_type == entity_type,
                ReportArtifact.entity_id == entity_id,
                ReportArtifact.params_hash == params_hash,
            )
            .order_by(ReportArtifact.created_at.desc())
            .first()
        )


def delete_artifact(artifact_id: UUID) -> bool:
    """Delete an artifact by ID. Returns True if deleted."""
    with db_session() as db:
        artifact = db.query(ReportArtifact).filter(ReportArtifact.id == artifact_id).first()
        if not artifact:
            return False
        db.delete(artifact)
        db.commit()
        return True

