"""Entity versioning service for provenance tracking."""

import hashlib
import json
from typing import Any, Dict, List, Optional
from uuid import UUID

from sqlalchemy.orm import Session
from sqlalchemy import func

from amprenta_rag.database.models import EntityVersion
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)

# Entities that support versioning
VERSIONABLE_ENTITIES = ["dataset", "experiment"]


def compute_snapshot_checksum(data: dict) -> str:
    """Compute SHA256 checksum of snapshot data."""
    json_str = json.dumps(data, sort_keys=True, default=str)
    return hashlib.sha256(json_str.encode()).hexdigest()


def entity_to_snapshot(entity) -> dict:
    """Generic serializer for SQLAlchemy models."""
    return {c.name: getattr(entity, c.name) for c in entity.__table__.columns}


def create_version(
    db: Session,
    entity_type: str,
    entity_id: UUID,
    data: Dict[str, Any],
    user_id: Optional[UUID] = None,
    change_summary: Optional[str] = None,
) -> EntityVersion:
    """Create a new version snapshot."""
    if entity_type not in VERSIONABLE_ENTITIES:
        raise ValueError(f"Entity type '{entity_type}' is not versionable")
    
    # Get next version number
    max_version = (
        db.query(func.max(EntityVersion.version_number))
        .filter(EntityVersion.entity_type == entity_type, EntityVersion.entity_id == entity_id)
        .scalar()
    ) or 0
    
    checksum = compute_snapshot_checksum(data)
    
    version = EntityVersion(
        entity_type=entity_type,
        entity_id=entity_id,
        version_number=max_version + 1,
        data_snapshot=data,
        checksum_sha256=checksum,
        created_by=user_id,
        change_summary=change_summary,
    )
    db.add(version)
    db.commit()
    db.refresh(version)
    
    logger.info(f"Created version {version.version_number} for {entity_type}/{entity_id}")
    return version


def get_versions(
    db: Session,
    entity_type: str,
    entity_id: UUID,
) -> List[EntityVersion]:
    """Get all versions for an entity, newest first."""
    return (
        db.query(EntityVersion)
        .filter(EntityVersion.entity_type == entity_type, EntityVersion.entity_id == entity_id)
        .order_by(EntityVersion.version_number.desc())
        .all()
    )


def get_version(db: Session, version_id: UUID) -> Optional[EntityVersion]:
    """Get a specific version by ID."""
    return db.query(EntityVersion).filter(EntityVersion.id == version_id).first()


def get_latest_version(
    db: Session,
    entity_type: str,
    entity_id: UUID,
) -> Optional[EntityVersion]:
    """Get the latest version for an entity."""
    return (
        db.query(EntityVersion)
        .filter(EntityVersion.entity_type == entity_type, EntityVersion.entity_id == entity_id)
        .order_by(EntityVersion.version_number.desc())
        .first()
    )


def compare_versions(v1_data: dict, v2_data: dict) -> dict:
    """Compare two version snapshots, return diff."""
    diff = {"added": {}, "removed": {}, "changed": {}}
    all_keys = set(v1_data.keys()) | set(v2_data.keys())
    
    for key in all_keys:
        if key not in v1_data:
            diff["added"][key] = v2_data[key]
        elif key not in v2_data:
            diff["removed"][key] = v1_data[key]
        elif v1_data[key] != v2_data[key]:
            diff["changed"][key] = {"old": v1_data[key], "new": v2_data[key]}
    
    return diff


def get_version_history(
    db: Session,
    entity_type: str,
    entity_id: UUID,
    include_diffs: bool = False,
) -> List[Dict[str, Any]]:
    """Get version history with optional diffs between consecutive versions."""
    versions = get_versions(db, entity_type, entity_id)
    
    history = []
    for i, version in enumerate(versions):
        version_info = {
            "id": version.id,
            "version_number": version.version_number,
            "created_at": version.created_at,
            "created_by": version.created_by,
            "change_summary": version.change_summary,
            "checksum": version.checksum_sha256,
        }
        
        if include_diffs and i < len(versions) - 1:
            # Compare with previous version (next in list since ordered desc)
            prev_version = versions[i + 1]
            version_info["diff"] = compare_versions(
                prev_version.data_snapshot,
                version.data_snapshot
            )
        
        history.append(version_info)
    
    return history


def rollback_to_version(
    db: Session,
    version_id: UUID,
    user_id: Optional[UUID] = None,
) -> EntityVersion:
    """Create a new version by rolling back to a previous version's data."""
    source_version = get_version(db, version_id)
    if not source_version:
        raise ValueError(f"Version {version_id} not found")
    
    # Create new version with the old data
    new_version = create_version(
        db=db,
        entity_type=source_version.entity_type,
        entity_id=source_version.entity_id,
        data=source_version.data_snapshot,
        user_id=user_id,
        change_summary=f"Rolled back to version {source_version.version_number}",
    )
    
    logger.info(
        f"Rolled back {source_version.entity_type}/{source_version.entity_id} "
        f"to version {source_version.version_number} (new version {new_version.version_number})"
    )
    
    return new_version


def delete_version(db: Session, version_id: UUID) -> bool:
    """Delete a specific version (use with caution)."""
    version = get_version(db, version_id)
    if not version:
        return False
    
    # Check if this is the only version
    other_versions = (
        db.query(EntityVersion)
        .filter(
            EntityVersion.entity_type == version.entity_type,
            EntityVersion.entity_id == version.entity_id,
            EntityVersion.id != version_id,
        )
        .count()
    )
    
    if other_versions == 0:
        raise ValueError("Cannot delete the only version of an entity")
    
    db.delete(version)
    db.commit()
    
    logger.warning(
        f"Deleted version {version.version_number} of {version.entity_type}/{version.entity_id}"
    )
    
    return True
