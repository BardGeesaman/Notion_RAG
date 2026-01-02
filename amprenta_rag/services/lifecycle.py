"""Data lifecycle management service.

Provides:
- Cascade impact calculation for deletions
- Status updates with audit trail
- Bulk operations with confirmation workflow
"""

from __future__ import annotations

import hashlib
import json
import zipfile
from datetime import datetime, timezone
from io import BytesIO
from typing import Any, Dict, List, Optional, Tuple
from uuid import UUID

from sqlalchemy import func, select
from sqlalchemy.orm import Session

from amprenta_rag.database.models import (
    Dataset, Experiment, Signature, AuditLog, LifecycleStatus,
    dataset_feature_assoc, dataset_signature_assoc, experiment_dataset_assoc,
    signature_feature_assoc,
)
from amprenta_rag.models.chemistry import Compound
from amprenta_rag.database.session import db_session
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)

# Entity type to model mapping
ENTITY_MODELS = {
    "dataset": Dataset,
    "experiment": Experiment,
    "compound": Compound,
    "signature": Signature,
}


def get_entity_by_id(
    entity_type: str,
    entity_id: UUID,
    db: Optional[Session] = None
) -> Optional[Any]:
    """Get entity by type and ID."""
    model = ENTITY_MODELS.get(entity_type)
    if not model:
        return None
    
    def _query(session: Session):
        return session.query(model).filter(model.id == entity_id).first()
    
    if db:
        return _query(db)
    
    with db_session() as session:
        entity = _query(session)
        if entity:
            session.expunge(entity)
        return entity


def calculate_deletion_impact(
    entity_type: str,
    entity_id: UUID,
) -> Dict[str, Any]:
    """
    Calculate cascade impact of deleting an entity.
    
    Returns counts of related entities and blocking references.
    """
    impact = {
        "entity": {"type": entity_type, "id": str(entity_id), "name": None},
        "impact": {},
        "blocking_references": [],
        "can_delete": True,
    }
    
    with db_session() as db:
        entity = get_entity_by_id(entity_type, entity_id, db)
        if not entity:
            impact["error"] = "Entity not found"
            impact["can_delete"] = False
            return impact
        
        impact["entity"]["name"] = getattr(entity, "name", None) or getattr(entity, "title", str(entity_id))
        
        if entity_type == "dataset":
            # Count related entities
            impact["impact"]["features"] = db.query(func.count()).select_from(
                dataset_feature_assoc
            ).filter(dataset_feature_assoc.c.dataset_id == entity_id).scalar() or 0
            
            impact["impact"]["signatures"] = db.query(func.count()).select_from(
                dataset_signature_assoc
            ).filter(dataset_signature_assoc.c.dataset_id == entity_id).scalar() or 0
            
            impact["impact"]["experiments"] = db.query(func.count()).select_from(
                experiment_dataset_assoc
            ).filter(experiment_dataset_assoc.c.dataset_id == entity_id).scalar() or 0
            
            # Check for blocking references (ML models using this dataset)
            from amprenta_rag.database.models import MLModelMetadata
            blocking_models = db.query(MLModelMetadata).filter(
                MLModelMetadata.training_dataset_id == entity_id
            ).all()
            
            for model in blocking_models:
                impact["blocking_references"].append({
                    "type": "ml_model",
                    "id": str(model.id),
                    "name": model.model_name,
                    "reason": "training_dataset_id",
                })
                impact["can_delete"] = False
        
        elif entity_type == "experiment":
            impact["impact"]["datasets"] = db.query(func.count()).select_from(
                experiment_dataset_assoc
            ).filter(experiment_dataset_assoc.c.experiment_id == entity_id).scalar() or 0
        
        elif entity_type == "compound":
            # Check for HTS results, docking poses, etc.
            from amprenta_rag.database.models import ChEMBLActivity, PubChemBioassay
            impact["impact"]["chembl_activities"] = db.query(ChEMBLActivity).filter(
                ChEMBLActivity.compound_id == entity_id
            ).count()
            impact["impact"]["pubchem_bioassays"] = db.query(PubChemBioassay).filter(
                PubChemBioassay.compound_id == entity_id
            ).count()
        
        elif entity_type == "signature":
            impact["impact"]["datasets"] = db.query(func.count()).select_from(
                dataset_signature_assoc
            ).filter(dataset_signature_assoc.c.signature_id == entity_id).scalar() or 0
    
    if impact["blocking_references"]:
        impact["blocking_reason"] = f"{len(impact['blocking_references'])} blocking reference(s) found"
    
    return impact


def update_lifecycle_status(
    entity_type: str,
    entity_id: UUID,
    new_status: str,
    reason: Optional[str] = None,
    actor_id: Optional[UUID] = None,
) -> Tuple[bool, str]:
    """
    Update entity lifecycle status with audit trail.
    
    Returns (success, message).
    """
    if new_status not in [s.value for s in LifecycleStatus]:
        return False, f"Invalid status: {new_status}"
    
    model = ENTITY_MODELS.get(entity_type)
    if not model:
        return False, f"Unknown entity type: {entity_type}"
    
    with db_session() as db:
        entity = db.query(model).filter(model.id == entity_id).first()
        if not entity:
            return False, "Entity not found"
        
        old_status = entity.lifecycle_status
        if old_status == new_status:
            return True, "Status unchanged"
        
        # Update status
        entity.lifecycle_status = new_status
        
        # Update archived_at if transitioning to/from archived
        if hasattr(entity, 'archived_at'):
            if new_status == LifecycleStatus.ARCHIVED.value:
                entity.archived_at = datetime.now(timezone.utc)
            elif old_status == LifecycleStatus.ARCHIVED.value:
                entity.archived_at = None
        
        # Create audit log entry
        audit_entry = AuditLog(
            entity_type=entity_type,
            entity_id=entity_id,
            action="lifecycle_status_change",
            old_value={"lifecycle_status": old_status},
            new_value={"lifecycle_status": new_status, "reason": reason},
            user_id=actor_id,
        )
        db.add(audit_entry)
        
        db.commit()
        
        logger.info(
            f"Lifecycle status changed: {entity_type}:{entity_id} "
            f"{old_status} → {new_status} (reason: {reason})"
        )
        
        return True, f"Status updated: {old_status} → {new_status}"


def bulk_update_status(
    entity_type: str,
    entity_ids: List[UUID],
    new_status: str,
    reason: Optional[str] = None,
    actor_id: Optional[UUID] = None,
) -> Dict[str, Any]:
    """Bulk update lifecycle status for multiple entities."""
    results = {
        "total": len(entity_ids),
        "success": 0,
        "failed": 0,
        "errors": [],
    }
    
    for entity_id in entity_ids:
        success, message = update_lifecycle_status(
            entity_type, entity_id, new_status, reason, actor_id
        )
        if success:
            results["success"] += 1
        else:
            results["failed"] += 1
            results["errors"].append({"id": str(entity_id), "error": message})
    
    return results


def bulk_delete_preview(
    entity_type: str,
    entity_ids: List[UUID],
) -> Dict[str, Any]:
    """Preview cascade impact for bulk deletion (dry run)."""
    preview = {
        "entity_type": entity_type,
        "entity_count": len(entity_ids),
        "total_impact": {},
        "blocking_entities": [],
        "can_proceed": True,
    }
    
    for entity_id in entity_ids:
        impact = calculate_deletion_impact(entity_type, entity_id)
        
        # Aggregate impact counts
        for key, count in impact.get("impact", {}).items():
            preview["total_impact"][key] = preview["total_impact"].get(key, 0) + count
        
        # Collect blocking references
        if impact.get("blocking_references"):
            preview["blocking_entities"].append({
                "entity_id": str(entity_id),
                "blocking": impact["blocking_references"],
            })
            preview["can_proceed"] = False
    
    return preview


def execute_bulk_archive(
    entity_type: str,
    entity_ids: List[UUID],
    reason: str,
    actor_id: Optional[UUID] = None,
) -> Dict[str, Any]:
    """
    Execute bulk archive (soft delete) operation.
    
    Sets lifecycle_status to 'archived' for all entities.
    """
    return bulk_update_status(
        entity_type=entity_type,
        entity_ids=entity_ids,
        new_status=LifecycleStatus.ARCHIVED.value,
        reason=reason,
        actor_id=actor_id,
    )


def execute_bulk_delete(
    entity_type: str,
    entity_ids: List[UUID],
    reason: str,
    actor_id: Optional[UUID] = None,
    dry_run: bool = False,
) -> Dict[str, Any]:
    """
    Execute permanent deletion of entities.
    
    Args:
        dry_run: If True, return impact preview without deleting
    """
    results = {
        "total": len(entity_ids),
        "deleted": 0,
        "failed": 0,
        "errors": [],
        "dry_run": dry_run,
    }
    
    if dry_run:
        # Return aggregated impact preview
        return bulk_delete_preview(entity_type, entity_ids)
    
    model = ENTITY_MODELS.get(entity_type)
    if not model:
        results["errors"].append({"error": f"Unknown entity type: {entity_type}"})
        return results
    
    with db_session() as db:
        for entity_id in entity_ids:
            try:
                entity = db.query(model).filter(model.id == entity_id).first()
                if not entity:
                    results["failed"] += 1
                    results["errors"].append({"id": str(entity_id), "error": "Not found"})
                    continue
                
                # Get entity name for audit
                entity_name = getattr(entity, "name", None) or str(entity_id)
                
                # Cascade: Remove from association tables first
                _cascade_delete_associations(db, entity_type, entity_id)
                
                # Delete main entity
                db.delete(entity)
                
                # Audit log entry
                audit_entry = AuditLog(
                    entity_type=entity_type,
                    entity_id=entity_id,
                    action="permanent_delete",
                    old_value={"name": entity_name},
                    new_value={"reason": reason, "deleted_by": str(actor_id) if actor_id else None},
                    user_id=actor_id,
                )
                db.add(audit_entry)
                
                results["deleted"] += 1
                
            except Exception as e:
                results["failed"] += 1
                results["errors"].append({"id": str(entity_id), "error": str(e)})
        
        db.commit()
    
    logger.info(f"Bulk delete: {results['deleted']}/{results['total']} {entity_type}s deleted")
    return results


def _cascade_delete_associations(db: Session, entity_type: str, entity_id: UUID) -> None:
    """Remove entity from all association tables before deletion."""
    if entity_type == "dataset":
        db.execute(dataset_feature_assoc.delete().where(dataset_feature_assoc.c.dataset_id == entity_id))
        db.execute(dataset_signature_assoc.delete().where(dataset_signature_assoc.c.dataset_id == entity_id))
        db.execute(experiment_dataset_assoc.delete().where(experiment_dataset_assoc.c.experiment_id == entity_id))
    elif entity_type == "experiment":
        db.execute(experiment_dataset_assoc.delete().where(experiment_dataset_assoc.c.experiment_id == entity_id))
    elif entity_type == "signature":
        db.execute(dataset_signature_assoc.delete().where(dataset_signature_assoc.c.signature_id == entity_id))
        db.execute(signature_feature_assoc.delete().where(signature_feature_assoc.c.signature_id == entity_id))
    # Compound associations handled by FK cascade or separate logic


def find_orphaned_entities() -> Dict[str, Any]:
    """
    Find orphaned entities not linked to any parent.
    
    Returns counts and IDs for each orphan type.
    """
    orphans = {
        "features": {"count": 0, "ids": []},
        "signatures": {"count": 0, "ids": []},
        "embeddings": {"count": 0, "ids": []},
    }
    
    with db_session() as db:
        from amprenta_rag.database.models import Feature
        
        # Orphaned features (not linked to any dataset)
        try:
            from amprenta_rag.database.models import Feature, dataset_feature_assoc
            orphan_features = db.execute(
                select(Feature.id).outerjoin(
                    dataset_feature_assoc, Feature.id == dataset_feature_assoc.c.feature_id
                ).where(dataset_feature_assoc.c.dataset_id.is_(None))
            ).fetchall()
            orphans["features"]["count"] = len(orphan_features)
            orphans["features"]["ids"] = [str(row[0]) for row in orphan_features[:100]]
        except Exception as e:
            orphans["features"]["count"] = 0
            orphans["features"]["error"] = str(e)
        
        # Orphaned signatures (not linked to any dataset)
        try:
            from amprenta_rag.database.models import Signature, dataset_signature_assoc
            orphan_sigs = db.execute(
                select(Signature.id).outerjoin(
                    dataset_signature_assoc, Signature.id == dataset_signature_assoc.c.signature_id
                ).where(dataset_signature_assoc.c.dataset_id.is_(None))
            ).fetchall()
            orphans["signatures"]["count"] = len(orphan_sigs)
            orphans["signatures"]["ids"] = [str(row[0]) for row in orphan_sigs[:100]]
        except Exception as e:
            orphans["signatures"]["count"] = 0
            orphans["signatures"]["error"] = str(e)
        
        # Embeddings check (optional)
        orphans["embeddings"]["count"] = 0
        orphans["embeddings"]["note"] = "Embeddings check skipped (table schema pending)"
    
    return orphans


def cleanup_orphans(dry_run: bool = True) -> Dict[str, Any]:
    """
    Clean up orphaned entities.
    
    Args:
        dry_run: If True, only report what would be deleted
    """
    orphans = find_orphaned_entities()
    
    results = {
        "dry_run": dry_run,
        "deleted": {"features": 0, "signatures": 0, "embeddings": 0},
        "errors": [],
    }
    
    if dry_run:
        results["would_delete"] = {
            "features": orphans["features"]["count"],
            "signatures": orphans["signatures"]["count"],
            "embeddings": orphans["embeddings"]["count"],
        }
        return results
    
    with db_session() as db:
        from amprenta_rag.database.models import Feature
        
        # Delete orphaned features
        if orphans["features"]["ids"]:
            from uuid import UUID
            feature_ids = [UUID(fid) for fid in orphans["features"]["ids"]]
            deleted = db.query(Feature).filter(Feature.id.in_(feature_ids)).delete(synchronize_session=False)
            results["deleted"]["features"] = deleted
        
        # Delete orphaned signatures
        if orphans["signatures"]["ids"]:
            from uuid import UUID
            sig_ids = [UUID(sid) for sid in orphans["signatures"]["ids"]]
            deleted = db.query(Signature).filter(Signature.id.in_(sig_ids)).delete(synchronize_session=False)
            results["deleted"]["signatures"] = deleted
        
        db.commit()
    
    logger.info(f"Orphan cleanup: {results['deleted']}")
    return results


def enforce_retention_policies(dry_run: bool = True) -> Dict[str, Any]:
    """
    Enforce active retention policies.
    
    Finds entities older than policy retention_days and applies action.
    """
    results = {
        "dry_run": dry_run,
        "policies_checked": 0,
        "entities_affected": 0,
        "actions": [],
        "errors": [],
    }
    
    try:
        from amprenta_rag.models.misc import RetentionPolicy
        from datetime import datetime, timedelta, timezone
        
        with db_session() as db:
            # Get active policies
            policies = db.query(RetentionPolicy).filter(
                RetentionPolicy.is_active.is_(True)
            ).all()
            
            results["policies_checked"] = len(policies)
            
            for policy in policies:
                model = ENTITY_MODELS.get(policy.entity_type)
                if not model:
                    results["errors"].append(f"Unknown entity type: {policy.entity_type}")
                    continue
                
                cutoff = datetime.now(timezone.utc) - timedelta(days=policy.retention_days)
                
                # Find entities older than cutoff, not exempt
                query = db.query(model).filter(
                    model.created_at < cutoff,
                    model.lifecycle_status == "active",  # Only active entities
                )
                
                # Add retention_exempt filter if column exists
                if hasattr(model, 'retention_exempt'):
                    query = query.filter(model.retention_exempt.is_(False))
                
                affected_entities = query.all()
                
                action_record = {
                    "policy": policy.name,
                    "entity_type": policy.entity_type,
                    "action": policy.action,
                    "count": len(affected_entities),
                    "entity_ids": [str(e.id) for e in affected_entities[:10]],  # Sample
                }
                results["actions"].append(action_record)
                results["entities_affected"] += len(affected_entities)
                
                if dry_run:
                    continue
                
                # Execute action
                for entity in affected_entities:
                    try:
                        if policy.action == "archive":
                            update_lifecycle_status(
                                policy.entity_type, entity.id, "archived",
                                reason=f"Auto-archived by policy: {policy.name}"
                            )
                        elif policy.action == "delete":
                            # Use bulk delete for permanent removal
                            execute_bulk_delete(
                                policy.entity_type, [entity.id],
                                reason=f"Auto-deleted by policy: {policy.name}",
                                dry_run=False
                            )
                        elif policy.action == "notify":
                            # Log notification (actual notification system integration TBD)
                            logger.info(f"Retention notify: {policy.entity_type}:{entity.id}")
                    except Exception as e:
                        results["errors"].append(f"{entity.id}: {str(e)}")
            
            if not dry_run:
                db.commit()
    
    except ImportError:
        results["errors"].append("RetentionPolicy model not found - skipping enforcement")
    except Exception as e:
        results["errors"].append(f"Retention enforcement failed: {str(e)}")
    
    logger.info(f"Retention enforcement: {results['entities_affected']} entities affected")
    return results


def export_entity_data(
    entity_type: str,
    entity_id: UUID,
) -> Dict[str, Any]:
    """
    Export entity and related data for GDPR compliance.
    
    Returns dict with:
    - data: serialized entity + relationships
    - audit_trail: all audit log entries for this entity
    - checksum: SHA256 of the data
    - format: JSON
    """
    result = {
        "entity_type": entity_type,
        "entity_id": str(entity_id),
        "data": None,
        "audit_trail": [],
        "related_entities": {},
        "checksum": None,
        "exported_at": datetime.now(timezone.utc).isoformat(),
    }
    
    model = ENTITY_MODELS.get(entity_type)
    if not model:
        result["error"] = f"Unknown entity type: {entity_type}"
        return result
    
    with db_session() as db:
        entity = db.query(model).filter(model.id == entity_id).first()
        if not entity:
            result["error"] = "Entity not found"
            return result
        
        # Serialize main entity (exclude SQLAlchemy internals)
        entity_data = {}
        for column in model.__table__.columns:
            value = getattr(entity, column.name, None)
            if isinstance(value, (datetime,)):
                value = value.isoformat()
            elif isinstance(value, UUID):
                value = str(value)
            entity_data[column.name] = value
        
        result["data"] = entity_data
        
        # Get related entities based on type
        if entity_type == "dataset":
            # Features linked to dataset
            result["related_entities"]["features"] = [
                str(f.id) for f in entity.features
            ] if hasattr(entity, 'features') else []
            # Signatures linked to dataset
            result["related_entities"]["signatures"] = [
                str(s.id) for s in entity.signatures
            ] if hasattr(entity, 'signatures') else []
        
        elif entity_type == "experiment":
            result["related_entities"]["datasets"] = [
                str(d.id) for d in entity.datasets
            ] if hasattr(entity, 'datasets') else []
        
        # Get audit trail
        audit_entries = db.query(AuditLog).filter(
            AuditLog.entity_type == entity_type,
            AuditLog.entity_id == entity_id,
        ).order_by(AuditLog.timestamp.desc()).all()
        
        result["audit_trail"] = [
            {
                "timestamp": e.timestamp.isoformat() if e.timestamp else None,
                "action": e.action,
                "old_value": e.old_value,
                "new_value": e.new_value,
                "user_id": str(e.user_id) if e.user_id else None,
            }
            for e in audit_entries
        ]
    
    # Calculate checksum
    data_str = json.dumps(result["data"], sort_keys=True, default=str)
    result["checksum"] = hashlib.sha256(data_str.encode()).hexdigest()
    
    return result


def create_export_package(
    entity_type: str,
    entity_id: UUID,
) -> bytes:
    """
    Create a downloadable ZIP package with entity export.
    
    Returns ZIP file bytes.
    """
    export_data = export_entity_data(entity_type, entity_id)
    
    # Create ZIP in memory
    zip_buffer = BytesIO()
    with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zf:
        # Main data as JSON
        zf.writestr(
            f"{entity_type}_{entity_id}_data.json",
            json.dumps(export_data["data"], indent=2, default=str)
        )
        
        # Audit trail
        zf.writestr(
            f"{entity_type}_{entity_id}_audit.json",
            json.dumps(export_data["audit_trail"], indent=2, default=str)
        )
        
        # Manifest with checksum
        manifest = {
            "entity_type": entity_type,
            "entity_id": str(entity_id),
            "exported_at": export_data["exported_at"],
            "checksum": export_data["checksum"],
            "files": [
                f"{entity_type}_{entity_id}_data.json",
                f"{entity_type}_{entity_id}_audit.json",
            ]
        }
        zf.writestr("manifest.json", json.dumps(manifest, indent=2))
    
    zip_buffer.seek(0)
    return zip_buffer.getvalue()
