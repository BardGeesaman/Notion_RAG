"""Data lifecycle management service.

Provides:
- Cascade impact calculation for deletions
- Status updates with audit trail
- Bulk operations with confirmation workflow
"""

from __future__ import annotations

from datetime import datetime, timezone
from typing import Any, Dict, List, Optional, Tuple
from uuid import UUID

from sqlalchemy import func
from sqlalchemy.orm import Session

from amprenta_rag.database.models import (
    Dataset, Experiment, Signature, AuditLog, LifecycleStatus,
    dataset_feature_assoc, dataset_signature_assoc, experiment_dataset_assoc,
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
