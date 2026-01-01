"""Auto-conflict resolution for sync jobs."""

from datetime import datetime, timezone
from enum import Enum
from typing import Any, Dict, Optional
from uuid import UUID

from amprenta_rag.database.models import SyncConflict
from amprenta_rag.database.session import db_session
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)

class ResolutionStrategy(str, Enum):
    PREFER_EXTERNAL = "prefer_external"
    PREFER_LOCAL = "prefer_local"
    PREFER_NEWEST = "prefer_newest"
    MANUAL_REQUIRED = "manual_required"

# Default strategies by conflict type
DEFAULT_STRATEGIES: Dict[str, ResolutionStrategy] = {
    "schema_mismatch": ResolutionStrategy.PREFER_EXTERNAL,
    "value_conflict": ResolutionStrategy.PREFER_NEWEST,
    "missing_field": ResolutionStrategy.PREFER_EXTERNAL,
    "type_mismatch": ResolutionStrategy.MANUAL_REQUIRED,
}

def resolve_conflict(
    conflict_id: UUID,
    strategy: Optional[ResolutionStrategy] = None,
) -> bool:
    """Resolve a sync conflict using the specified strategy."""
    with db_session() as db:
        conflict = db.query(SyncConflict).filter(SyncConflict.id == conflict_id).first()
        if not conflict:
            logger.warning(f"Conflict {conflict_id} not found")
            return False
        
        if conflict.resolution_status == "resolved":
            logger.debug(f"Conflict {conflict_id} already resolved")
            return True
        
        # Determine strategy
        strat = strategy or DEFAULT_STRATEGIES.get(
            conflict.conflict_type, ResolutionStrategy.MANUAL_REQUIRED
        )
        
        if strat == ResolutionStrategy.MANUAL_REQUIRED:
            logger.info(f"Conflict {conflict_id} requires manual resolution")
            return False
        
        # Apply resolution
        resolved_value = _apply_strategy(strat, conflict.local_value, conflict.external_value)
        
        conflict.resolution_status = "resolved"
        conflict.resolved_at = datetime.now(timezone.utc)
        conflict.resolution_strategy = strat.value
        conflict.resolved_value = resolved_value
        
        db.commit()
        logger.info(f"Resolved conflict {conflict_id} using {strat.value}")
        return True

def _apply_strategy(
    strategy: ResolutionStrategy,
    local_value: Dict[str, Any],
    external_value: Dict[str, Any],
) -> Dict[str, Any]:
    """Apply resolution strategy and return resolved value."""
    if strategy == ResolutionStrategy.PREFER_EXTERNAL:
        return external_value
    elif strategy == ResolutionStrategy.PREFER_LOCAL:
        return local_value
    elif strategy == ResolutionStrategy.PREFER_NEWEST:
        local_ts = local_value.get("updated_at") or local_value.get("synced_at") or ""
        external_ts = external_value.get("updated_at") or external_value.get("synced_at") or ""
        return external_value if external_ts >= local_ts else local_value
    else:
        return local_value  # Fallback

def auto_resolve_pending_conflicts(source: Optional[str] = None, limit: int = 100) -> int:
    """Auto-resolve pending conflicts for a source."""
    resolved_count = 0
    
    with db_session() as db:
        query = db.query(SyncConflict).filter(SyncConflict.resolution_status == "pending")
        if source:
            # Join through SyncRecord to filter by source
            from amprenta_rag.database.models import SyncRecord
            query = query.join(SyncRecord).filter(SyncRecord.source == source)
        
        conflicts = query.limit(limit).all()
        conflict_ids = [c.id for c in conflicts]
        for c in conflicts:
            db.expunge(c)
    
    for conflict_id in conflict_ids:
        if resolve_conflict(conflict_id):
            resolved_count += 1
    
    logger.info(f"Auto-resolved {resolved_count}/{len(conflict_ids)} conflicts")
    return resolved_count
