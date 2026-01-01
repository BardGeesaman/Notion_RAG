"""KEGG cache refresh adapter for expiring mappings."""

from datetime import datetime, timedelta, timezone
from typing import AsyncIterator, Dict, Any, Optional
from uuid import UUID
import time

from amprenta_rag.sync.adapters.base import BaseSyncAdapter
from amprenta_rag.database.models import IDMapping
from amprenta_rag.database.session import db_session
from amprenta_rag.analysis.id_mapping import (
    map_gene_to_kegg,
    map_metabolite_to_kegg,
)
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)

class KEGGRefreshAdapter(BaseSyncAdapter):
    """Adapter to refresh expiring KEGG mappings."""
    
    source = "kegg_refresh"
    REFRESH_WINDOW_DAYS = 7  # Refresh mappings expiring within 7 days
    RATE_LIMIT_SECONDS = 0.35  # ~3 req/sec max
    
    async def fetch_records(self, since: Optional[datetime] = None) -> AsyncIterator[Dict[str, Any]]:
        """Find mappings approaching expiration and refresh them."""
        cutoff = datetime.now(timezone.utc) + timedelta(days=self.REFRESH_WINDOW_DAYS)
        
        with db_session() as db:
            expiring = db.query(IDMapping).filter(
                IDMapping.target_type.like("kegg_%"),
                IDMapping.expires_at < cutoff,
                IDMapping.expires_at > datetime.now(timezone.utc)
            ).limit(100).all()  # Process in batches
            
            for mapping in expiring:
                db.expunge(mapping)
        
        logger.info(f"[KEGG-REFRESH] Found {len(expiring)} expiring mappings to refresh")
        
        for mapping in expiring:
            time.sleep(self.RATE_LIMIT_SECONDS)
            
            # Re-fetch from KEGG API
            new_target_id = None
            if mapping.target_type == "kegg_gene":
                new_target_id = map_gene_to_kegg(mapping.source_id)
            elif mapping.target_type == "kegg_compound":
                new_target_id = map_metabolite_to_kegg(mapping.source_id)
            
            yield {
                "external_id": f"{mapping.source_type}:{mapping.source_id}",
                "mapping_id": str(mapping.id),
                "source_type": mapping.source_type,
                "source_id": mapping.source_id,
                "target_type": mapping.target_type,
                "old_target_id": mapping.target_id,
                "new_target_id": new_target_id,
                "refreshed": new_target_id is not None,
            }
    
    def compute_checksum(self, record: Dict[str, Any]) -> str:
        import hashlib
        import json
        content = json.dumps({
            "source": record.get("source_id"),
            "target": record.get("new_target_id"),
        }, sort_keys=True)
        return hashlib.md5(content.encode()).hexdigest()
    
    def map_to_entity(self, record: Dict[str, Any], db_session) -> tuple[str, UUID | None]:
        return ("id_mapping", None)
