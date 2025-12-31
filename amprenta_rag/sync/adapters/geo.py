"""GEO sync adapter for incremental harvesting."""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from typing import Any, AsyncIterator, Dict, Optional
from uuid import UUID

from amprenta_rag.ingestion.repositories.geo import GEORepository
from amprenta_rag.sync.adapters.base import BaseSyncAdapter
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


class GEOSyncAdapter(BaseSyncAdapter):
    """Sync adapter for Gene Expression Omnibus (GEO)."""
    
    source = "geo"
    BATCH_SIZE = 50  # Smaller than ChEMBL due to heavier metadata
    
    def __init__(self, geo_repo: Optional[GEORepository] = None):
        """Initialize with GEORepository for rate limiting."""
        self._repo = geo_repo or GEORepository()
    
    async def fetch_records(
        self, 
        since: Optional[datetime] = None
    ) -> AsyncIterator[Dict[str, Any]]:
        """
        Fetch GEO studies modified since given date.
        
        Uses MDAT (modification date) filter for true incremental sync.
        """
        # Build query with MDAT filter for incremental
        query = '"Homo sapiens"[Organism] AND gse[Entry Type]'
        
        if since:
            # MDAT = modification date (includes updates, not just new)
            date_str = since.strftime('%Y/%m/%d')
            query += f" AND {date_str}:3000[MDAT]"
            logger.info(f"[GEO-SYNC] Incremental sync since {date_str}")
        else:
            logger.info("[GEO-SYNC] Full sync (no since date)")
        
        # Use GEORepository's search with rate limiting
        gse_ids = self._repo._search_geo(query, max_results=self.BATCH_SIZE)
        logger.info(f"[GEO-SYNC] Found {len(gse_ids)} studies to sync")
        
        for gse_id in gse_ids:
            try:
                metadata = self._repo.fetch_study_metadata(gse_id)
                if metadata:
                    record = self._to_sync_record(gse_id, metadata)
                    yield record
            except Exception as e:
                logger.error(f"[GEO-SYNC] Error fetching {gse_id}: {e}")
                continue
    
    def compute_checksum(self, record: Dict[str, Any]) -> str:
        """Compute MD5 hash for change detection."""
        # Use the data field for checksum computation
        data = record.get("data", {})
        json_str = json.dumps(data, sort_keys=True, default=str)
        return hashlib.md5(json_str.encode()).hexdigest()
    
    def map_to_entity(self, record: Dict[str, Any], db_session) -> tuple[str, UUID | None]:
        """Map external record to local entity. Returns (entity_type, entity_id)."""
        # For GEO studies, we map to a study entity
        # This would need to be implemented based on the actual entity model
        # For now, return the entity type and None for entity_id (new entity)
        return "study", None
    
    def _to_sync_record(self, gse_id: str, metadata) -> Dict[str, Any]:
        """Convert StudyMetadata to sync record format."""
        # Compute checksum for change detection
        raw_data = getattr(metadata, 'raw_metadata', None) or {}
        checksum = self._compute_checksum_internal(raw_data)
        
        return {
            "external_id": gse_id,
            "checksum": checksum,
            "data": {
                "study_id": gse_id,
                "title": getattr(metadata, 'title', ''),
                "summary": getattr(metadata, 'summary', ''),
                "organism": self._normalize_organism(getattr(metadata, 'organism', [])),
                "platform": getattr(metadata, 'platform', ''),
                "num_samples": getattr(metadata, 'num_samples', 0),
                "omics_type": getattr(metadata, 'omics_type', ''),
                "disease": getattr(metadata, 'disease', ''),
                "repository": "GEO",
            },
            "entity_type": "study",
        }
    
    def _compute_checksum_internal(self, data: dict) -> str:
        """Compute SHA256 checksum for change detection."""
        json_str = json.dumps(data, sort_keys=True, default=str)
        return hashlib.sha256(json_str.encode()).hexdigest()[:16]
    
    def _normalize_organism(self, organisms) -> list:
        """Normalize organism names (inline, minimal)."""
        # Handle both list and single string cases
        if isinstance(organisms, str):
            organisms = [organisms]
        elif not isinstance(organisms, list):
            organisms = []
            
        mapping = {
            "homo sapiens": "human",
            "mus musculus": "mouse", 
            "rattus norvegicus": "rat",
            "danio rerio": "zebrafish",
            "drosophila melanogaster": "fruit fly",
            "caenorhabditis elegans": "c. elegans",
            "saccharomyces cerevisiae": "yeast",
        }
        normalized = []
        for org in organisms:
            org_lower = str(org).lower().strip()
            normalized.append(mapping.get(org_lower, org))
        return normalized
