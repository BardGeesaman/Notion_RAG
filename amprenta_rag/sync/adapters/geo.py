"""GEO sync adapter for incremental harvesting."""

from __future__ import annotations

import hashlib
import json
import time
from datetime import datetime
from typing import Any, AsyncIterator, Dict, Optional
from uuid import UUID

from amprenta_rag.ingestion.repositories.geo import GEORepository
from amprenta_rag.sync.adapters.base import BaseSyncAdapter
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


class GEOSyncAdapter(BaseSyncAdapter):
    """Sync adapter for Gene Expression Omnibus (GEO)."""
    
    source = "geo"
    PAGE_SIZE = 50  # Records per page
    MAX_PAGES = 20  # Max 1000 total records per sync
    PAGE_RATE_LIMIT = 0.5  # P1 FIX: 0.5s between pages
    
    def __init__(self, geo_repo: Optional[GEORepository] = None):
        """Initialize with GEORepository for rate limiting."""
        self._repo = geo_repo or GEORepository()
    
    async def fetch_records(
        self, 
        since: Optional[datetime] = None
    ) -> AsyncIterator[Dict[str, Any]]:
        """Fetch GEO studies with pagination and rate limiting."""
        query = '"Homo sapiens"[Organism] AND gse[Entry Type]'
        
        if since:
            date_str = since.strftime('%Y/%m/%d')
            query += f" AND {date_str}:3000[MDAT]"
            logger.info(f"[GEO-SYNC] Incremental sync since {date_str}")
        else:
            logger.info("[GEO-SYNC] Full sync (no since date)")
        
        total_yielded = 0
        
        for page in range(self.MAX_PAGES):
            # P1 FIX: Rate limit between pages
            if page > 0:
                time.sleep(self.PAGE_RATE_LIMIT)
            
            retstart = page * self.PAGE_SIZE
            gse_ids = self._repo._search_geo(
                query, 
                max_results=self.PAGE_SIZE,
                retstart=retstart  # Add pagination parameter
            )
            
            if not gse_ids:
                logger.info(f"[GEO-SYNC] No more results at page {page}")
                break
            
            logger.info(f"[GEO-SYNC] Page {page+1}: Processing {len(gse_ids)} studies")
            
            for gse_id in gse_ids:
                try:
                    metadata = self._repo.fetch_study_metadata(gse_id)
                    if metadata:
                        record = self._to_sync_record(gse_id, metadata)
                        yield record
                        total_yielded += 1
                except Exception as e:
                    logger.error(f"[GEO-SYNC] Error fetching {gse_id}: {e}")
                    continue
            
            # Stop if we got fewer results than page size (last page)
            if len(gse_ids) < self.PAGE_SIZE:
                break
        
        logger.info(f"[GEO-SYNC] Completed: {total_yielded} studies synced")
    
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
