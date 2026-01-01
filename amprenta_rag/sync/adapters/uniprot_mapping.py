"""UniProt ID mapping synchronization adapter."""

from __future__ import annotations

import asyncio
import gzip
import logging
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import AsyncIterator

import httpx
from uuid import UUID

from amprenta_rag.database.models import IDMapping, MappingRefreshLog
from amprenta_rag.database.session import db_session
from amprenta_rag.sync.adapters.base import BaseSyncAdapter


logger = logging.getLogger(__name__)


class UniProtMappingAdapter(BaseSyncAdapter):
    """Adapter for syncing UniProt ID mappings."""
    
    source = "uniprot_mapping"
    
    # UniProt ID mapping file for human
    MAPPING_URL = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping_selected.tab.gz"
    
    # Column mappings for the tab file
    # Columns: UniProt-AC, UniProt-ID, GeneID (Entrez), RefSeq, GI, PDB, GO, UniRef100, UniRef90, UniRef50, UniParc, PIR, NCBI-taxon, MIM, UniGene, PubMed, EMBL, EMBL-CDS, Ensembl, Ensembl_TRS, Ensembl_PRO, Additional PubMed
    COLUMNS = [
        "uniprot_ac", "uniprot_id", "gene_id", "refseq", "gi", "pdb", "go", 
        "uniref100", "uniref90", "uniref50", "uniparc", "pir", "ncbi_taxon", 
        "mim", "unigene", "pubmed", "embl", "embl_cds", "ensembl", 
        "ensembl_trs", "ensembl_pro", "additional_pubmed"
    ]
    
    async def fetch_records(self, since: datetime | None = None) -> AsyncIterator[dict]:
        """Download and parse UniProt mapping file with conditional GET."""
        logger.info("Checking UniProt ID mapping file for updates...")
        
        # Check last successful sync for ETag/Last-Modified
        last_etag = None
        last_modified = None
        with db_session() as db:
            last_log = db.query(MappingRefreshLog).filter(
                MappingRefreshLog.source == "uniprot_mapping",
                MappingRefreshLog.status == "success"
            ).order_by(MappingRefreshLog.completed_at.desc()).first()
            if last_log and last_log.metadata_:
                last_etag = last_log.metadata_.get("etag")
                last_modified = last_log.metadata_.get("last_modified")
        
        # Build conditional headers
        headers = {}
        if last_etag:
            headers["If-None-Match"] = last_etag
        if last_modified:
            headers["If-Modified-Since"] = last_modified
        
        # Download the gzipped file
        async with httpx.AsyncClient(timeout=300.0) as client:
            response = await client.get(self.MAPPING_URL, headers=headers)
            
            # 304 Not Modified - file unchanged
            if response.status_code == 304:
                logger.info("UniProt mapping file unchanged (304 Not Modified)")
                return  # Empty iterator - nothing to sync
            
            response.raise_for_status()
            
            # Store ETag/Last-Modified for next sync
            self._last_etag = response.headers.get("ETag")
            self._last_modified = response.headers.get("Last-Modified")
            
            # Write to temporary file
            with tempfile.NamedTemporaryFile(suffix=".tab.gz", delete=False) as temp_file:
                temp_file.write(response.content)
                temp_path = Path(temp_file.name)
        
        try:
            logger.info(f"Parsing UniProt mapping file: {temp_path}")
            
            # Parse the gzipped tab file
            with gzip.open(temp_path, 'rt', encoding='utf-8') as f:
                line_count = 0
                for line in f:
                    line_count += 1
                    if line_count % 10000 == 0:
                        logger.info(f"Processed {line_count} lines...")
                    
                    # Parse tab-separated line
                    parts = line.strip().split('\t')
                    if len(parts) < 3:  # Need at least UniProt AC and one other ID
                        continue
                    
                    # Create mapping dict
                    mapping_data = {}
                    for i, value in enumerate(parts):
                        if i < len(self.COLUMNS) and value.strip():
                            mapping_data[self.COLUMNS[i]] = value.strip()
                    
                    if not mapping_data.get("uniprot_ac"):
                        continue  # Skip if no UniProt AC
                    
                    yield mapping_data
                    
            logger.info(f"Completed parsing {line_count} lines from UniProt mapping file")
            
        finally:
            # Clean up temporary file
            temp_path.unlink(missing_ok=True)
    
    async def transform_record(self, record: dict) -> list[dict]:
        """Transform UniProt mapping record into IDMapping records."""
        mappings = []
        uniprot_ac = record.get("uniprot_ac")
        
        if not uniprot_ac:
            return mappings
        
        # Map gene symbols (if available via gene_id -> gene symbol lookup)
        if record.get("gene_id"):
            mappings.append({
                "source_type": "gene",
                "source_id": record["gene_id"],
                "target_type": "uniprot",
                "target_id": uniprot_ac,
                "organism": "human",
                "confidence": 1.0,
                "expires_at": None,  # Permanent mapping
            })
        
        # Map Ensembl gene IDs
        if record.get("ensembl"):
            mappings.append({
                "source_type": "ensembl",
                "source_id": record["ensembl"],
                "target_type": "uniprot",
                "target_id": uniprot_ac,
                "organism": "human",
                "confidence": 1.0,
                "expires_at": None,
            })
        
        # Map RefSeq IDs
        if record.get("refseq"):
            mappings.append({
                "source_type": "refseq",
                "source_id": record["refseq"],
                "target_type": "uniprot",
                "target_id": uniprot_ac,
                "organism": "human",
                "confidence": 1.0,
                "expires_at": None,
            })
        
        # Reverse mappings (UniProt -> other IDs)
        if record.get("gene_id"):
            mappings.append({
                "source_type": "uniprot",
                "source_id": uniprot_ac,
                "target_type": "gene",
                "target_id": record["gene_id"],
                "organism": "human",
                "confidence": 1.0,
                "expires_at": None,
            })
        
        if record.get("ensembl"):
            mappings.append({
                "source_type": "uniprot",
                "source_id": uniprot_ac,
                "target_type": "ensembl",
                "target_id": record["ensembl"],
                "organism": "human",
                "confidence": 1.0,
                "expires_at": None,
            })
        
        return mappings
    
    async def save_records(self, records: list[dict]) -> int:
        """Save IDMapping records to database."""
        if not records:
            return 0
        
        saved_count = 0
        
        with db_session() as db:
            for record_data in records:
                try:
                    # Check if mapping already exists
                    existing = db.query(IDMapping).filter(
                        IDMapping.source_type == record_data["source_type"],
                        IDMapping.source_id == record_data["source_id"],
                        IDMapping.target_type == record_data["target_type"],
                        IDMapping.organism == record_data["organism"]
                    ).first()
                    
                    if existing:
                        # Update existing record
                        existing.target_id = record_data["target_id"]
                        existing.confidence = record_data.get("confidence")
                        existing.expires_at = record_data.get("expires_at")
                        existing.updated_at = datetime.now(timezone.utc)
                    else:
                        # Create new record
                        mapping = IDMapping(**record_data)
                        db.add(mapping)
                    
                    saved_count += 1
                    
                except Exception as e:
                    logger.warning(f"Failed to save mapping {record_data}: {e}")
                    continue
            
            db.commit()
        
        return saved_count
    
    async def sync_batch(self, records: list[dict]) -> int:
        """Process a batch of UniProt mapping records."""
        all_mappings = []
        
        for record in records:
            mappings = await self.transform_record(record)
            all_mappings.extend(mappings)
        
        return await self.save_records(all_mappings)
    
    def compute_checksum(self, record: dict) -> str:
        """Compute MD5 hash for change detection."""
        import hashlib
        import json
        
        # Create a stable representation for hashing
        stable_data = {
            "uniprot_ac": record.get("uniprot_ac", ""),
            "gene_id": record.get("gene_id", ""),
            "ensembl": record.get("ensembl", ""),
            "refseq": record.get("refseq", "")
        }
        
        content = json.dumps(stable_data, sort_keys=True)
        return hashlib.md5(content.encode()).hexdigest()
    
    def map_to_entity(self, record: dict, db_session) -> tuple[str, UUID | None]:
        """Map external record to local entity. Returns (entity_type, entity_id)."""
        # For ID mappings, we don't map to a specific entity
        # This is a cross-reference table, not entity-specific
        return ("id_mapping", None)
    
    def get_sync_metadata(self) -> dict:
        """Return metadata to store in MappingRefreshLog."""
        return {
            "etag": getattr(self, '_last_etag', None),
            "last_modified": getattr(self, '_last_modified', None),
        }


async def sync_uniprot_mappings() -> int:
    """Sync UniProt ID mappings."""
    adapter = UniProtMappingAdapter()
    total_saved = 0
    batch_size = 1000
    batch = []
    
    async for record in adapter.fetch_records():
        batch.append(record)
        
        if len(batch) >= batch_size:
            saved = await adapter.sync_batch(batch)
            total_saved += saved
            logger.info(f"Saved {saved} mappings in batch (total: {total_saved})")
            batch = []
    
    # Process final batch
    if batch:
        saved = await adapter.sync_batch(batch)
        total_saved += saved
        logger.info(f"Saved {saved} mappings in final batch (total: {total_saved})")
    
    logger.info(f"UniProt mapping sync complete. Total mappings saved: {total_saved}")
    return total_saved


if __name__ == "__main__":
    asyncio.run(sync_uniprot_mappings())
