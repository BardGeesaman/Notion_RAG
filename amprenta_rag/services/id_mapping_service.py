"""
ID Mapping service for managing cross-reference mappings between databases.

This service handles:
- Getting cached ID mappings from the database
- Batch lookups for performance
- Saving new mappings with TTL support
- Refreshing mappings from external sources
- Migration from legacy gene_protein_map table
"""

import os
from datetime import datetime, timedelta, timezone
from typing import Dict, List, Optional

from sqlalchemy import func, text

from amprenta_rag.database.models import IDMapping
from amprenta_rag.database.session import db_session
from amprenta_rag.logging_utils import get_logger
from amprenta_rag.sync.adapters.uniprot_mapping import sync_uniprot_mappings

logger = get_logger(__name__)

# Settings (read from environment)
ID_MAPPING_FALLBACK_ENABLED: bool = os.getenv("ID_MAPPING_FALLBACK_ENABLED", "true").lower() in ("true", "1", "yes")
ID_MAPPING_CACHE_TTL_DAYS: int = int(os.getenv("ID_MAPPING_CACHE_TTL_DAYS", "90"))


def get_mapping(
    source_type: str,
    source_id: str,
    target_type: str,
    organism: str = "human",
    fallback: bool = True
) -> Optional[str]:
    """
    Get a single ID mapping from the cache.
    
    Args:
        source_type: Type of source ID (gene, protein, ensembl, etc.)
        source_id: Source identifier
        target_type: Type of target ID (uniprot, kegg_gene, etc.)
        organism: Target organism (default: human)
        fallback: Whether to use fallback strategies if not found
        
    Returns:
        Target ID if found, None otherwise
    """
    try:
        with db_session() as db:
            # Query for exact mapping
            mapping = db.query(IDMapping).filter(
                IDMapping.source_type == source_type,
                IDMapping.source_id == source_id,
                IDMapping.target_type == target_type,
                IDMapping.organism == organism
            ).first()
            
            if mapping:
                # Check if mapping is expired
                if mapping.expires_at and mapping.expires_at < datetime.now(timezone.utc):
                    logger.debug(f"Mapping {source_type}:{source_id} -> {target_type} expired")
                    return None
                
                return mapping.target_id
            
            # Fallback strategies if enabled
            if fallback and ID_MAPPING_FALLBACK_ENABLED:
                # Try case-insensitive match
                mapping = db.query(IDMapping).filter(
                    IDMapping.source_type == source_type,
                    func.lower(IDMapping.source_id) == source_id.lower(),
                    IDMapping.target_type == target_type,
                    IDMapping.organism == organism
                ).first()
                
                if mapping and (not mapping.expires_at or mapping.expires_at >= datetime.now(timezone.utc)):
                    logger.debug(f"Found case-insensitive mapping for {source_id}")
                    return mapping.target_id
            
            return None
            
    except Exception as e:
        logger.error(f"Error getting mapping {source_type}:{source_id} -> {target_type}: {e}")
        return None


def get_mappings_batch(
    ids: List[str],
    source_type: str,
    target_type: str,
    organism: str = "human"
) -> Dict[str, Optional[str]]:
    """
    Get multiple ID mappings in a single database query.
    
    Args:
        ids: List of source identifiers
        source_type: Type of source IDs
        target_type: Type of target IDs
        organism: Target organism
        
    Returns:
        Dict mapping source_id -> target_id (or None if not found)
    """
    result = {id_: None for id_ in ids}
    
    if not ids:
        return result
    
    try:
        with db_session() as db:
            # Query all mappings at once
            mappings = db.query(IDMapping).filter(
                IDMapping.source_type == source_type,
                IDMapping.source_id.in_(ids),
                IDMapping.target_type == target_type,
                IDMapping.organism == organism
            ).all()
            
            now = datetime.now(timezone.utc)
            
            for mapping in mappings:
                # Check expiration
                if mapping.expires_at and mapping.expires_at < now:
                    continue
                
                result[mapping.source_id] = mapping.target_id
            
            return result
            
    except Exception as e:
        logger.error(f"Error getting batch mappings {source_type} -> {target_type}: {e}")
        return result


def save_mapping(
    source_type: str,
    source_id: str,
    target_type: str,
    target_id: str,
    organism: str = "human",
    ttl_days: Optional[int] = None,
    confidence: Optional[float] = None
) -> Optional[IDMapping]:
    """
    Save or update an ID mapping.
    
    Args:
        source_type: Type of source ID
        source_id: Source identifier
        target_type: Type of target ID
        target_id: Target identifier
        organism: Target organism
        ttl_days: Time-to-live in days (None for permanent)
        confidence: Confidence score (0.0-1.0)
        
    Returns:
        Created or updated IDMapping, or None on error
    """
    try:
        with db_session() as db:
            # Check for existing mapping
            existing = db.query(IDMapping).filter(
                IDMapping.source_type == source_type,
                IDMapping.source_id == source_id,
                IDMapping.target_type == target_type,
                IDMapping.organism == organism
            ).first()
            
            expires_at = None
            if ttl_days is not None:
                expires_at = datetime.now(timezone.utc) + timedelta(days=ttl_days)
            
            if existing:
                # Update existing mapping
                existing.target_id = target_id
                existing.confidence = confidence
                existing.expires_at = expires_at
                existing.updated_at = datetime.now(timezone.utc)
                mapping = existing
            else:
                # Create new mapping
                mapping = IDMapping(
                    source_type=source_type,
                    source_id=source_id,
                    target_type=target_type,
                    target_id=target_id,
                    organism=organism,
                    confidence=confidence,
                    expires_at=expires_at
                )
                db.add(mapping)
            
            db.commit()
            db.refresh(mapping)
            return mapping
            
    except Exception as e:
        logger.error(f"Error saving mapping {source_type}:{source_id} -> {target_type}:{target_id}: {e}")
        return None


async def refresh_uniprot_mappings() -> int:
    """
    Refresh UniProt ID mappings from external source.
    
    Returns:
        Number of mappings upserted
    """
    try:
        logger.info("Starting UniProt mappings refresh...")
        count = await sync_uniprot_mappings()
        logger.info(f"UniProt mappings refresh completed. Upserted {count} mappings")
        return count
        
    except Exception as e:
        logger.error(f"Error refreshing UniProt mappings: {e}")
        return 0


def get_mapping_stats() -> Dict[str, any]:
    """
    Get statistics about ID mappings coverage and freshness.
    
    Returns:
        Dict with statistics
    """
    try:
        with db_session() as db:
            # Total mappings
            total_mappings = db.query(IDMapping).count()
            
            # Mappings by source type
            source_stats = db.query(
                IDMapping.source_type,
                func.count(IDMapping.id).label('count')
            ).group_by(IDMapping.source_type).all()
            
            # Mappings by target type
            target_stats = db.query(
                IDMapping.target_type,
                func.count(IDMapping.id).label('count')
            ).group_by(IDMapping.target_type).all()
            
            # Expired mappings
            now = datetime.now(timezone.utc)
            expired_count = db.query(IDMapping).filter(
                IDMapping.expires_at < now
            ).count()
            
            # Permanent mappings (no expiration)
            permanent_count = db.query(IDMapping).filter(
                IDMapping.expires_at.is_(None)
            ).count()
            
            return {
                "total_mappings": total_mappings,
                "permanent_mappings": permanent_count,
                "expired_mappings": expired_count,
                "by_source_type": {row.source_type: row.count for row in source_stats},
                "by_target_type": {row.target_type: row.count for row in target_stats},
                "last_updated": datetime.now(timezone.utc).isoformat()
            }
            
    except Exception as e:
        logger.error(f"Error getting mapping statistics: {e}")
        return {
            "total_mappings": 0,
            "error": str(e),
            "last_updated": datetime.now(timezone.utc).isoformat()
        }


def migrate_gene_protein_map() -> int:
    """
    One-time migration from legacy gene_protein_map table to IDMapping.
    
    Returns:
        Number of records migrated
    """
    migrated_count = 0
    
    try:
        with db_session() as db:
            # Query legacy table
            legacy_mappings = db.execute(
                text("SELECT gene_symbol, uniprot_id FROM gene_protein_map")
            ).fetchall()
            
            logger.info(f"Found {len(legacy_mappings)} records in gene_protein_map to migrate")
            
            for gene_symbol, uniprot_id in legacy_mappings:
                if not gene_symbol or not uniprot_id:
                    continue
                
                # Check if mapping already exists
                existing = db.query(IDMapping).filter(
                    IDMapping.source_type == "gene",
                    IDMapping.source_id == gene_symbol,
                    IDMapping.target_type == "uniprot",
                    IDMapping.organism == "human"
                ).first()
                
                if not existing:
                    # Create new mapping
                    mapping = IDMapping(
                        source_type="gene",
                        source_id=gene_symbol,
                        target_type="uniprot",
                        target_id=uniprot_id,
                        organism="human",
                        confidence=0.9,  # High confidence for curated mappings
                        expires_at=None  # Permanent
                    )
                    db.add(mapping)
                    migrated_count += 1
                
                # Also create reverse mapping
                existing_reverse = db.query(IDMapping).filter(
                    IDMapping.source_type == "uniprot",
                    IDMapping.source_id == uniprot_id,
                    IDMapping.target_type == "gene",
                    IDMapping.organism == "human"
                ).first()
                
                if not existing_reverse:
                    reverse_mapping = IDMapping(
                        source_type="uniprot",
                        source_id=uniprot_id,
                        target_type="gene",
                        target_id=gene_symbol,
                        organism="human",
                        confidence=0.9,
                        expires_at=None
                    )
                    db.add(reverse_mapping)
                    migrated_count += 1
            
            db.commit()
            logger.info(f"Migration completed. Migrated {migrated_count} mappings from gene_protein_map")
            
    except Exception as e:
        logger.error(f"Error migrating gene_protein_map: {e}")
        
    return migrated_count


def cleanup_expired_mappings() -> int:
    """
    Remove expired mappings from the database.
    
    Returns:
        Number of mappings removed
    """
    try:
        with db_session() as db:
            now = datetime.now(timezone.utc)
            
            # Delete expired mappings
            deleted_count = db.query(IDMapping).filter(
                IDMapping.expires_at < now
            ).delete()
            
            db.commit()
            
            if deleted_count > 0:
                logger.info(f"Cleaned up {deleted_count} expired mappings")
            
            return deleted_count
            
    except Exception as e:
        logger.error(f"Error cleaning up expired mappings: {e}")
        return 0
