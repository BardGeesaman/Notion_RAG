"""
Postgres-based signature metadata extraction.

Extracts signature metadata (short_id, biomarker_role, phenotype_axes, data_ownership)
from Postgres Signature models instead of Notion pages.
"""

from __future__ import annotations

from typing import Any, Dict, List
from uuid import UUID

from amprenta_rag.database.base import get_db
from amprenta_rag.database.models import Signature as SignatureModel
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def get_signature_metadata_from_postgres(signature: SignatureModel) -> Dict[str, Any]:
    """
    Extract metadata from a Postgres Signature model.
    
    Args:
        signature: SignatureModel instance from Postgres
        
    Returns:
        Dictionary with signature metadata:
        - short_id: str | None
        - biomarker_role: List[str]
        - phenotype_axes: List[str]
        - data_ownership: str | None
    """
    return {
        "short_id": signature.short_id,
        "biomarker_role": signature.biomarker_role or [],
        "phenotype_axes": signature.phenotype_axes or [],
        "data_ownership": signature.data_ownership,
    }


def collect_signature_metadata_from_postgres(
    signature_ids: List[UUID],
) -> Dict[str, Any]:
    """
    Fetch Signature models from Postgres and collect their metadata.
    
    Similar to the Notion-based collect_signature_metadata() but uses Postgres.
    
    Args:
        signature_ids: List of Postgres UUIDs for signature models
        
    Returns:
        Dictionary with aggregated metadata:
        - sig_short_ids: List[str] - Short IDs from signatures
        - sig_roles: List[str] - Biomarker roles (aggregated)
        - sig_axes: List[str] - Phenotype axes (aggregated)
        - sig_ownership: List[str] - Data ownership values (aggregated)
    """
    sig_short_ids: List[str] = []
    sig_roles: List[str] = []
    sig_axes: List[str] = []
    sig_ownership: List[str] = []
    
    if not signature_ids:
        return {
            "sig_short_ids": [],
            "sig_roles": [],
            "sig_axes": [],
            "sig_ownership": [],
        }
    
    db = next(get_db())
    
    try:
        signatures = (
            db.query(SignatureModel)
            .filter(SignatureModel.id.in_(signature_ids))
            .all()
        )
        
        for signature in signatures:
            # Extract short ID
            if signature.short_id:
                sig_short_ids.append(signature.short_id)
            elif signature.name:
                # Fallback to name if short_id not set
                sig_short_ids.append(signature.name)
            
            # Extract biomarker roles
            if signature.biomarker_role:
                sig_roles.extend(signature.biomarker_role)
            
            # Extract phenotype axes
            if signature.phenotype_axes:
                sig_axes.extend(signature.phenotype_axes)
            
            # Extract data ownership
            if signature.data_ownership:
                sig_ownership.append(signature.data_ownership)
        
        return {
            "sig_short_ids": list(set(sig_short_ids)),  # Deduplicate
            "sig_roles": list(set(sig_roles)),  # Deduplicate
            "sig_axes": list(set(sig_axes)),  # Deduplicate
            "sig_ownership": list(set(sig_ownership)),  # Deduplicate
        }
    except Exception as e:
        logger.error(
            "[METADATA] Error fetching signature metadata from Postgres: %r",
            e,
        )
        return {
            "sig_short_ids": [],
            "sig_roles": [],
            "sig_axes": [],
            "sig_ownership": [],
        }
    finally:
        db.close()


def find_signatures_by_short_id(short_id: str) -> List[SignatureModel]:
    """
    Find signatures in Postgres by short ID.
    
    Args:
        short_id: Short identifier to search for
        
    Returns:
        List of SignatureModel instances matching the short_id
    """
    db = next(get_db())
    try:
        signatures = (
            db.query(SignatureModel)
            .filter(SignatureModel.short_id == short_id)
            .all()
        )
        return signatures
    finally:
        db.close()


def find_signature_by_notion_id(notion_page_id: str) -> SignatureModel | None:
    """
    Find a signature in Postgres by Notion page ID.
    
    Args:
        notion_page_id: Notion page ID (with or without dashes)
        
    Returns:
        SignatureModel instance or None if not found
    """
    db = next(get_db())
    try:
        # Clean up the notion_page_id (remove dashes for comparison)
        clean_id = notion_page_id.replace("-", "")
        
        signature = (
            db.query(SignatureModel)
            .filter(SignatureModel.notion_page_id == clean_id)
            .first()
        )
        return signature
    finally:
        db.close()

