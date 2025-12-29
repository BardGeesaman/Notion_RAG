"""
IP and patent tracking service layer.

Provides functions for managing invention disclosures, patents, and entity links.
"""

from __future__ import annotations

from datetime import datetime
from typing import List, Optional
from uuid import UUID

from sqlalchemy.orm import Session

from amprenta_rag.models.ip import (
    DisclosureInventor,
    InventionDisclosure,
    IPLink,
    PatentApplication,
    PatentClaim,
)
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)

# Valid status transitions
STATUS_TRANSITIONS = {
    "draft": ["submitted"],
    "submitted": ["under_review"],
    "under_review": ["filed", "rejected"],
    "filed": ["granted", "abandoned"],
}


def create_disclosure(
    title: str,
    description: str,
    technical_field: Optional[str],
    inventors: List[dict],
    user_id: UUID,
    company_id: UUID,
    db: Session,
) -> InventionDisclosure:
    """
    Create new invention disclosure with inventors.

    Args:
        title: Invention title
        description: Technical description
        technical_field: Field of invention
        inventors: List of dicts with user_id, contribution_percentage, is_primary
        user_id: Creating user UUID
        company_id: Company UUID
        db: Database session

    Returns:
        Created InventionDisclosure
    """
    disclosure = InventionDisclosure(
        title=title,
        description=description,
        technical_field=technical_field,
        status="draft",
        created_by_id=user_id,
        company_id=company_id,
    )
    
    db.add(disclosure)
    db.flush()  # Get disclosure.id
    
    # Add inventors
    for inventor_data in inventors:
        inventor = DisclosureInventor(
            disclosure_id=disclosure.id,
            user_id=inventor_data["user_id"],
            contribution_percentage=inventor_data.get("contribution_percentage"),
            is_primary=inventor_data.get("is_primary", False),
        )
        db.add(inventor)
    
    db.commit()
    db.refresh(disclosure)
    
    logger.info("[IP] Created disclosure: %s", title)
    
    return disclosure


def update_disclosure_status(
    disclosure_id: UUID,
    new_status: str,
    user_id: UUID,
    db: Session,
) -> InventionDisclosure:
    """
    Update disclosure status with workflow validation.

    Args:
        disclosure_id: Disclosure UUID
        new_status: Target status
        user_id: User making the change
        db: Database session

    Returns:
        Updated InventionDisclosure

    Raises:
        ValueError: If transition is invalid
    """
    disclosure = db.query(InventionDisclosure).filter(InventionDisclosure.id == disclosure_id).first()
    
    if not disclosure:
        raise ValueError(f"Disclosure {disclosure_id} not found")
    
    # Validate transition
    current_status = disclosure.status
    allowed_transitions = STATUS_TRANSITIONS.get(current_status, [])
    
    if new_status not in allowed_transitions:
        raise ValueError(
            f"Invalid status transition: {current_status} -> {new_status}. "
            f"Allowed: {allowed_transitions}"
        )
    
    disclosure.status = new_status
    db.commit()
    db.refresh(disclosure)
    
    logger.info("[IP] Updated disclosure %s status: %s -> %s", disclosure_id, current_status, new_status)
    
    return disclosure


def file_patent(
    disclosure_id: UUID,
    application_number: str,
    jurisdiction: str,
    filing_date: datetime,
    db: Session,
) -> PatentApplication:
    """
    File patent application from disclosure.

    Args:
        disclosure_id: Disclosure UUID
        application_number: Patent application number
        jurisdiction: Jurisdiction code (US, EP, PCT, etc.)
        filing_date: Filing date
        db: Database session

    Returns:
        Created PatentApplication
    """
    disclosure = db.query(InventionDisclosure).filter(InventionDisclosure.id == disclosure_id).first()
    
    if not disclosure:
        raise ValueError(f"Disclosure {disclosure_id} not found")
    
    patent = PatentApplication(
        disclosure_id=disclosure_id,
        application_number=application_number,
        filing_date=filing_date,
        jurisdiction=jurisdiction,
        status="pending",
    )
    
    db.add(patent)
    db.commit()
    db.refresh(patent)
    
    logger.info("[IP] Filed patent %s for disclosure %s", application_number, disclosure_id)
    
    return patent


def add_patent_claim(
    patent_id: UUID,
    claim_number: int,
    claim_text: str,
    claim_type: str,
    parent_claim_id: Optional[UUID],
    db: Session,
) -> PatentClaim:
    """
    Add claim to patent application.

    Args:
        patent_id: Patent UUID
        claim_number: Claim number
        claim_text: Claim text
        claim_type: independent or dependent
        parent_claim_id: Parent claim UUID (for dependent claims)
        db: Database session

    Returns:
        Created PatentClaim
    """
    claim = PatentClaim(
        patent_id=patent_id,
        claim_number=claim_number,
        claim_text=claim_text,
        claim_type=claim_type,
        parent_claim_id=parent_claim_id,
    )
    
    db.add(claim)
    db.commit()
    db.refresh(claim)
    
    logger.info("[IP] Added claim %d to patent %s", claim_number, patent_id)
    
    return claim


def link_entity_to_disclosure(
    entity_type: str,
    entity_id: UUID,
    disclosure_id: UUID,
    link_type: str,
    notes: Optional[str],
    db: Session,
) -> IPLink:
    """
    Link experiment/compound/dataset to disclosure.

    Args:
        entity_type: Type of entity (compound, experiment, dataset)
        entity_id: Entity UUID
        disclosure_id: Disclosure UUID
        link_type: Link type (embodiment, prior_art, enablement)
        notes: Optional notes
        db: Database session

    Returns:
        Created IPLink
    """
    link = IPLink(
        entity_type=entity_type,
        entity_id=entity_id,
        disclosure_id=disclosure_id,
        link_type=link_type,
        notes=notes,
    )
    
    db.add(link)
    db.commit()
    db.refresh(link)
    
    logger.info("[IP] Linked %s %s to disclosure %s", entity_type, entity_id, disclosure_id)
    
    return link


def get_disclosure_portfolio(
    company_id: UUID,
    status_filter: Optional[str],
    db: Session,
) -> List[InventionDisclosure]:
    """
    Get all disclosures for company with optional status filter.

    Args:
        company_id: Company UUID
        status_filter: Optional status to filter by
        db: Database session

    Returns:
        List of InventionDisclosure records
    """
    query = db.query(InventionDisclosure).filter(InventionDisclosure.company_id == company_id)
    
    if status_filter:
        query = query.filter(InventionDisclosure.status == status_filter)
    
    return query.order_by(InventionDisclosure.created_at.desc()).all()


def get_patent_timeline(
    disclosure_id: UUID,
    db: Session,
) -> List[PatentApplication]:
    """
    Get patent applications timeline for disclosure.

    Args:
        disclosure_id: Disclosure UUID
        db: Database session

    Returns:
        List of PatentApplication records ordered by filing date
    """
    return (
        db.query(PatentApplication)
        .filter(PatentApplication.disclosure_id == disclosure_id)
        .order_by(PatentApplication.filing_date.asc())
        .all()
    )


def get_entity_ip_links(
    entity_type: str,
    entity_id: UUID,
    db: Session,
) -> List[IPLink]:
    """
    Get IP links for a specific entity.

    Args:
        entity_type: Type of entity
        entity_id: Entity UUID
        db: Database session

    Returns:
        List of IPLink records
    """
    return (
        db.query(IPLink)
        .filter(IPLink.entity_type == entity_type, IPLink.entity_id == entity_id)
        .all()
    )

