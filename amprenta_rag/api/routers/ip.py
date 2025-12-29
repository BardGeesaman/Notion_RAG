"""
IP and patent tracking API endpoints.

Provides endpoints for managing invention disclosures, patents, and entity links.
"""

from __future__ import annotations

from datetime import datetime
from typing import List, Optional
from uuid import UUID

from fastapi import APIRouter, Depends, Query
from pydantic import BaseModel
from sqlalchemy.orm import Session

from amprenta_rag.api.dependencies import get_database_session
from amprenta_rag.models.ip import InventionDisclosure, PatentApplication
from amprenta_rag.services.ip_service import (
    create_disclosure,
    file_patent,
    get_disclosure_portfolio,
    get_entity_ip_links,
    get_patent_timeline,
    link_entity_to_disclosure,
    update_disclosure_status,
)

router = APIRouter()


class InventorRequest(BaseModel):
    """Inventor info for disclosure."""

    user_id: UUID
    contribution_percentage: int | None = None
    is_primary: bool = False


class DisclosureCreate(BaseModel):
    """Request for creating disclosure."""

    title: str
    description: str
    technical_field: str | None = None
    inventors: List[InventorRequest]
    user_id: UUID
    company_id: UUID


class DisclosureResponse(BaseModel):
    """Response for disclosure."""

    id: UUID
    title: str
    description: str
    status: str
    technical_field: str | None
    priority_date: str | None
    created_at: str


class StatusUpdate(BaseModel):
    """Request for status update."""

    new_status: str
    user_id: UUID


class PatentFileRequest(BaseModel):
    """Request for filing patent."""

    application_number: str
    jurisdiction: str
    filing_date: datetime


class PatentResponse(BaseModel):
    """Response for patent application."""

    id: UUID
    application_number: str
    jurisdiction: str
    status: str
    filing_date: str
    grant_date: str | None


class LinkCreate(BaseModel):
    """Request for creating IP link."""

    entity_type: str
    entity_id: UUID
    disclosure_id: UUID
    link_type: str
    notes: str | None = None


class LinkResponse(BaseModel):
    """Response for IP link."""

    id: UUID
    entity_type: str
    entity_id: UUID
    disclosure_id: UUID
    link_type: str
    notes: str | None


@router.post("/disclosures", response_model=DisclosureResponse, status_code=201)
async def create_disclosure_endpoint(
    request: DisclosureCreate,
    db: Session = Depends(get_database_session),
) -> DisclosureResponse:
    """Create invention disclosure."""
    inventors_data = [inv.dict() for inv in request.inventors]
    
    disclosure = create_disclosure(
        request.title,
        request.description,
        request.technical_field,
        inventors_data,
        request.user_id,
        request.company_id,
        db,
    )
    
    return DisclosureResponse(
        id=disclosure.id,
        title=disclosure.title,
        description=disclosure.description,
        status=disclosure.status,
        technical_field=disclosure.technical_field,
        priority_date=disclosure.priority_date.isoformat() if disclosure.priority_date else None,
        created_at=disclosure.created_at.isoformat() if disclosure.created_at else "",
    )


@router.get("/disclosures", response_model=List[DisclosureResponse])
async def list_disclosures(
    company_id: UUID,
    status: str | None = Query(None),
    db: Session = Depends(get_database_session),
) -> List[DisclosureResponse]:
    """List disclosures for company."""
    disclosures = get_disclosure_portfolio(company_id, status, db)
    
    return [
        DisclosureResponse(
            id=d.id,
            title=d.title,
            description=d.description,
            status=d.status,
            technical_field=d.technical_field,
            priority_date=d.priority_date.isoformat() if d.priority_date else None,
            created_at=d.created_at.isoformat() if d.created_at else "",
        )
        for d in disclosures
    ]


@router.get("/disclosures/{disclosure_id}", response_model=DisclosureResponse)
async def get_disclosure_details(
    disclosure_id: UUID,
    db: Session = Depends(get_database_session),
) -> DisclosureResponse:
    """Get disclosure details."""
    disclosure = db.query(InventionDisclosure).filter(InventionDisclosure.id == disclosure_id).first()
    
    if not disclosure:
        from fastapi import HTTPException
        raise HTTPException(status_code=404, detail="Disclosure not found")
    
    return DisclosureResponse(
        id=disclosure.id,
        title=disclosure.title,
        description=disclosure.description,
        status=disclosure.status,
        technical_field=disclosure.technical_field,
        priority_date=disclosure.priority_date.isoformat() if disclosure.priority_date else None,
        created_at=disclosure.created_at.isoformat() if disclosure.created_at else "",
    )


@router.put("/disclosures/{disclosure_id}/status", response_model=DisclosureResponse)
async def update_status(
    disclosure_id: UUID,
    request: StatusUpdate,
    db: Session = Depends(get_database_session),
) -> DisclosureResponse:
    """Update disclosure status."""
    disclosure = update_disclosure_status(disclosure_id, request.new_status, request.user_id, db)
    
    return DisclosureResponse(
        id=disclosure.id,
        title=disclosure.title,
        description=disclosure.description,
        status=disclosure.status,
        technical_field=disclosure.technical_field,
        priority_date=disclosure.priority_date.isoformat() if disclosure.priority_date else None,
        created_at=disclosure.created_at.isoformat() if disclosure.created_at else "",
    )


@router.post("/disclosures/{disclosure_id}/file-patent", response_model=PatentResponse, status_code=201)
async def file_patent_endpoint(
    disclosure_id: UUID,
    request: PatentFileRequest,
    db: Session = Depends(get_database_session),
) -> PatentResponse:
    """File patent from disclosure."""
    patent = file_patent(
        disclosure_id,
        request.application_number,
        request.jurisdiction,
        request.filing_date,
        db,
    )
    
    return PatentResponse(
        id=patent.id,
        application_number=patent.application_number,
        jurisdiction=patent.jurisdiction,
        status=patent.status,
        filing_date=patent.filing_date.isoformat() if patent.filing_date else "",
        grant_date=patent.grant_date.isoformat() if patent.grant_date else None,
    )


@router.get("/patents", response_model=List[PatentResponse])
async def list_patents(
    disclosure_id: UUID | None = Query(None),
    status: str | None = Query(None),
    db: Session = Depends(get_database_session),
) -> List[PatentResponse]:
    """List patents."""
    query = db.query(PatentApplication)
    
    if disclosure_id:
        query = query.filter(PatentApplication.disclosure_id == disclosure_id)
    if status:
        query = query.filter(PatentApplication.status == status)
    
    patents = query.all()
    
    return [
        PatentResponse(
            id=p.id,
            application_number=p.application_number,
            jurisdiction=p.jurisdiction,
            status=p.status,
            filing_date=p.filing_date.isoformat() if p.filing_date else "",
            grant_date=p.grant_date.isoformat() if p.grant_date else None,
        )
        for p in patents
    ]


@router.get("/patents/{patent_id}/claims")
async def get_patent_claims_endpoint(
    patent_id: UUID,
    db: Session = Depends(get_database_session),
):
    """Get claims for patent."""
    from amprenta_rag.models.ip import PatentClaim
    
    claims = db.query(PatentClaim).filter(PatentClaim.patent_id == patent_id).order_by(PatentClaim.claim_number).all()
    
    return [
        {
            "id": str(c.id),
            "claim_number": c.claim_number,
            "claim_text": c.claim_text,
            "claim_type": c.claim_type,
        }
        for c in claims
    ]


@router.post("/links", response_model=LinkResponse, status_code=201)
async def create_link(
    request: LinkCreate,
    db: Session = Depends(get_database_session),
) -> LinkResponse:
    """Create IP link."""
    link = link_entity_to_disclosure(
        request.entity_type,
        request.entity_id,
        request.disclosure_id,
        request.link_type,
        request.notes,
        db,
    )
    
    return LinkResponse(
        id=link.id,
        entity_type=link.entity_type,
        entity_id=link.entity_id,
        disclosure_id=link.disclosure_id,
        link_type=link.link_type,
        notes=link.notes,
    )


@router.get("/links/{entity_type}/{entity_id}", response_model=List[LinkResponse])
async def get_links(
    entity_type: str,
    entity_id: UUID,
    db: Session = Depends(get_database_session),
) -> List[LinkResponse]:
    """Get IP links for entity."""
    links = get_entity_ip_links(entity_type, entity_id, db)
    
    return [
        LinkResponse(
            id=link.id,
            entity_type=link.entity_type,
            entity_id=link.entity_id,
            disclosure_id=link.disclosure_id,
            link_type=link.link_type,
            notes=link.notes,
        )
        for link in links
    ]

