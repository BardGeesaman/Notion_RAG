"""
Share links API endpoints for Voila dashboard sharing.

Provides endpoints for creating, validating, and managing share links.
"""

from __future__ import annotations

from typing import List
from uuid import UUID

from fastapi import APIRouter, Depends, Request
from pydantic import BaseModel
from sqlalchemy.orm import Session

from amprenta_rag.api.dependencies import get_current_user, get_database_session
from amprenta_rag.models.auth import User
from amprenta_rag.services.share_link_service import (
    generate_share_link,
    get_share_link_stats,
    get_user_share_links,
    revoke_share_link,
    validate_share_link,
)
from amprenta_rag.api.rate_limit import limiter

router = APIRouter()


class ShareLinkCreate(BaseModel):
    """Request for creating share link."""

    dashboard_path: str
    context: dict | None = None
    expires_in_hours: int
    max_views: int | None = None
    permissions: str = "view"


class ShareLinkResponse(BaseModel):
    """Response for share link."""

    id: UUID
    token: str
    dashboard_path: str
    expires_at: str
    max_views: int | None
    view_count: int
    is_active: bool
    permissions: str
    created_at: str


class ShareLinkValidateResponse(BaseModel):
    """Response for share link validation."""

    valid: bool
    dashboard_path: str | None = None
    context: dict | None = None
    permissions: str | None = None
    voila_url: str | None = None


class ShareLinkStatsResponse(BaseModel):
    """Response for share link statistics."""

    view_count: int
    max_views: int | None
    last_accessed_at: str | None
    time_remaining_hours: float
    is_active: bool


@router.post("", response_model=ShareLinkResponse, status_code=201)
async def create_share_link(
    request: ShareLinkCreate,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_database_session),
) -> ShareLinkResponse:
    """Create share link for dashboard."""
    link = generate_share_link(
        request.dashboard_path,
        request.context,
        current_user.id,
        current_user.company_id,
        request.expires_in_hours,
        request.max_views,
        request.permissions,
        db,
    )
    
    return ShareLinkResponse(
        id=link.id,
        token=link.token,
        dashboard_path=link.dashboard_path,
        expires_at=link.expires_at.isoformat() if link.expires_at else "",
        max_views=link.max_views,
        view_count=link.view_count or 0,
        is_active=link.is_active,
        permissions=link.permissions,
        created_at=link.created_at.isoformat() if link.created_at else "",
    )


@router.get("", response_model=List[ShareLinkResponse])
async def list_share_links(
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_database_session),
) -> List[ShareLinkResponse]:
    """List user's share links."""
    links = get_user_share_links(current_user.id, db)
    
    return [
        ShareLinkResponse(
            id=link.id,
            token=link.token,
            dashboard_path=link.dashboard_path,
            expires_at=link.expires_at.isoformat() if link.expires_at else "",
            max_views=link.max_views,
            view_count=link.view_count or 0,
            is_active=link.is_active,
            permissions=link.permissions,
            created_at=link.created_at.isoformat() if link.created_at else "",
        )
        for link in links
    ]


@router.delete("/{link_id}")
async def revoke_link(
    link_id: UUID,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_database_session),
):
    """Revoke share link."""
    success = revoke_share_link(link_id, current_user.id, db)
    
    if not success:
        from fastapi import HTTPException
        raise HTTPException(status_code=404, detail="Link not found or not authorized")
    
    return {"revoked": True}


@router.get("/{token}/validate", response_model=ShareLinkValidateResponse)
@limiter.limit("30/minute")
async def validate_token(
    request: Request,
    token: str,
    db: Session = Depends(get_database_session),
) -> ShareLinkValidateResponse:
    """
    Validate share link token (PUBLIC endpoint - no auth required).

    Returns dashboard info if valid, or invalid response.
    """
    result = validate_share_link(token, db)
    
    if not result:
        return ShareLinkValidateResponse(valid=False)
    
    # Construct Voila URL
    voila_url = f"/voila/render{result['dashboard_path']}"
    if result.get("context"):
        # Add context as query params
        import urllib.parse
        params = urllib.parse.urlencode(result["context"])
        voila_url += f"?{params}"
    
    return ShareLinkValidateResponse(
        valid=True,
        dashboard_path=result["dashboard_path"],
        context=result.get("context"),
        permissions=result.get("permissions"),
        voila_url=voila_url,
    )


@router.get("/{link_id}/stats", response_model=ShareLinkStatsResponse)
async def get_stats(
    link_id: UUID,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_database_session),
) -> ShareLinkStatsResponse:
    """Get share link statistics."""
    stats = get_share_link_stats(link_id, db)
    
    if not stats:
        from fastapi import HTTPException
        raise HTTPException(status_code=404, detail="Link not found")
    
    return ShareLinkStatsResponse(**stats)

