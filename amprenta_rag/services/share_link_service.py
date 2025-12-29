"""
Share link service for Voila dashboard sharing.

Provides secure link generation, validation, and lifecycle management.
"""

from __future__ import annotations

import secrets
from datetime import datetime, timedelta
from typing import List, Optional
from uuid import UUID

from sqlalchemy.orm import Session

from amprenta_rag.models.share_links import ShareLink
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def generate_share_link(
    dashboard_path: str,
    context: Optional[dict],
    user_id: UUID,
    company_id: UUID,
    expires_in_hours: int,
    max_views: Optional[int],
    permissions: str,
    db: Session,
) -> ShareLink:
    """
    Create share link with secure token.

    Args:
        dashboard_path: Path to Voila dashboard
        context: Dashboard context parameters
        user_id: Creating user UUID
        company_id: Company UUID
        expires_in_hours: Hours until expiration
        max_views: Maximum view count (None = unlimited)
        permissions: Permission level (view, edit)
        db: Database session

    Returns:
        Created ShareLink
    """
    # Generate secure token (48 bytes -> 64 URL-safe chars)
    token = secrets.token_urlsafe(48)
    
    expires_at = datetime.utcnow() + timedelta(hours=expires_in_hours)
    
    # Serialize context
    import json
    context_json = json.dumps(context) if context else None
    
    link = ShareLink(
        token=token,
        dashboard_path=dashboard_path,
        context_json=context_json,
        created_by_id=user_id,
        company_id=company_id,
        expires_at=expires_at,
        max_views=max_views,
        permissions=permissions,
    )
    
    db.add(link)
    db.commit()
    db.refresh(link)
    
    logger.info("[SHARE] Created share link: %s (expires: %s)", token[:16], expires_at)
    
    return link


def validate_share_link(token: str, db: Session) -> Optional[dict]:
    """
    Validate token and return dashboard info.

    Args:
        token: Share link token
        db: Database session

    Returns:
        Dict with dashboard_path, context, permissions if valid, None otherwise
    """
    link = db.query(ShareLink).filter(ShareLink.token == token).first()
    
    if not link:
        logger.warning("[SHARE] Token not found: %s", token[:16])
        return None
    
    # Check if active
    if not link.is_active:
        logger.warning("[SHARE] Link inactive: %s", token[:16])
        return None
    
    # Check expiration
    if datetime.utcnow() > link.expires_at:
        logger.warning("[SHARE] Link expired: %s", token[:16])
        link.is_active = False
        db.commit()
        return None
    
    # Check max views
    if link.max_views and link.view_count >= link.max_views:
        logger.warning("[SHARE] Max views exceeded: %s", token[:16])
        link.is_active = False
        db.commit()
        return None
    
    # Valid - increment view count
    link.view_count = (link.view_count or 0) + 1
    link.last_accessed_at = datetime.utcnow()
    db.commit()
    
    # Deserialize context
    import json
    context = json.loads(link.context_json) if link.context_json else {}
    
    logger.info("[SHARE] Validated link: %s (views: %d)", token[:16], link.view_count)
    
    return {
        "dashboard_path": link.dashboard_path,
        "context": context,
        "permissions": link.permissions,
    }


def revoke_share_link(link_id: UUID, user_id: UUID, db: Session) -> bool:
    """
    Revoke share link.

    Args:
        link_id: ShareLink UUID
        user_id: User UUID (must be owner)
        db: Database session

    Returns:
        True if revoked, False if not found or not owner
    """
    link = db.query(ShareLink).filter(ShareLink.id == link_id).first()
    
    if not link:
        return False
    
    # Check ownership
    if link.created_by_id != user_id:
        logger.warning("[SHARE] Revoke denied - not owner: %s", link_id)
        return False
    
    link.is_active = False
    db.commit()
    
    logger.info("[SHARE] Revoked link: %s", link.id)
    
    return True


def get_user_share_links(user_id: UUID, db: Session) -> List[ShareLink]:
    """
    Get all share links created by user.

    Args:
        user_id: User UUID
        db: Database session

    Returns:
        List of ShareLink records
    """
    return (
        db.query(ShareLink)
        .filter(ShareLink.created_by_id == user_id)
        .order_by(ShareLink.created_at.desc())
        .all()
    )


def cleanup_expired_links(db: Session) -> int:
    """
    Deactivate expired share links.

    Args:
        db: Database session

    Returns:
        Count of deactivated links
    """
    count = (
        db.query(ShareLink)
        .filter(
            ShareLink.expires_at < datetime.utcnow(),
            ShareLink.is_active.is_(True),
        )
        .update({"is_active": False}, synchronize_session=False)
    )
    
    db.commit()
    
    logger.info("[SHARE] Cleaned up %d expired links", count)
    
    return count


def get_share_link_stats(link_id: UUID, db: Session) -> Optional[dict]:
    """
    Get share link usage statistics.

    Args:
        link_id: ShareLink UUID
        db: Database session

    Returns:
        Dict with view_count, last_accessed_at, time_remaining
    """
    link = db.query(ShareLink).filter(ShareLink.id == link_id).first()
    
    if not link:
        return None
    
    time_remaining = (link.expires_at - datetime.utcnow()).total_seconds() / 3600  # Hours
    
    return {
        "view_count": link.view_count or 0,
        "max_views": link.max_views,
        "last_accessed_at": link.last_accessed_at.isoformat() if link.last_accessed_at else None,
        "time_remaining_hours": max(0, time_remaining),
        "is_active": link.is_active,
    }

