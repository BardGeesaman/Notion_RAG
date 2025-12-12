"""Feature visibility permissions by role."""
from __future__ import annotations

from typing import Dict, List

from amprenta_rag.database.models import FeaturePermission
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def get_visible_pages(role: str, db) -> List[str]:
    """
    Get list of visible pages for a role.
    
    Args:
        role: User role (admin, researcher, viewer)
        db: Database session
        
    Returns:
        List of page names that should be visible to this role.
        If no permissions are set, returns all pages (default: all visible).
    """
    # Check if any permissions exist for this role
    has_permissions = db.query(FeaturePermission).filter(
        FeaturePermission.role == role
    ).first() is not None
    
    if not has_permissions:
        # No permissions set - default to all pages visible
        # Return all pages from ALL_PAGES
        from scripts.run_dashboard import ALL_PAGES
        return ALL_PAGES
    
    # Get visible pages for this role
    permissions = db.query(FeaturePermission).filter(
        FeaturePermission.role == role,
        FeaturePermission.is_visible == True
    ).all()
    
    return [perm.page_name for perm in permissions]


def set_page_visibility(role: str, page_name: str, visible: bool, db) -> bool:
    """
    Set visibility for a page for a specific role.
    
    Args:
        role: User role (admin, researcher, viewer)
        page_name: Name of the page
        visible: Whether the page should be visible
        db: Database session
        
    Returns:
        True if successful, False otherwise
    """
    try:
        # Check if permission already exists
        existing = db.query(FeaturePermission).filter(
            FeaturePermission.role == role,
            FeaturePermission.page_name == page_name
        ).first()
        
        if existing:
            existing.is_visible = visible
        else:
            permission = FeaturePermission(
                role=role,
                page_name=page_name,
                is_visible=visible,
            )
            db.add(permission)
        
        db.commit()
        return True
    except Exception as e:
        logger.error("[PERMISSIONS] Error setting page visibility: %s", e)
        db.rollback()
        return False


def get_all_permissions(db) -> Dict[str, Dict[str, bool]]:
    """
    Get all permissions organized by role.
    
    Args:
        db: Database session
        
    Returns:
        Dictionary mapping role -> {page_name: is_visible, ...}
    """
    permissions = db.query(FeaturePermission).all()
    
    result: Dict[str, Dict[str, bool]] = {}
    
    for perm in permissions:
        if perm.role not in result:
            result[perm.role] = {}
        result[perm.role][perm.page_name] = perm.is_visible
    
    return result


def initialize_default_permissions(db) -> None:
    """
    Initialize default permissions if none exist.
    
    Defaults:
    - Admin: all pages visible
    - Researcher: all pages except admin-only pages
    - Viewer: read-only pages only
    
    Args:
        db: Database session
    """
    from scripts.run_dashboard import ALL_PAGES, ADMIN_PAGES
    
    # Check if any permissions exist
    existing_count = db.query(FeaturePermission).count()
    if existing_count > 0:
        logger.info("[PERMISSIONS] Permissions already exist, skipping initialization")
        return
    
    logger.info("[PERMISSIONS] Initializing default permissions")
    
    # Read-only pages (viewer can access)
    read_only_pages = [
        "Overview",
        "Getting Started",
        "Datasets",
        "Programs",
        "Experiments",
        "Features",
        "Signatures",
        "Literature",
        "Emails",
        "RAG Chunks",
        "RAG Query",
        "Cross-Omics",
        "Timeline",
        "Compare",
        "Data Quality",
        "Literature Analysis",
        "Candidate Selection",
        "Data Lineage",
    ]
    
    # Admin-only pages
    admin_only_pages = ADMIN_PAGES
    
    # All pages
    all_pages = ALL_PAGES
    
    try:
        # Admin: all pages visible
        for page_name in all_pages:
            perm = FeaturePermission(
                role="admin",
                page_name=page_name,
                is_visible=True,
            )
            db.add(perm)
        
        # Researcher: all pages except admin-only
        researcher_pages = [p for p in all_pages if p not in admin_only_pages]
        for page_name in researcher_pages:
            perm = FeaturePermission(
                role="researcher",
                page_name=page_name,
                is_visible=True,
            )
            db.add(perm)
        
        # Researcher: admin pages hidden
        for page_name in admin_only_pages:
            perm = FeaturePermission(
                role="researcher",
                page_name=page_name,
                is_visible=False,
            )
            db.add(perm)
        
        # Viewer: read-only pages visible
        for page_name in read_only_pages:
            if page_name in all_pages:
                perm = FeaturePermission(
                    role="viewer",
                    page_name=page_name,
                    is_visible=True,
                )
                db.add(perm)
        
        # Viewer: other pages hidden
        viewer_hidden_pages = [p for p in all_pages if p not in read_only_pages]
        for page_name in viewer_hidden_pages:
            perm = FeaturePermission(
                role="viewer",
                page_name=page_name,
                is_visible=False,
            )
            db.add(perm)
        
        db.commit()
        logger.info("[PERMISSIONS] Default permissions initialized successfully")
    except Exception as e:
        logger.error("[PERMISSIONS] Error initializing default permissions: %s", e)
        db.rollback()
        raise

