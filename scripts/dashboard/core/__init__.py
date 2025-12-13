"""Dashboard core package."""

from .config import (
    AUTH_DISABLED,
    DISCOVERY_PAGES,
    ANALYSIS_PAGES,
    ELN_PAGES,
    ADMIN_PAGES,
    ALL_PAGES,
    PAGE_REGISTRY,
)
from .auth import check_authentication, handle_session_timeout, get_mock_user
from .sidebar import render_sidebar
from .favorites import get_user_favorites, toggle_favorite, update_recent_pages
from .routing import route_to_page

__all__ = [
    "AUTH_DISABLED",
    "DISCOVERY_PAGES",
    "ANALYSIS_PAGES",
    "ELN_PAGES",
    "ADMIN_PAGES",
    "ALL_PAGES",
    "PAGE_REGISTRY",
    "check_authentication",
    "handle_session_timeout",
    "get_mock_user",
    "render_sidebar",
    "get_user_favorites",
    "toggle_favorite",
    "update_recent_pages",
    "route_to_page",
]

