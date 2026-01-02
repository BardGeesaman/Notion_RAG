"""HTML sanitization utilities for user-generated content.

Uses nh3 (Rust-backed, maintained by Mozilla) instead of deprecated bleach.
"""
import nh3
from typing import Set

# Allowed tags for rich text (comments, descriptions)
ALLOWED_TAGS: Set[str] = {'b', 'i', 'u', 'em', 'strong', 'a', 'br', 'p', 'span'}
ALLOWED_ATTRIBUTES: dict[str, Set[str]] = {
    'a': {'href', 'title'},
    'span': {'style'},
}


def sanitize_html(content: str) -> str:
    """Sanitize HTML content, removing dangerous tags/attributes.
    
    Args:
        content: User-provided HTML content
    
    Returns:
        Sanitized HTML safe for rendering
    """
    if not content:
        return content
    return nh3.clean(
        content,
        tags=ALLOWED_TAGS,
        attributes=ALLOWED_ATTRIBUTES,
    )


def escape_html(content: str) -> str:
    """Escape all HTML in content (for plain text display).
    
    Args:
        content: Content that may contain HTML
        
    Returns:
        Content with all HTML escaped
    """
    if not content:
        return content
    return nh3.clean(content, tags=set())


def escape_for_html_display(content: str) -> str:
    """Escape content for safe display in HTML context.
    
    Use this when rendering user content in st.markdown with unsafe_allow_html.
    
    Args:
        content: User-provided content
        
    Returns:
        Content safe for HTML display
    """
    return escape_html(content)
