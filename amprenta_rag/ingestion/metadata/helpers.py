"""
Helper functions for metadata extraction.

Provides utilities for fetching Notion pages and extracting property values.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional

from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def get_select_name(props: Dict[str, Any], name: str) -> Optional[str]:
    """Extract a single select value from Notion properties."""
    sel = props.get(name, {}).get("select")
    if sel and sel.get("name"):
        return sel["name"]
    return None


def get_multi_names(props: Dict[str, Any], name: str) -> List[str]:
    """Extract multi-select values from Notion properties."""
    ms = props.get(name, {}).get("multi_select", []) or []
    return [x["name"] for x in ms if x.get("name")]


def get_relation_ids(props: Dict[str, Any], name: str) -> List[str]:
    """Extract page IDs from a Notion relation property."""
    rel = props.get(name, {}).get("relation", []) or []
    return [r.get("id", "").replace("-", "") for r in rel if r.get("id")]

