"""
Compatibility wrapper for the refactored mwTab extraction pipeline.

This module maintains backward compatibility by re-exporting all public APIs
from the new modular structure.

Scripts and other code importing from mwtab_extraction.py will continue to work
without changes.
"""

from __future__ import annotations

from amprenta_rag.ingestion.mwtab import (
    extract_metadata_from_mwtab,
    extract_mwtab_from_page_content,
    extract_study_id_from_page_properties,
    fetch_mwtab_from_api,
)

__all__ = [
    "extract_mwtab_from_page_content",
    "extract_metadata_from_mwtab",
    "extract_study_id_from_page_properties",
    "fetch_mwtab_from_api",
]
