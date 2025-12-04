"""
mwTab extraction and parsing utilities.

This package provides functions for extracting mwTab JSON data from various sources
and extracting metadata from mwTab structures.

Maintains backward compatibility by re-exporting all public functions.
"""

from __future__ import annotations

from amprenta_rag.ingestion.mwtab.api import fetch_mwtab_from_api
from amprenta_rag.ingestion.mwtab.extraction import extract_mwtab_from_page_content
from amprenta_rag.ingestion.mwtab.metadata import (
    extract_metadata_from_mwtab,
    extract_study_id_from_page_properties,
)

__all__ = [
    "extract_mwtab_from_page_content",
    "extract_metadata_from_mwtab",
    "extract_study_id_from_page_properties",
    "fetch_mwtab_from_api",
]

