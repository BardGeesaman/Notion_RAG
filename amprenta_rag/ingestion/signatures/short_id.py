"""
Short ID generation for signatures.

This module provides utilities for creating consistent, URL-safe short identifiers
for signatures that can be used for idempotent Notion page creation and lookup.

Key Functions:
    - generate_signature_short_id: Creates a normalized short ID from a signature name
      and optional version, handling special characters and length limits.

Example:
    >>> from amprenta_rag.ingestion.signatures.short_id import generate_signature_short_id
    >>> generate_signature_short_id("ALS-CSF-Core-6Ceramides")
    'ALS-CSF-Core-6Ceramides'
    >>> generate_signature_short_id("Test Signature", version="1.0")
    'Test-Signature-v1.0'
"""

from __future__ import annotations

import re
from typing import Optional


def generate_signature_short_id(signature_name: str, version: Optional[str] = None) -> str:
    """
    Generate a deterministic Short ID for a signature.

    Args:
        signature_name: Signature name
        version: Optional version string

    Returns:
        Deterministic Short ID (e.g., "ALS-CSF-Core-6Cer-v1")
    """
    # Normalize name for ID generation
    normalized = re.sub(r"[^a-zA-Z0-9]+", "-", signature_name).strip("-")
    normalized = normalized[:50]  # Limit length

    if version:
        return f"{normalized}-v{version}"

    # Generate version-less ID
    return normalized

