"""
Short ID generation for signatures.

Generates deterministic Short IDs for signature pages in Notion.
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

