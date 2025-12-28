"""
Metadata utility functions.

This module provides metadata sanitization for vector store compatibility.
"""

from __future__ import annotations

from typing import Any, Dict


def sanitize_metadata(meta: Dict[str, Any]) -> Dict[str, Any]:
    """
    Ensure all metadata values are compatible with vector stores:
      - string
      - number (int/float)
      - boolean
      - list of strings

    Drop keys with None or empty lists.
    Convert non-string, non-primitive list items to strings.
    """
    cleaned: Dict[str, Any] = {}

    for key, value in meta.items():
        if value is None:
            # Vector stores do not allow null metadata
            continue

        # Primitive types are fine
        if isinstance(value, (str, bool, int, float)):
            cleaned[key] = value
            continue

        # Lists need to be list of strings
        if isinstance(value, list):
            string_list = [str(v) for v in value if v is not None]
            if string_list:
                cleaned[key] = string_list
            continue

        # Anything else (e.g. dicts) – convert to string for now
        cleaned[key] = str(value)

    return cleaned

