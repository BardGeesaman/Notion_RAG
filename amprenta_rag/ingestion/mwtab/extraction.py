"""
mwTab extraction from page content and API.

This module provides functions for extracting mwTab JSON data from various sources.
"""

from __future__ import annotations

import json
from typing import Any, Dict, Optional

from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def extract_mwtab_from_page_content(page_content: str) -> Optional[Dict[str, Any]]:
    """
    Extract mwTab JSON from page content.

    Searches all code blocks and text sections for mwTab JSON data.
    More robust: checks code blocks, headings, and text content.

    Args:
        page_content: Full text content from extract_page_content()

    Returns:
        Parsed mwTab JSON dict or None if not found
    """
    if not page_content:
        logger.debug("[INGEST][MWTAB] Page content is empty")
        return None

    # Strategy 1: Look for code blocks in the text
    # Split by potential code block markers
    page_content.split("\n")

    # Look for JSON-like content (starts with {, contains common mwTab keys)
    mwtab_indicators = [
        "METABOLOMICS WORKBENCH",
        "MS_METABOLITE_DATA",
        "STUDY_ID",
        "SUBJECT_SAMPLE_FACTORS",
    ]

    # Try to find JSON content in the text
    text_content = page_content
    first_brace = text_content.find("{")

    if first_brace >= 0:
        logger.debug(
            "[INGEST][MWTAB] Found code block with '{', attempting JSON parse..."
        )

        # Try to extract complete JSON object
        # Start from first brace and try different end positions
        remaining_text = text_content[first_brace:]

        # Try parsing with different end positions (work backwards from end)
        for end_offset in range(len(remaining_text), 100, -100):
            try:
                json_str = remaining_text[:end_offset].strip()
                # Clean up: remove trailing incomplete JSON
                json_str = json_str.rstrip().rstrip(",")

                data = json.loads(json_str)
                if isinstance(data, dict):
                    # Check if it looks like mwTab data
                    if any(indicator in str(data) for indicator in mwtab_indicators):
                        logger.info(
                            "[INGEST][MWTAB] Successfully parsed mwTab JSON from page content"
                        )
                        return data
            except (json.JSONDecodeError, ValueError):
                continue

        # Try brace-counting method for more precise extraction
        brace_count = 0
        json_end = first_brace
        for i, char in enumerate(text_content[first_brace:], start=first_brace):
            if char == "{":
                brace_count += 1
            elif char == "}":
                brace_count -= 1
                if brace_count == 0:
                    json_end = i + 1
                    break

        if json_end > first_brace:
            try:
                json_str = text_content[first_brace:json_end].strip()
                data = json.loads(json_str)
                if isinstance(data, dict):
                    if any(indicator in str(data) for indicator in mwtab_indicators):
                        logger.info(
                            "[INGEST][MWTAB] Successfully parsed mwTab JSON from page content (brace-counting)"
                        )
                        return data
            except (json.JSONDecodeError, ValueError) as e:
                logger.debug(
                    "[INGEST][MWTAB] Failed JSON parse on block: %r", str(e)[:100]
                )

    logger.debug("[INGEST][MWTAB] No valid mwTab JSON found in page content")
    return None

