"""
mwTab API fetching utilities.

This module provides functions for fetching mwTab data from the Metabolomics Workbench API.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any, Dict, Optional

from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def fetch_mwtab_from_api(study_id: str) -> Optional[Dict[str, Any]]:
    """
    Fetch mwTab data from Metabolomics Workbench API and parse as JSON.

    Optimized: Direct JSON parse first, fallback to extraction only if needed.

    Args:
        study_id: MW study ID (e.g., "ST001111")

    Returns:
        Parsed mwTab JSON dict or None if fetch/parse failed
    """
    try:
        # Import fetch function from harvest script
        scripts_path = Path(__file__).parent.parent.parent.parent / "scripts"
        if str(scripts_path) not in sys.path:
            sys.path.insert(0, str(scripts_path))
        from harvest_mw_studies import fetch_mw_mwtab

        mwtab_text = fetch_mw_mwtab(study_id)
        if not mwtab_text:
            logger.warning("[INGEST][MWTAB] MW API returned empty response")
            return None

        logger.debug("[INGEST][MWTAB] Fetched mwTab (%d chars)", len(mwtab_text))

        # OPTIMIZATION: Try direct JSON parse first (fastest path)
        try:
            mwtab_data = json.loads(mwtab_text)
            if isinstance(mwtab_data, dict):
                logger.info("[INGEST][MWTAB] Direct JSON parse successful")
                return mwtab_data
        except json.JSONDecodeError:
            pass  # Fall through to extraction methods

        # Find JSON object in response
        first_brace = mwtab_text.find("{")
        if first_brace < 0:
            logger.warning("[INGEST][MWTAB] No JSON object found in response")
            return None

        # OPTIMIZATION: Try parsing from first brace directly
        try:
            mwtab_data = json.loads(mwtab_text[first_brace:])
            if isinstance(mwtab_data, dict):
                logger.info("[INGEST][MWTAB] JSON parse from first brace successful")
                return mwtab_data
        except json.JSONDecodeError:
            pass

        # Fallback: Brace-counting for malformed responses
        brace_count = 0
        json_end = first_brace
        for i, char in enumerate(mwtab_text[first_brace:], start=first_brace):
            if char == "{":
                brace_count += 1
            elif char == "}":
                brace_count -= 1
                if brace_count == 0:
                    json_end = i + 1
                    break

        if json_end > first_brace:
            try:
                mwtab_data = json.loads(mwtab_text[first_brace:json_end])
                if isinstance(mwtab_data, dict):
                    logger.info("[INGEST][MWTAB] Brace-counting parse successful")
                    return mwtab_data
            except json.JSONDecodeError as e:
                logger.warning("[INGEST][MWTAB] Brace-counting parse failed: %s", e)

        logger.warning("[INGEST][MWTAB] All parsing methods failed")
        return None

    except Exception as fetch_error:
        logger.warning("[INGEST][MWTAB] Fetch failed: %r", fetch_error)
        # Try to extract HTTP status if available
        if (
            hasattr(fetch_error, "response")
            and fetch_error.response is not None
        ):
            logger.warning(
                "[INGEST][MWTAB] HTTP Status: %s, Response: %s",
                fetch_error.response.status_code,
                (
                    fetch_error.response.text[:500]
                    if hasattr(fetch_error.response, "text")
                    else "N/A"
                ),
            )
        return None

