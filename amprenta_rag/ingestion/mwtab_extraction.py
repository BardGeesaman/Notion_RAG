"""
mwTab extraction and parsing utilities.

This module handles extraction of mwTab JSON data from various sources:
- Notion page content
- Metabolomics Workbench API
- Metadata extraction from mwTab structures
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path
from typing import Any, Dict, Optional

from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)

__all__ = [
    "extract_mwtab_from_page_content",
    "extract_metadata_from_mwtab",
    "extract_study_id_from_page_properties",
    "fetch_mwtab_from_api",
]


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
    lines = page_content.split("\n")

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
            except (json.JSONDecodeError, ValueError) as e:
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


def extract_metadata_from_mwtab(mwtab_data: Dict[str, Any]) -> Dict[str, Any]:
    """
    Extract scientific metadata from mwTab JSON data.

    Args:
        mwtab_data: Parsed mwTab JSON dictionary

    Returns:
        Dictionary with extracted metadata:
        - model_systems: List[str]
        - disease_terms: List[str]
        - matrix_terms: List[str]
        - methods: str
        - summary: str
        - results: str
        - conclusions: str
        - data_origin: str | None
        - dataset_source_type: str | None
        - source_url: str | None
    """
    metadata = {
        "model_systems": [],
        "disease_terms": [],
        "matrix_terms": [],
        "methods": "",
        "summary": "",
        "results": "",
        "conclusions": "",
        "data_origin": None,
        "dataset_source_type": None,
        "source_url": None,
    }

    # Extract model systems from SUBJECT
    if "SUBJECT" in mwtab_data:
        subject = mwtab_data["SUBJECT"]
        if isinstance(subject, dict):
            species = subject.get("SUBJECT_SPECIES") or subject.get("species") or ""
            if species:
                metadata["model_systems"] = [species.strip()]

            # Extract matrix from SUBJECT_TYPE
            subject_type = (
                subject.get("SUBJECT_TYPE") or subject.get("subject_type") or ""
            )
            if subject_type:
                metadata["matrix_terms"] = [subject_type.strip()]

    # Extract summary from STUDY
    if "STUDY" in mwtab_data:
        study = mwtab_data["STUDY"]
        if isinstance(study, dict):
            study_summary = study.get("STUDY_SUMMARY") or study.get("summary") or ""
            if study_summary:
                metadata["summary"] = study_summary.strip()

            # Try to extract disease terms from STUDY_TITLE or STUDY_SUMMARY
            study_title = study.get("STUDY_TITLE") or study.get("title") or ""
            if study_title:
                # Look for common disease terms in title/summary
                combined_text = (study_title + " " + study_summary).lower()

                # Common disease patterns (expandable)
                disease_patterns = [
                    r"\bfxs\b",
                    r"\bfragile x\b",
                    r"\bfragile-x\b",
                    r"\bals\b",
                    r"\bamyotrophic\b",
                    r"\balzheimer\b",
                    r"\bad\b",
                    r"\bparkinson\b",
                    r"\bpd\b",
                    r"\bhuntington\b",
                ]

                disease_terms = []
                if "fxs" in combined_text or "fragile x" in combined_text:
                    disease_terms.append("Fragile X syndrome")
                if "als" in combined_text or "amyotrophic" in combined_text:
                    disease_terms.append("ALS")
                if "alzheimer" in combined_text or " ad " in combined_text:
                    disease_terms.append("Alzheimer's disease")

                if disease_terms:
                    metadata["disease_terms"] = disease_terms

    # Extract methods from TREATMENT, SAMPLEPREP, CHROMATOGRAPHY sections
    methods_parts = []

    if "TREATMENT" in mwtab_data:
        treatment = mwtab_data["TREATMENT"]
        if isinstance(treatment, dict):
            treatment_summary = (
                treatment.get("TREATMENT_SUMMARY") or treatment.get("summary") or ""
            )
            if treatment_summary:
                methods_parts.append(f"Treatment: {treatment_summary.strip()}")

    if "SAMPLEPREP" in mwtab_data:
        sampleprep = mwtab_data["SAMPLEPREP"]
        if isinstance(sampleprep, dict):
            sampleprep_summary = (
                sampleprep.get("SAMPLEPREP_SUMMARY") or sampleprep.get("summary") or ""
            )
            if sampleprep_summary:
                methods_parts.append(
                    f"Sample Preparation: {sampleprep_summary.strip()}"
                )

    if "CHROMATOGRAPHY" in mwtab_data:
        chrom = mwtab_data["CHROMATOGRAPHY"]
        if isinstance(chrom, dict):
            chrom_type = chrom.get("CHROMATOGRAPHY_TYPE") or chrom.get("type") or ""
            column = chrom.get("COLUMN_NAME") or chrom.get("column") or ""
            if chrom_type or column:
                methods_parts.append(f"Chromatography: {chrom_type} ({column})")

    if methods_parts:
        metadata["methods"] = "\n\n".join(methods_parts)

    # Results: Can be enhanced with LLM summarization later
    # For now, leave empty or add basic info about metabolite data
    has_processed_table = False
    if "MS_METABOLITE_DATA" in mwtab_data:
        ms_data = mwtab_data["MS_METABOLITE_DATA"]
        if isinstance(ms_data, dict):
            data_array = ms_data.get("Data", [])
            if isinstance(data_array, list) and len(data_array) > 0:
                metadata["results"] = (
                    f"Metabolite profiling data with {len(data_array)} metabolites detected."
                )
                has_processed_table = True

    # Extract conclusions (empty for now, can be enhanced with LLM later)
    metadata["conclusions"] = ""

    # Extract data_origin: Check if this is a Metabolomics Workbench dataset
    study_id = None
    if "METABOLOMICS WORKBENCH" in mwtab_data:
        mw = mwtab_data["METABOLOMICS WORKBENCH"]
        if isinstance(mw, dict):
            study_id_raw = mw.get("STUDY_ID") or ""
            if isinstance(study_id_raw, str) and study_id_raw.strip():
                study_id = study_id_raw.strip()
                metadata["data_origin"] = "External – Open Dataset"
    # Also check top-level STUDY_ID if METABOLOMICS WORKBENCH section missing
    if not metadata["data_origin"]:
        study_id_top = mwtab_data.get("STUDY_ID") or ""
        if isinstance(study_id_top, str) and study_id_top.strip().startswith("ST0"):
            study_id = study_id_top.strip()
            metadata["data_origin"] = "External – Open Dataset"

    # Extract dataset_source_type: Check if processed metabolite table exists
    if has_processed_table:
        metadata["dataset_source_type"] = "Processed table"

    # Extract source_url: Look for URL fields in mwTab
    url_candidates = []

    # Check STUDY section for URL/DOI fields
    if "STUDY" in mwtab_data and isinstance(mwtab_data["STUDY"], dict):
        study = mwtab_data["STUDY"]
        # Check for URL-like fields (case-insensitive matching)
        study_keys_lower = {k.lower(): k for k in study.keys()}
        for field_pattern in ["study_website", "study_doi", "study_url", "doi", "url"]:
            if field_pattern in study_keys_lower:
                key = study_keys_lower[field_pattern]
                value = study[key]
                if isinstance(value, str) and value.strip():
                    url_candidates.append(value.strip())

    # Check PROJECT section
    if "PROJECT" in mwtab_data and isinstance(mwtab_data["PROJECT"], dict):
        project = mwtab_data["PROJECT"]
        project_keys_lower = {k.lower(): k for k in project.keys()}
        for field_pattern in [
            "project_website",
            "project_doi",
            "project_url",
            "doi",
            "url",
        ]:
            if field_pattern in project_keys_lower:
                key = project_keys_lower[field_pattern]
                value = project[key]
                if isinstance(value, str) and value.strip():
                    url_candidates.append(value.strip())

    # If we have a STUDY_ID from Metabolomics Workbench, construct MW URL
    if study_id and isinstance(study_id, str) and study_id.strip().startswith("ST0"):
        url_candidates.append(
            f"https://www.metabolomicsworkbench.org/study/index.php?study_id={study_id.strip()}"
        )

    # Find first valid URL
    for candidate in url_candidates:
        candidate = candidate.strip()
        if not candidate:
            continue
        # Check if it's a valid URL
        if candidate.startswith("http://") or candidate.startswith("https://"):
            metadata["source_url"] = candidate
            break
        # Handle DOI format
        elif candidate.startswith("10."):
            metadata["source_url"] = f"https://doi.org/{candidate}"
            break
        elif candidate.lower().startswith("doi:"):
            metadata["source_url"] = f"https://doi.org/{candidate[4:].strip()}"
            break

    return metadata


def extract_study_id_from_page_properties(
    page_properties: Dict[str, Any], page_content: str = ""
) -> Optional[str]:
    """
    Extract STUDY_ID from Notion page properties or content.

    Searches multiple properties and page content for STUDY_ID pattern (ST######).

    Args:
        page_properties: Notion page properties dictionary
        page_content: Optional page content text to search

    Returns:
        STUDY_ID string (e.g., "ST001111") or None if not found
    """
    study_id = None

    # 1. Check Summary property
    summary_prop = page_properties.get("Summary", {})
    summary_rich = summary_prop.get("rich_text", []) or []
    summary_text = "".join(rt.get("plain_text", "") for rt in summary_rich)
    study_id_match = re.search(r"ST\d{6}", summary_text)
    if study_id_match:
        study_id = study_id_match.group(0)
        return study_id

    # 2. Check Title/Name property
    name_prop = (
        page_properties.get("Name", {})
        or page_properties.get("Experiment Name", {})
        or {}
    )
    if name_prop.get("title"):
        title_text = "".join(
            t.get("plain_text", "") for t in name_prop["title"]
        )
        study_id_match = re.search(r"ST\d{6}", title_text)
        if study_id_match:
            study_id = study_id_match.group(0)
            return study_id

    # 3. Check Source URL / DOI property
    url_prop = page_properties.get("Source URL / DOI", {}) or {}
    url_text = url_prop.get("url", "") or ""
    study_id_match = re.search(r"ST\d{6}", url_text)
    if study_id_match:
        study_id = study_id_match.group(0)
        return study_id

    # 4. Check page content itself
    if page_content:
        study_id_match = re.search(r"ST\d{6}", page_content)
        if study_id_match:
            study_id = study_id_match.group(0)
            return study_id

    return None


def fetch_mwtab_from_api(study_id: str) -> Optional[Dict[str, Any]]:
    """
    Fetch mwTab data from Metabolomics Workbench API and parse as JSON.

    Args:
        study_id: MW study ID (e.g., "ST001111")

    Returns:
        Parsed mwTab JSON dict or None if fetch/parse failed
    """
    try:
        # Import fetch function from harvest script
        scripts_path = Path(__file__).parent.parent.parent / "scripts"
        if str(scripts_path) not in sys.path:
            sys.path.insert(0, str(scripts_path))
        from harvest_mw_studies import fetch_mw_mwtab

        mwtab_text = fetch_mw_mwtab(study_id)
        if not mwtab_text:
            logger.warning(
                "[INGEST][MWTAB] MW API fallback fetch returned empty response"
            )
            return None

        logger.debug(
            "[INGEST][MWTAB] Fetched mwTab text from MW API (%d chars)",
            len(mwtab_text),
        )

        # Parse the fetched mwTab JSON
        first_brace = mwtab_text.find("{")
        if first_brace < 0:
            logger.warning(
                "[INGEST][MWTAB] MW API response does not contain JSON (no '{' found)"
            )
            return None

        # Try parsing from first brace
        remaining = mwtab_text[first_brace:]
        mwtab_data = None

        for end_offset in range(len(remaining), 100, -100):
            try:
                json_str = (
                    remaining[:end_offset]
                    .strip()
                    .rstrip(",")
                    .rstrip()
                )
                mwtab_data = json.loads(json_str)
                if isinstance(mwtab_data, dict):
                    logger.info(
                        "[INGEST][MWTAB] MW API fallback fetch successful. Parsed mwTab JSON."
                    )
                    return mwtab_data
            except (json.JSONDecodeError, ValueError):
                continue

        # If chunked parsing failed, try brace-counting
        if not mwtab_data:
            brace_count = 0
            json_end = first_brace
            for i, char in enumerate(
                mwtab_text[first_brace:], start=first_brace
            ):
                if char == "{":
                    brace_count += 1
                elif char == "}":
                    brace_count -= 1
                    if brace_count == 0:
                        json_end = i + 1
                        break

            if json_end > first_brace:
                try:
                    json_str = mwtab_text[first_brace:json_end]
                    mwtab_data = json.loads(json_str)
                    if isinstance(mwtab_data, dict):
                        logger.info(
                            "[INGEST][MWTAB] MW API fallback fetch successful. Parsed mwTab JSON (brace-counting)."
                        )
                        return mwtab_data
                except (
                    json.JSONDecodeError,
                    ValueError,
                ) as parse_err:
                    logger.warning(
                        "[INGEST][MWTAB] Failed to parse mwTab JSON from MW API response: %r",
                        str(parse_err)[:200],
                    )

        return None

    except Exception as fetch_error:
        logger.warning(
            "[INGEST][MWTAB] MW API fallback fetch FAILED: %r",
            fetch_error,
        )
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

