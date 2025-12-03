# amprenta_rag/ingestion/dataset_ingestion.py
"""
Dataset ingestion module.

Handles ingestion of experimental datasets from Notion into Pinecone.
Supports mwTab data extraction, species extraction, signature matching,
and automatic metadata population.
"""

from __future__ import annotations

import json
import re
import textwrap
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional

import requests

from amprenta_rag.clients.notion_client import notion_headers
from amprenta_rag.clients.pinecone_client import get_pinecone_index
from amprenta_rag.config import get_config
from amprenta_rag.ingestion.feature_extraction import (
    extract_features_from_mwtab, link_features_to_notion_items)
from amprenta_rag.ingestion.metadata_semantic import \
    get_dataset_semantic_metadata
from amprenta_rag.ingestion.notion_pages import extract_page_content
from amprenta_rag.ingestion.pinecone_utils import sanitize_metadata
from amprenta_rag.ingestion.signature_integration import \
    detect_and_ingest_signatures_from_content
from amprenta_rag.ingestion.signature_matching import (
    find_matching_signatures_for_dataset, map_raw_lipid_to_canonical_species,
    update_dataset_with_signature_matches)
from amprenta_rag.ingestion.zotero_ingest import _chunk_text, _embed_texts
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def _fetch_notion_page(page_id: str) -> Dict[str, Any]:
    """
    Fetch a Notion page by ID.

    Args:
        page_id: Notion page ID (with or without dashes)

    Returns:
        Notion page JSON response

    Raises:
        Exception: If page fetch fails
    """
    from amprenta_rag.ingestion.metadata_semantic import \
        _fetch_notion_page as _fetch

    try:
        return _fetch(page_id)
    except Exception as e:
        logger.error("[DATASET][NOTION] Error fetching page %s: %r", page_id, e)
        raise


def _update_notion_page_with_embeddings(
    page_id: str,
    embedding_ids: List[str],
    embedding_count: int,
) -> None:
    """
    Update the Notion Dataset page with embedding metadata after successful Pinecone upsert.

    Args:
        page_id: Notion page ID (with or without dashes)
        embedding_ids: List of all vector/embedding IDs
        embedding_count: Number of chunks/vectors created
    """
    cfg = get_config()

    # Format embedding IDs as newline-separated string
    embedding_ids_text = "\n".join(embedding_ids)

    # Get current timestamp for Last Embedded
    last_embedded_iso = datetime.now(timezone.utc).isoformat()

    # Build properties payload
    # Note: Database has "Embedding IDs" and "Last Embedded" but not "Embedding Count"
    props: Dict[str, Any] = {
        "Embedding IDs": {
            "rich_text": [
                {
                    "type": "text",
                    "text": {
                        "content": embedding_ids_text,
                    },
                }
            ]
        },
        "Last Embedded": {
            "date": {
                "start": last_embedded_iso,
            },
        },
    }

    payload = {"properties": props}

    # Use page_id with dashes for Notion API
    url = f"{cfg.notion.base_url}/pages/{page_id}"

    try:
        resp = requests.patch(
            url,
            headers=notion_headers(),
            json=payload,
            timeout=30,
        )
        resp.raise_for_status()
        logger.info(
            "[INGEST][DATASET] Updated Notion page %s with %d embedding IDs.",
            page_id,
            embedding_count,
        )
    except Exception as e:
        if hasattr(e, "response") and e.response is not None:
            logger.error(
                "[INGEST][DATASET] Error updating Notion page %s: %s - Response: %s",
                page_id,
                str(e),
                e.response.text,
            )
        else:
            logger.error(
                "[INGEST][DATASET] Error updating Notion page %s: %r",
                page_id,
                e,
            )
        # Don't raise - ingestion succeeded, Notion update is metadata only


def _extract_mwtab_from_page_content(page_content: str) -> Optional[Dict[str, Any]]:
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


def _extract_metadata_from_mwtab(mwtab_data: Dict[str, Any]) -> Dict[str, Any]:
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


def _update_experimental_data_asset_metadata(
    page_id: str,
    metadata: Dict[str, Any],
) -> None:
    """
    Update Experimental Data Assets page with scientific metadata.

    Sets only these properties (does NOT touch Embedding IDs or Last Embedded):
    - Model Systems (multi_select)
    - Disease (multi_select)
    - Matrix (multi_select)
    - Methods (rich_text)
    - Summary (rich_text)
    - Results (rich_text)
    - Conclusions (rich_text)
    - Data Origin (select)
    - Dataset Source Type (select)
    - Source URL / DOI (url)

    Args:
        page_id: Notion page ID (with dashes)
        metadata: Dictionary with model_systems, disease_terms, matrix_terms, methods, summary,
                 results, conclusions, data_origin, dataset_source_type, source_url
    """
    cfg = get_config()

    props: Dict[str, Any] = {}

    # Model Systems (multi_select)
    if metadata.get("model_systems"):
        props["Model Systems"] = {
            "multi_select": [{"name": ms} for ms in metadata["model_systems"] if ms]
        }

    # Disease (multi_select)
    if metadata.get("disease_terms"):
        props["Disease"] = {
            "multi_select": [{"name": dt} for dt in metadata["disease_terms"] if dt]
        }

    # Matrix (multi_select)
    if metadata.get("matrix_terms"):
        props["Matrix"] = {
            "multi_select": [{"name": mt} for mt in metadata["matrix_terms"] if mt]
        }

    # Methods (rich_text)
    if metadata.get("methods"):
        props["Methods"] = {
            "rich_text": [{"type": "text", "text": {"content": metadata["methods"]}}]
        }

    # Summary (rich_text)
    if metadata.get("summary"):
        props["Summary"] = {
            "rich_text": [{"type": "text", "text": {"content": metadata["summary"]}}]
        }

    # Results (rich_text)
    if metadata.get("results"):
        props["Results"] = {
            "rich_text": [{"type": "text", "text": {"content": metadata["results"]}}]
        }

    # Conclusions (rich_text)
    if metadata.get("conclusions"):
        props["Conclusions"] = {
            "rich_text": [
                {"type": "text", "text": {"content": metadata["conclusions"]}}
            ]
        }

    # Data Origin (select)
    if metadata.get("data_origin"):
        props["Data Origin"] = {"select": {"name": metadata["data_origin"]}}

    # Dataset Source Type (select)
    if metadata.get("dataset_source_type"):
        props["Dataset Source Type"] = {
            "select": {"name": metadata["dataset_source_type"]}
        }

    # Source URL / DOI (url)
    if metadata.get("source_url"):
        props["Source URL / DOI"] = {"url": metadata["source_url"]}

    if not props:
        logger.info("[INGEST][DATASET] No metadata to update for page %s", page_id)
        return

    payload = {"properties": props}
    url = f"{cfg.notion.base_url}/pages/{page_id}"

    try:
        resp = requests.patch(
            url,
            headers=notion_headers(),
            json=payload,
            timeout=30,
        )
        resp.raise_for_status()
        logger.info(
            "[INGEST][DATASET] Updated scientific metadata for page %s",
            page_id,
        )
    except Exception as e:
        if hasattr(e, "response") and e.response is not None:
            logger.warning(
                "[INGEST][DATASET] Error updating scientific metadata for page %s: %s - Response: %s",
                page_id,
                str(e),
                e.response.text,
            )
        else:
            logger.warning(
                "[INGEST][DATASET] Error updating scientific metadata for page %s: %r",
                page_id,
                e,
            )
        # Don't raise - metadata update is non-critical


def ingest_dataset(page_id: str, force: bool = False) -> None:
    """
    Ingest a dataset page from Notion into Pinecone.

    This function performs the complete dataset ingestion pipeline:
    1. Fetches the dataset page from Notion
    2. Extracts mwTab data (from page content or MW API)
    3. Extracts full text content
    4. Extracts semantic metadata (diseases, matrix, signatures, etc.)
    5. Chunks and embeds the text
    6. Upserts vectors to Pinecone with metadata
    7. Updates Notion page with embedding IDs
    8. Extracts scientific metadata and updates Notion
    9. Extracts metabolite features and links them
    10. Detects and ingests signatures from content
    11. Matches dataset against lipid signatures and updates Notion

    Args:
        page_id: Notion page ID (with or without dashes)
        force: If True, re-ingest even if already embedded

    Raises:
        Exception: If ingestion fails at any step
    """
    """
    Ingest a single Dataset page from Notion into Pinecone.

    Args:
        page_id: Notion page ID (with or without dashes) for the Dataset page
        force: If True, re-ingest even if already embedded
    """
    logger.info("[INGEST][DATASET] Ingesting dataset page %s", page_id)

    try:
        page = _fetch_notion_page(page_id)
    except Exception as e:
        logger.error(
            "[INGEST][DATASET] Error fetching Notion page %s: %r",
            page_id,
            e,
        )
        raise

    # Get the canonical page ID from the fetched page (with dashes)
    canonical_page_id = page.get("id", page_id)

    props = page.get("properties", {}) or {}

    # Get title from "Name" property
    title_prop = props.get("Name", {})
    title_parts = title_prop.get("title", []) or []
    dataset_title = (
        title_parts[0].get("plain_text", "").strip()
        if title_parts
        else "(untitled dataset)"
    )

    # Extract full content from blocks (using shared helper)
    try:
        # Normalize page ID (remove dashes for extract_page_content)
        page_id_clean = page_id.replace("-", "")
        full_text = extract_page_content(page_id_clean)
    except Exception as e:
        logger.error(
            "[INGEST][DATASET] Error extracting content for %s: %r",
            page_id,
            e,
        )
        raise

    if not full_text or len(full_text.strip()) < 50:
        logger.info(
            "[INGEST][DATASET] Dataset %s has very little text; skipping.", page_id
        )
        return

    # Semantic + signature metadata
    try:
        base_meta = get_dataset_semantic_metadata(page)
        logger.info(
            "[INGEST][DATASET] Lipid signatures for %s: %r",
            page_id,
            base_meta.get("lipid_signatures"),
        )
    except Exception as e:
        logger.error(
            "[INGEST][DATASET] Error reading semantic metadata for %s: %r",
            page_id,
            e,
        )
        raise

    # Chunk and embed
    chunks = _chunk_text(full_text)
    if not chunks:
        logger.info("[INGEST][DATASET] No chunks produced for %s; skipping.", page_id)
        return

    logger.info(
        "[INGEST][DATASET] Generated %d chunk(s) for dataset %s",
        len(chunks),
        page_id,
    )

    try:
        embeddings = _embed_texts(chunks)
    except Exception as e:
        logger.error(
            "[INGEST][DATASET] Error embedding chunks for %s: %r",
            page_id,
            e,
        )
        raise

    index = get_pinecone_index()
    cfg = get_config()

    vectors: List[Dict[str, Any]] = []
    embedding_ids: List[str] = []

    for order, (chunk, emb) in enumerate(zip(chunks, embeddings)):
        chunk_id = f"{page_id.replace('-', '')}_chunk_{order:03d}"
        embedding_ids.append(chunk_id)

        snippet = textwrap.shorten(chunk, width=300)

        meta: Dict[str, Any] = {
            **base_meta,
            "source": "Dataset",
            "source_type": "Dataset",
            "dataset_page_id": page_id.replace("-", ""),
            "title": dataset_title,
            "snippet": snippet,
        }

        vectors.append(
            {
                "id": chunk_id,
                "values": emb,
                "metadata": sanitize_metadata(meta),
            }
        )

    if not vectors:
        logger.info("[INGEST][DATASET] No vectors to upsert for %s; skipping.", page_id)
        return

    logger.info(
        "[INGEST][DATASET] Upserting %d vectors into Pinecone for dataset %s",
        len(vectors),
        page_id,
    )

    # Batch upserts to avoid Pinecone size limits (2MB per request)
    # Rough estimate: ~15-20KB per vector (1536 dims * 4 bytes + metadata)
    # Use batches of ~100 vectors to stay well under 2MB limit
    batch_size = 100
    try:
        for i in range(0, len(vectors), batch_size):
            batch = vectors[i : i + batch_size]
            batch_num = (i // batch_size) + 1
            total_batches = (len(vectors) + batch_size - 1) // batch_size

            logger.debug(
                "[INGEST][DATASET] Upserting batch %d/%d (%d vectors) for dataset %s",
                batch_num,
                total_batches,
                len(batch),
                page_id,
            )

            index.upsert(vectors=batch, namespace=cfg.pinecone.namespace)

            if batch_num % 10 == 0 or batch_num == total_batches:
                logger.info(
                    "[INGEST][DATASET] Completed batch %d/%d for dataset %s",
                    batch_num,
                    total_batches,
                    page_id,
                )
    except Exception as e:
        logger.error(
            "[INGEST][DATASET] Error upserting vectors for %s: %r",
            page_id,
            e,
        )
        raise

    # Update Notion page with embedding metadata after successful upsert
    try:
        _update_notion_page_with_embeddings(
            page_id=canonical_page_id,
            embedding_ids=embedding_ids,
            embedding_count=len(embedding_ids),
        )
    except Exception as e:
        # Log but don't raise - Pinecone ingestion succeeded
        logger.warning(
            "[INGEST][DATASET] Failed to update Notion page metadata for %s (ingestion succeeded): %r",
            canonical_page_id,
            e,
        )

    # Extract mwTab data once for both metadata and feature extraction
    mwtab_data: Optional[Dict[str, Any]] = None
    try:
        mwtab_data = _extract_mwtab_from_page_content(full_text)

        # Fallback: If extraction from page content failed, try fetching from MW API
        if not mwtab_data:
            logger.info(
                "[INGEST][MWTAB] Using fallback MW API fetch (page content extraction failed)"
            )

            # Try to extract STUDY_ID from multiple sources
            import re

            study_id = None

            # 1. Check Summary property
            summary_prop = props.get("Summary", {})
            summary_rich = summary_prop.get("rich_text", []) or []
            summary_text = "".join(rt.get("plain_text", "") for rt in summary_rich)
            study_id_match = re.search(r"ST\d{6}", summary_text)
            if study_id_match:
                study_id = study_id_match.group(0)

            # 2. Check Title/Name property
            if not study_id:
                name_prop = (
                    props.get("Name", {}) or props.get("Experiment Name", {}) or {}
                )
                if name_prop.get("title"):
                    title_text = "".join(
                        t.get("plain_text", "") for t in name_prop["title"]
                    )
                    study_id_match = re.search(r"ST\d{6}", title_text)
                    if study_id_match:
                        study_id = study_id_match.group(0)

            # 3. Check Source URL / DOI property
            if not study_id:
                url_prop = props.get("Source URL / DOI", {}) or {}
                url_text = url_prop.get("url", "") or ""
                study_id_match = re.search(r"ST\d{6}", url_text)
                if study_id_match:
                    study_id = study_id_match.group(0)

            # 4. Check page content itself
            if not study_id:
                study_id_match = re.search(r"ST\d{6}", full_text)
                if study_id_match:
                    study_id = study_id_match.group(0)

            if study_id:
                logger.info(
                    "[INGEST][MWTAB] Attempting MW API fallback fetch for STUDY_ID: %s",
                    study_id,
                )
                try:
                    # Import fetch function from harvest script
                    import sys
                    from pathlib import Path

                    scripts_path = Path(__file__).parent.parent.parent / "scripts"
                    if str(scripts_path) not in sys.path:
                        sys.path.insert(0, str(scripts_path))
                    from harvest_mw_studies import fetch_mw_mwtab

                    mwtab_text = fetch_mw_mwtab(study_id)
                    if mwtab_text:
                        logger.debug(
                            "[INGEST][MWTAB] Fetched mwTab text from MW API (%d chars)",
                            len(mwtab_text),
                        )

                        # Parse the fetched mwTab JSON
                        first_brace = mwtab_text.find("{")
                        if first_brace >= 0:
                            # Try parsing from first brace
                            remaining = mwtab_text[first_brace:]
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
                                        break
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
                                    except (
                                        json.JSONDecodeError,
                                        ValueError,
                                    ) as parse_err:
                                        logger.warning(
                                            "[INGEST][MWTAB] Failed to parse mwTab JSON from MW API response: %r",
                                            str(parse_err)[:200],
                                        )
                        else:
                            logger.warning(
                                "[INGEST][MWTAB] MW API response does not contain JSON (no '{' found)"
                            )
                    else:
                        logger.warning(
                            "[INGEST][MWTAB] MW API fallback fetch returned empty response"
                        )

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
            else:
                logger.info(
                    "[INGEST][MWTAB] Could not extract STUDY_ID from page properties or content for MW API fallback"
                )
    except Exception as e:
        logger.warning(
            "[INGEST][DATASET] Error extracting mwTab data for %s: %r",
            page_id,
            e,
        )
        mwtab_data = None

    # Final check: if mwTab data is still None, log explicit error
    if not mwtab_data:
        logger.warning(
            "[INGEST][MWTAB] ERROR: Could not extract mwTab for dataset %s. Signature matching will be skipped.",
            page_id,
        )

    # Extract and update scientific metadata from mwTab data
    if mwtab_data:
        try:
            logger.info("[INGEST][DATASET] Extracted mwTab data for %s", page_id)
            metadata = _extract_metadata_from_mwtab(mwtab_data)

            # Only update if we have at least one field to set
            has_metadata = any(
                [
                    metadata.get("model_systems"),
                    metadata.get("disease_terms"),
                    metadata.get("matrix_terms"),
                    metadata.get("methods"),
                    metadata.get("summary"),
                    metadata.get("results"),
                    metadata.get("conclusions"),
                    metadata.get("data_origin"),
                    metadata.get("dataset_source_type"),
                    metadata.get("source_url"),
                ]
            )

            if has_metadata:
                _update_experimental_data_asset_metadata(
                    page_id=canonical_page_id,
                    metadata=metadata,
                )
            else:
                logger.info(
                    "[INGEST][DATASET] No extractable metadata found in mwTab for %s",
                    page_id,
                )
        except Exception as e:
            logger.warning(
                "[INGEST][DATASET] Error extracting/updating scientific metadata for %s: %r",
                page_id,
                e,
            )

    # Extract and link metabolite features from mwTab
    if mwtab_data:
        try:
            feature_names = extract_features_from_mwtab(mwtab_data)
            if feature_names:
                link_features_to_notion_items(
                    feature_names=feature_names,
                    item_page_id=canonical_page_id,
                    item_type="dataset",
                )
            else:
                logger.debug(
                    "[INGEST][DATASET] No metabolite features extracted from mwTab for %s",
                    page_id,
                )
        except Exception as e:
            logger.warning(
                "[INGEST][DATASET] Error extracting/linking features for %s: %r",
                page_id,
                e,
            )

    # Detect and ingest signatures from dataset content
    try:
        # Combine page text and mwTab JSON as text for signature detection
        all_text_content = full_text
        if mwtab_data:
            # Convert mwTab JSON to text representation for signature detection
            mwtab_text = json.dumps(mwtab_data, indent=2)
            all_text_content = f"{full_text}\n\n{mwtab_text}"

        # Build source metadata for signature inference
        source_metadata = {
            "diseases": base_meta.get("diseases", []),
            "matrix": base_meta.get("matrix", []),
            "model_systems": base_meta.get("model_systems", []),
        }

        # No attached files for datasets typically (mwTab is in page content)
        attachment_paths: List[Path] = []

        detect_and_ingest_signatures_from_content(
            all_text_content=all_text_content,
            attachment_paths=attachment_paths,
            source_page_id=canonical_page_id,
            source_type="dataset",
            source_metadata=source_metadata,
            source_name=dataset_title,
        )
    except Exception as e:
        logger.warning(
            "[INGEST][DATASET] Error detecting/ingesting signatures for %s: %r",
            page_id,
            e,
        )
        # Non-blocking - continue

    # Match dataset against existing signatures
    cfg = get_config()
    if cfg.pipeline.enable_signature_scoring:
        logger.info(
            "[INGEST][SIGNATURE-MATCH] Signature scoring enabled, mwTab data present: %s",
            mwtab_data is not None,
        )
    if cfg.pipeline.enable_signature_scoring and mwtab_data:
        try:
            # Extract dataset species from mwTab metabolite data
            dataset_species_set: set[str] = set()

            # Extract from mwTab MS_METABOLITE_DATA if available
            metabolite_sections = [
                "MS_METABOLITE_DATA",
                "GC_METABOLITE_DATA",
                "LC_METABOLITE_DATA",
                "METABOLITE_DATA",
            ]
            for section_key in metabolite_sections:
                if section_key in mwtab_data:
                    section = mwtab_data.get(section_key, {})
                    data_array = section.get("Data", [])
                    if isinstance(data_array, list):
                        for row in data_array:
                            if isinstance(row, dict):
                                # Look for metabolite name keys
                                for key in row.keys():
                                    if key.lower() in [
                                        "metabolite",
                                        "metabolite_name",
                                        "compound",
                                        "name",
                                    ]:
                                        raw_name = row.get(key)
                                        if raw_name and isinstance(raw_name, str):
                                            raw_name = raw_name.strip()
                                            if raw_name:
                                                # Map to canonical species if possible
                                                canonical = (
                                                    map_raw_lipid_to_canonical_species(
                                                        raw_name
                                                    )
                                                )
                                                dataset_species_set.add(
                                                    canonical if canonical else raw_name
                                                )

            # Find matching signatures
            if dataset_species_set:
                logger.info(
                    "[INGEST][SIGNATURE-MATCH] Matching dataset %s against signatures (%d species)",
                    page_id,
                    len(dataset_species_set),
                )

                matches = find_matching_signatures_for_dataset(
                    dataset_species=dataset_species_set,
                    overlap_threshold=cfg.pipeline.signature_overlap_threshold,
                )

                if matches:
                    logger.info(
                        "[INGEST][SIGNATURE-MATCH] Found %d matching signature(s) for dataset %s",
                        len(matches),
                        page_id,
                    )

                    # Update dataset page with matches
                    update_dataset_with_signature_matches(
                        dataset_page_id=canonical_page_id,
                        matches=matches,
                    )
                else:
                    logger.info(
                        "[INGEST][SIGNATURE-MATCH] No signature matches found for dataset %s",
                        page_id,
                    )
            else:
                logger.info(
                    "[INGEST][SIGNATURE-MATCH] No species extracted from dataset %s for matching",
                    page_id,
                )

        except Exception as e:
            logger.warning(
                "[INGEST][SIGNATURE-MATCH] Error matching signatures for dataset %s: %r",
                page_id,
                e,
            )
            # Non-blocking - continue

    logger.info("[INGEST][DATASET] Ingestion complete for dataset %s", page_id)
