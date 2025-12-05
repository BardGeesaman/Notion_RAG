"""
Dataset Notion utilities.

This module handles Notion operations for dataset pages:
- Fetching dataset pages
- Updating embedding metadata
- Updating scientific metadata
"""

from __future__ import annotations

from datetime import datetime, timezone
from typing import Any, Dict, List

import requests

from amprenta_rag.clients.notion_client import notion_headers
from amprenta_rag.config import get_config
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)

__all__ = [
    "fetch_dataset_page",
    "update_dataset_embedding_metadata",
    "update_dataset_scientific_metadata",
]


def fetch_dataset_page(page_id: str) -> Dict[str, Any]:
    """
    Fetch a Notion dataset page by ID.

    Args:
        page_id: Notion page ID (with or without dashes)

    Returns:
        Notion page JSON response

    Raises:
        Exception: If page fetch fails
    """
    from amprenta_rag.ingestion.metadata.helpers import fetch_notion_page

    try:
        return fetch_notion_page(page_id)
    except Exception as e:
        logger.error("[DATASET][NOTION] Error fetching page %s: %r", page_id, e)
        raise


def update_dataset_embedding_metadata(
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

    # Format embedding IDs as summary (Notion rich_text limit is 2000 chars)
    # Store: count, first few IDs, and note that full list is in Pinecone
    max_chars = 1900  # Leave buffer for formatting
    if len(embedding_ids) == 0:
        embedding_ids_text = "No embeddings created"
    elif len(embedding_ids) <= 10:
        # Small list: show all
        embedding_ids_text = f"Total: {embedding_count} embeddings\n\n" + "\n".join(embedding_ids)
    else:
        # Large list: show summary with first few examples
        sample_ids = embedding_ids[:5]
        sample_text = "\n".join(sample_ids)
        remaining = len(embedding_ids) - 5
        summary = f"Total: {embedding_count} embeddings\n\nFirst 5 IDs:\n{sample_text}\n\n... and {remaining} more (all IDs stored in Pinecone)"
        
        # Truncate if still too long
        if len(summary) > max_chars:
            # Just store count and note
            embedding_ids_text = f"Total: {embedding_count} embeddings created. All IDs stored in Pinecone (namespace: {cfg.pinecone.namespace})."
        else:
            embedding_ids_text = summary

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


def update_dataset_scientific_metadata(
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

