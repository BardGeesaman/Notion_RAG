# amprenta_rag/ingestion/dataset_ingestion.py

from __future__ import annotations

from typing import Dict, Any, List

import textwrap

from amprenta_rag.config import get_config
from amprenta_rag.clients.pinecone_client import get_pinecone_index
from amprenta_rag.ingestion.metadata_semantic import get_dataset_semantic_metadata
from amprenta_rag.ingestion.notion_pages import extract_page_content
from amprenta_rag.ingestion.pinecone_utils import sanitize_metadata
from amprenta_rag.ingestion.zotero_ingest import _chunk_text, _embed_texts
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def _fetch_notion_page(page_id: str) -> Dict[str, Any]:
    """Fetch a Notion page by ID."""
    from amprenta_rag.ingestion.metadata_semantic import _fetch_notion_page as _fetch
    
    try:
        return _fetch(page_id)
    except Exception as e:
        logger.error("[DATASET][NOTION] Error fetching page %s: %r", page_id, e)
        raise


def ingest_dataset(page_id: str, force: bool = False) -> None:
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
    
    props = page.get("properties", {}) or {}

    # Get title from "Name" property
    title_prop = props.get("Name", {})
    title_parts = title_prop.get("title", []) or []
    dataset_title = title_parts[0].get("plain_text", "").strip() if title_parts else "(untitled dataset)"

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
        logger.info("[INGEST][DATASET] Dataset %s has very little text; skipping.", page_id)
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
    for order, (chunk, emb) in enumerate(zip(chunks, embeddings)):
        chunk_id = f"{page_id.replace('-', '')}_chunk_{order:03d}"

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
    
    try:
        index.upsert(vectors=vectors, namespace=cfg.pinecone.namespace)
    except Exception as e:
        logger.error(
            "[INGEST][DATASET] Error upserting vectors for %s: %r",
            page_id,
            e,
        )
        raise
    
    logger.info("[INGEST][DATASET] Ingestion complete for dataset %s", page_id)

