# amprenta_rag/ingestion/experiments_ingestion.py

from __future__ import annotations

from typing import Dict, Any, List
from pathlib import Path

import textwrap

from amprenta_rag.config import get_config
from amprenta_rag.clients.pinecone_client import get_pinecone_index
from amprenta_rag.ingestion.metadata_semantic import get_experiment_semantic_metadata
from amprenta_rag.ingestion.notion_pages import extract_page_content
from amprenta_rag.ingestion.pinecone_utils import sanitize_metadata
from amprenta_rag.ingestion.zotero_ingest import _chunk_text, _embed_texts
from amprenta_rag.ingestion.feature_extraction import (
    extract_features_from_text,
    link_features_to_notion_items,
)
from amprenta_rag.ingestion.signature_integration import detect_and_ingest_signatures_from_content
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def _fetch_notion_page(page_id: str) -> Dict[str, Any]:
    """Fetch a Notion page by ID."""
    from amprenta_rag.ingestion.metadata_semantic import _fetch_notion_page as _fetch
    
    try:
        return _fetch(page_id)
    except Exception as e:
        logger.error("[EXPERIMENT][NOTION] Error fetching page %s: %r", page_id, e)
        raise


def ingest_experiment(exp_page_id: str, parent_type: str = "Experiment") -> None:
    """
    Ingest a single Experiment page from Notion into Pinecone.
    """
    logger.info("[INGEST][EXPERIMENT] Ingesting experiment page %s", exp_page_id)
    
    try:
        page = _fetch_notion_page(exp_page_id)
    except Exception as e:
        logger.error(
            "[INGEST][EXPERIMENT] Error fetching Notion page %s: %r",
            exp_page_id,
            e,
        )
        raise
    
    props = page.get("properties", {}) or {}

    title_prop = props.get("Name", {})
    title_parts = title_prop.get("title", []) or []
    exp_title = title_parts[0].get("plain_text", "").strip() if title_parts else "(untitled experiment)"

    # Extract full content from blocks (using shared helper)
    try:
        # Normalize page ID (remove dashes for extract_page_content)
        page_id_clean = exp_page_id.replace("-", "")
        full_text = extract_page_content(page_id_clean)
    except Exception as e:
        logger.error(
            "[INGEST][EXPERIMENT] Error extracting content for %s: %r",
            exp_page_id,
            e,
        )
        raise
    if not full_text or len(full_text.strip()) < 50:
        logger.info("[INGEST][EXPERIMENT] Experiment %s has very little text; skipping.", exp_page_id)
        return

    # Semantic + signature metadata
    try:
        base_meta = get_experiment_semantic_metadata(page)
        logger.info(
            "[INGEST][EXPERIMENT] Lipid signatures for %s: %r",
            exp_page_id,
            base_meta.get("lipid_signatures"),
        )
    except Exception as e:
        logger.error(
            "[INGEST][EXPERIMENT] Error reading semantic metadata for %s: %r",
            exp_page_id,
            e,
        )
        raise

    # Chunk and embed
    chunks = _chunk_text(full_text)
    if not chunks:
        logger.info("[INGEST][EXPERIMENT] No chunks produced for %s; skipping.", exp_page_id)
        return

    logger.info(
        "[INGEST][EXPERIMENT] Generated %d chunk(s) for experiment %s",
        len(chunks),
        exp_page_id,
    )

    try:
        embeddings = _embed_texts(chunks)
    except Exception as e:
        logger.error(
            "[INGEST][EXPERIMENT] Error embedding chunks for %s: %r",
            exp_page_id,
            e,
        )
        raise
    
    index = get_pinecone_index()
    cfg = get_config()

    vectors: List[Dict[str, Any]] = []
    for order, (chunk, emb) in enumerate(zip(chunks, embeddings)):
        chunk_id = f"{exp_page_id}_chunk_{order:03d}"

        snippet = textwrap.shorten(chunk, width=300)

        meta: Dict[str, Any] = {
            **base_meta,
            "source": "Experiments",
            "source_type": parent_type,
            "experiment_page_id": exp_page_id,
            "title": exp_title,
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
        logger.info("[INGEST][EXPERIMENT] No vectors to upsert for %s; skipping.", exp_page_id)
        return

    logger.info(
        "[INGEST][EXPERIMENT] Upserting %d vectors into Pinecone for experiment %s",
        len(vectors),
        exp_page_id,
    )
    
    try:
        index.upsert(vectors=vectors, namespace=cfg.pinecone.namespace)
    except Exception as e:
        logger.error(
            "[INGEST][EXPERIMENT] Error upserting vectors for %s: %r",
            exp_page_id,
            e,
        )
        raise
    
    # Extract and link metabolite features from experiment content
    try:
        feature_names = extract_features_from_text(full_text)
        if feature_names:
            logger.info(
                "[INGEST][EXPERIMENT] Extracted %d metabolite feature(s) from experiment %s",
                len(feature_names),
                exp_page_id,
            )
            # Use canonical page ID with dashes from fetched page
            canonical_page_id = page.get("id", exp_page_id)
            link_features_to_notion_items(
                feature_names=feature_names,
                item_page_id=canonical_page_id,
                item_type="experiment",
            )
        else:
            logger.debug(
                "[INGEST][EXPERIMENT] No metabolite features found in experiment %s",
                exp_page_id,
            )
    except Exception as e:
        # Log but don't raise - feature extraction is non-critical
        logger.warning(
            "[INGEST][EXPERIMENT] Feature extraction error for experiment %s: %r",
            exp_page_id,
            e,
        )
    
    # Detect and ingest signatures from experiment content
    try:
        # Get source metadata for signature inference
        source_metadata = {
            "diseases": base_meta.get("diseases", []),
            "matrix": base_meta.get("matrix", []),
            "model_systems": base_meta.get("model_systems", []),
        }
        
        # No attached files for experiments typically
        from pathlib import Path
        attachment_paths: List[Path] = []
        
        canonical_page_id = page.get("id", exp_page_id)
        detect_and_ingest_signatures_from_content(
            all_text_content=full_text,
            attachment_paths=attachment_paths,
            source_page_id=canonical_page_id,
            source_type="experiment",
            source_metadata=source_metadata,
            source_name=exp_title or f"Experiment {exp_page_id[:8]}",
        )
    except Exception as e:
        logger.warning(
            "[INGEST][EXPERIMENT] Error detecting/ingesting signatures for experiment %s: %r",
            exp_page_id,
            e,
        )
        # Non-blocking - continue
    
    logger.info("[INGEST][EXPERIMENT] Ingestion complete for experiment %s", exp_page_id)

