# amprenta_rag/ingestion/experiments_ingestion.py
"""
Experiment ingestion module.

Handles ingestion of experiment data into Pinecone.
Extracts experiment content, chunks it, embeds it, and upserts to Pinecone
with full metadata including lipid signatures and metabolite features.
"""

from __future__ import annotations

import textwrap
from typing import Any, Dict, List

from amprenta_rag.clients.vector_store import get_vector_store
from amprenta_rag.config import get_config
from amprenta_rag.ingestion.feature_extraction import (
    extract_features_from_text, link_features_to_notion_items)
from amprenta_rag.ingestion.metadata_semantic import \
    get_experiment_semantic_metadata
from amprenta_rag.ingestion.pinecone_utils import sanitize_metadata
from amprenta_rag.ingestion.signature_integration import \
    detect_and_ingest_signatures_from_content
from amprenta_rag.ingestion.text_embedding_utils import chunk_text, embed_texts
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def extract_page_content(page_id: str) -> str:
    """Stub: Notion support removed. Returns empty string."""
    logger.debug("[INGESTION][EXPERIMENTS] extract_page_content() is a no-op (Notion removed)")
    return ""


def _fetch_notion_page(page_id: str) -> Dict[str, Any]:
    """Fetch a Notion page by ID."""
    try:
        from amprenta_rag.ingestion import metadata_semantic as meta
        fetch_func = getattr(meta, "_fetch_notion_page", None)
    except Exception as e:
        logger.error("[EXPERIMENT][NOTION] Error importing notion fetcher: %r", e)
        fetch_func = None

    if callable(fetch_func):
        try:
            return fetch_func(page_id)
        except Exception as e:
            logger.error("[EXPERIMENT][NOTION] Error fetching page %s: %r", page_id, e)
    return {}


def _build_experiment_text_representation(
    page: Dict[str, Any], full_text: str
) -> str:
    """
    Build a structured text representation of an experiment ELN page.

    Includes properties and body content in a coherent format.

    Args:
        page: Notion page dictionary
        full_text: Full text content from blocks

    Returns:
        Structured text representation
    """
    props = page.get("properties", {}) or {}

    # Extract title
    title_prop = props.get("Name", {})
    title_parts = title_prop.get("title", []) or []
    exp_name = (
        title_parts[0].get("plain_text", "").strip()
        if title_parts
        else "(untitled experiment)"
    )

    text_parts = [f"Experiment: {exp_name}"]

    # Extract properties
    type_prop = props.get("Type", {})
    if type_prop.get("select"):
        text_parts.append(f"Type: {type_prop['select'].get('name', '')}")

    disease_prop = props.get("Disease", {})
    if disease_prop.get("multi_select"):
        diseases = [d.get("name", "") for d in disease_prop["multi_select"]]
        if diseases:
            text_parts.append(f"Disease: {', '.join(diseases)}")

    matrix_prop = props.get("Matrix", {})
    if matrix_prop.get("multi_select"):
        matrices = [m.get("name", "") for m in matrix_prop["multi_select"]]
        if matrices:
            text_parts.append(f"Matrix: {', '.join(matrices)}")

    model_systems_prop = props.get("Model Systems", {})
    if model_systems_prop.get("multi_select"):
        models = [m.get("name", "") for m in model_systems_prop["multi_select"]]
        if models:
            text_parts.append(f"Model Systems: {', '.join(models)}")

    programs_prop = props.get("Related Programs", {})
    if programs_prop.get("relation"):
        program_count = len(programs_prop["relation"])
        if program_count > 0:
            text_parts.append(f"Programs: {program_count} program(s) linked")

    signatures_prop = props.get("Readout Signatures", {})
    if signatures_prop.get("relation"):
        sig_count = len(signatures_prop["relation"])
        if sig_count > 0:
            text_parts.append(f"Readout Signatures: {sig_count} signature(s) linked")

    datasets_prop = props.get("Related Datasets", {})
    if datasets_prop.get("relation"):
        dataset_count = len(datasets_prop["relation"])
        if dataset_count > 0:
            text_parts.append(f"Related Datasets: {dataset_count} dataset(s) linked")

    text_parts.append("")  # Blank line separator

    # Add body content
    if full_text:
        text_parts.append(full_text)

    return "\n".join(text_parts)


def _update_experiment_embedding_metadata(
    page_id: str,
    embedding_ids: List[str],
    embedding_count: int,
) -> None:
    """
    Stub: Notion support removed. No-op.

    Previously updated the Experiment page with embedding metadata in Notion.
    """
    logger.debug(
        "[INGEST][EXPERIMENT] _update_experiment_embedding_metadata() is a no-op (Notion removed)"
    )


def ingest_experiment(exp_page_id: str, parent_type: str = "Experiment", force: bool = False) -> None:
    """
    Ingest a single Experiment page from Notion into Pinecone.

    This function performs the complete experiment ingestion pipeline:
    1. Fetches the experiment page from Notion
    2. Extracts full text content and builds structured representation
    3. Extracts semantic metadata (diseases, matrix, signatures, etc.)
    4. Chunks and embeds the text
    5. Upserts vectors to Pinecone with metadata
    6. Updates Embedding IDs and Last Embedded on Experiment page
    7. Extracts metabolite features and links them
    8. Detects and ingests signatures from content

    Args:
        exp_page_id: Notion page ID (with or without dashes)
        parent_type: Type label for the experiment (default: "Experiment")
        force: If True, bypass any caching logic (for future use)

    Raises:
        Exception: If ingestion fails at any step
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
    exp_title = (
        title_parts[0].get("plain_text", "").strip()
        if title_parts
        else "(untitled experiment)"
    )

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
    # Build structured text representation with properties
    structured_text = _build_experiment_text_representation(page, full_text)

    if not structured_text or len(structured_text.strip()) < 50:
        logger.info(
            "[INGEST][EXPERIMENT] Experiment %s has very little text; skipping.",
            exp_page_id,
        )
        return

    logger.info(
        "[INGEST][EXPERIMENT] Extracted %d characters of text from experiment %s",
        len(structured_text),
        exp_page_id,
    )

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

    # Chunk and embed structured text
    chunks = chunk_text(structured_text)
    if not chunks:
        logger.info(
            "[INGEST][EXPERIMENT] No chunks produced for %s; skipping.", exp_page_id
        )
        return

    logger.info(
        "[INGEST][EXPERIMENT] Generated %d chunk(s) for experiment %s",
        len(chunks),
        exp_page_id,
    )

    try:
        embeddings = embed_texts(chunks)
    except Exception as e:
        logger.error(
            "[INGEST][EXPERIMENT] Error embedding chunks for %s: %r",
            exp_page_id,
            e,
        )
        raise

    cfg = get_config()
    store = get_vector_store()

    vectors: List[Dict[str, Any]] = []
    canonical_page_id = page.get("id", exp_page_id)
    page_id_clean = canonical_page_id.replace("-", "")

    for order, (chunk, emb) in enumerate(zip(chunks, embeddings)):
        chunk_id = f"{page_id_clean}_exp_chunk_{order:03d}"

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
        logger.info(
            "[INGEST][EXPERIMENT] No vectors to upsert for %s; skipping.", exp_page_id
        )
        return

    logger.info(
        "[INGEST][EXPERIMENT] Upserting %d vectors into Pinecone for experiment %s",
        len(vectors),
        exp_page_id,
    )

    try:
        store.upsert(vectors=vectors, namespace=cfg.pinecone.namespace)
    except Exception as e:
        logger.error(
            "[INGEST][EXPERIMENT] Error upserting vectors for %s: %r",
            exp_page_id,
            e,
        )
        raise

    # Update Experiment page with embedding metadata
    canonical_page_id = page.get("id", exp_page_id)
    embedding_ids = [v["id"] for v in vectors]
    try:
        _update_experiment_embedding_metadata(
            page_id=canonical_page_id,
            embedding_ids=embedding_ids,
            embedding_count=len(embedding_ids),
        )
    except Exception as e:
        logger.warning(
            "[INGEST][EXPERIMENT] Error updating embedding metadata (non-critical): %r",
            e,
        )

    logger.info(
        "[INGEST][EXPERIMENT] Created %d chunk(s), upserted to Pinecone",
        len(chunks),
    )
    logger.info(
        "[INGEST][EXPERIMENT] Updated Embedding IDs and Last Embedded for Experiment %s",
        exp_page_id,
    )

    # Extract and link metabolite features from experiment content
    try:
        feature_names = extract_features_from_text(structured_text)
        if feature_names:
            logger.info(
                "[INGEST][EXPERIMENT] Extracted %d metabolite feature(s) from experiment %s",
                len(feature_names),
                exp_page_id,
            )
            if callable(link_features_to_notion_items):
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

        detect_and_ingest_signatures_from_content(
            all_text_content=structured_text,
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

    logger.info(
        "[INGEST][EXPERIMENT] Ingestion complete for experiment %s", exp_page_id
    )
