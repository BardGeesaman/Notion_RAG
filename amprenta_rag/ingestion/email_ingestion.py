# amprenta_rag/ingestion/email_ingestion.py

"""
Email & Note ingestion into the Amprenta RAG engine.

Refactored from backup_2025_11_30/email_ingestion.py to:
- Use amprenta_rag.config + clients (OpenAI, Pinecone, Notion)
- Reuse the existing RAG DB schema:
    - Chunk ID (title)
    - Chunk Text (rich_text)
    - Parent Type (select)
    - Parent Email/Note Item (relation) for Email/Note
    - Parent Item (relation) for other types
    - Order (number)
    - Embedding Vector ID (rich_text)
    - Last Embedded (date)
- Operate on Notion Email DB where:
    - "Embedding Status" = "Not Embedded"
"""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional

from amprenta_rag.clients.pinecone_client import get_pinecone_index
from amprenta_rag.config import get_config
from amprenta_rag.ingestion.feature_extraction import (
    extract_features_from_text, link_features_to_notion_items)
from amprenta_rag.ingestion.metadata_semantic import \
    get_email_semantic_metadata
from amprenta_rag.ingestion.email_cleanup import (
    cleanup_orphaned_chunks, delete_email_and_chunks)
# DEPRECATED: Notion imports removed - Postgres is now source of truth
# Stub functions for backward compatibility
def extract_page_content(page_id: str) -> str:
    """DEPRECATED: Notion support removed. Returns empty string."""
    logger.debug("[INGEST][EMAIL] extract_page_content() deprecated - Notion support removed")
    return ""

def create_rag_chunk_page(chunk_id: str, chunk_text: str, parent_type: str, parent_id: str, order: int, when_iso: str) -> Optional[str]:
    """DEPRECATED: Notion support removed. Returns None."""
    logger.debug("[INGEST][EMAIL] create_rag_chunk_page() deprecated - Notion support removed")
    return None

def fetch_not_embedded_emails() -> List[Dict[str, Any]]:
    """DEPRECATED: Notion support removed. Returns empty list."""
    logger.debug("[INGEST][EMAIL] fetch_not_embedded_emails() deprecated - Notion support removed")
    return []

def update_email_page(page_id: str, when_iso: str) -> None:
    """DEPRECATED: Notion support removed. Does nothing."""
    logger.debug("[INGEST][EMAIL] update_email_page() deprecated - Notion support removed")
    return

from amprenta_rag.ingestion.pinecone_utils import sanitize_metadata
from amprenta_rag.ingestion.signature_integration import \
    detect_and_ingest_signatures_from_content
from amprenta_rag.ingestion.text_embedding_utils import chunk_text, embed_texts
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)

# Backward compatibility wrappers for cleanup functions moved to email_cleanup.py
def cleanup_orphaned_chunks() -> None:
    """Backward compatibility wrapper for cleanup_orphaned_chunks."""
    from amprenta_rag.ingestion.email_cleanup import cleanup_orphaned_chunks as cleanup_func
    return cleanup_func()


def delete_email_and_chunks(email_page_id: str) -> None:
    """Backward compatibility wrapper for delete_email_and_chunks."""
    from amprenta_rag.ingestion.email_cleanup import delete_email_and_chunks as delete_func
    return delete_func(email_page_id)


# ----------------- Main ingestion ----------------- #


def ingest_email(email_page: Dict[str, Any], parent_type: str = "Email") -> None:
    """
    Ingest a single email/note page into RAG Engine + Pinecone.

    This function:
    1. Extracts metadata (title, from, tags, semantic metadata) from the email page
    2. Extracts full text content from Notion blocks
    3. Chunks the text using shared chunking logic
    4. Embeds chunks using OpenAI
    5. Creates RAG chunk pages in Notion
    6. Upserts vectors to Pinecone with rich metadata
    7. Updates the email page status to "Embedded"

    Metadata structure matches zotero_ingest.py for consistency:
    - Doc-level fields: doc_id, doc_source, doc_type, diseases, targets, etc.
    - Lipid-level fields: lipids, lipid_classes, lipid_signatures, etc.
    - Chunk-level fields: chunk_id, chunk_index, snippet, etc.

    Args:
        email_page: Notion page object from Email DB query
        parent_type: Type label ("Email" or "Note") for the parent
    """
    index = get_pinecone_index()
    cfg = get_config()

    page_id = email_page["id"].replace("-", "")
    props = email_page.get("properties", {})

    # Extract metadata
    title_prop = props.get("Title", {}).get("title", [])
    title = title_prop[0]["plain_text"] if title_prop else "Untitled"

    from_prop = props.get("From", {}).get("rich_text", [])
    from_sender = from_prop[0]["plain_text"] if from_prop else ""

    # Extract tags from multi-select property "Tags" (if present)
    tags_prop = props.get("Tags", {})
    tags: List[str] = []
    if "multi_select" in tags_prop:
        tags = [t["name"] for t in tags_prop["multi_select"] if t.get("name")]

    type_prop = props.get("Type", {}).get("select", {})
    item_type = type_prop.get("name", "Email")

    # Doc-level + lipid-level metadata from Email DB
    semantic_meta = get_email_semantic_metadata(email_page)

    # Identity / doc-level fields
    doc_meta: Dict[str, Any] = {
        "doc_id": f"EMAIL:{page_id}",
        "doc_source": "Email",
        "doc_source_subtype": "EmailBody",
        "doc_type": "InternalNote" if item_type == "Note" else "Email",
        "diseases": semantic_meta.get("diseases", []),
        "targets": semantic_meta.get("targets", []),
        "modality": semantic_meta.get("modality", []),
        "stage": semantic_meta.get("stage"),
        "model_systems": semantic_meta.get("model_systems", []),
        "biomarker_role": semantic_meta.get("biomarker_role", []),
        "importance": semantic_meta.get("importance"),
        "manual_tags": tags,
        # lipid-level
        "lipids_raw": semantic_meta.get("lipids_raw", []),
        "lipids": semantic_meta.get("lipids", []),
        "lipid_classes": semantic_meta.get("lipid_classes", []),
        "lipid_signatures": semantic_meta.get("lipid_signatures", []),
        "lipid_signature_role": semantic_meta.get("lipid_signature_role", []),
        "phenotype_axes": semantic_meta.get("phenotype_axes", []),
        "signature_ownership": semantic_meta.get("signature_ownership", []),
        "matrix": semantic_meta.get("matrix", []),
        "treatment_arms": semantic_meta.get("treatment_arms", []),
    }

    logger.info("[INGEST][EMAIL] Processing: %s (%s)", title[:60], item_type)

    # Extract full content from blocks
    try:
        full_text = extract_page_content(page_id)
    except Exception as e:
        logger.error(
            "[INGEST][EMAIL] Error extracting page content for email %s: %r",
            page_id,
            e,
        )
        raise

    # Add metadata header for more context
    header = f"Title: {title}\n"
    if from_sender:
        header += f"From: {from_sender}\n"
    header += f"Type: {item_type}\n\n"

    full_text = header + full_text

    if len(full_text.strip()) < 50:
        logger.info(
            "[INGEST][EMAIL] Skipping - content too short (%d chars) for email %s",
            len(full_text),
            page_id,
        )
        return

    logger.info(
        "[INGEST][EMAIL] Extracted %d characters for email %s", len(full_text), page_id
    )

    # Chunk and embed using shared helpers from zotero_ingest
    chunks = chunk_text(full_text)
    logger.info(
        "[INGEST][EMAIL] Generated %d chunks for email %s", len(chunks), page_id
    )

    if not chunks:
        logger.info("[INGEST][EMAIL] No chunks to embed for email %s", page_id)
        return

    logger.info("[INGEST][EMAIL] Embedding chunks with OpenAI for email %s", page_id)
    try:
        embeddings = embed_texts(chunks)
    except Exception as e:
        logger.error(
            "[INGEST][EMAIL] OpenAI API error embedding chunks for email %s: %r",
            page_id,
            e,
        )
        raise

    now = datetime.now(timezone.utc).isoformat()
    vectors = []

    for order, (chunk, emb) in enumerate(zip(chunks, embeddings)):
        chunk_id = f"{page_id}_chunk_{order:03d}"

        # Create Notion chunk page using shared helper from notion_pages
        # The helper automatically handles Email/Note relation property selection
        try:
            chunk_page_id = create_rag_chunk_page(
                chunk_id, chunk, parent_type, page_id, order, when_iso=now
            )
        except Exception as e:
            logger.error(
                "[INGEST][EMAIL] Notion API error creating chunk page %s: %r",
                chunk_id,
                e,
            )
            raise

        snippet = chunk[:200]

        # Build metadata for Pinecone
        metadata: Dict[str, Any] = {
            **doc_meta,  # doc-level + lipid-level from Email DB
            "chunk_id": chunk_id,
            "chunk_index": order,
            "section": "Body",
            "snippet": snippet,
            "notion_chunk_page_id": chunk_page_id,
            "parent_email_page_id": page_id,
            "source": "Email",
            "source_type": parent_type,  # "Email" or "Note"
            "title": title[:200],
            "from": (from_sender[:200] if from_sender else ""),
            "tags": tags,
        }

        vectors.append(
            {
                "id": chunk_id,
                "values": emb,
                "metadata": sanitize_metadata(metadata),
            }
        )

    logger.info(
        "[INGEST][EMAIL] Upserting %d vectors into Pinecone for email %s...",
        len(vectors),
        page_id,
    )
    try:
        index.upsert(vectors=vectors, namespace=cfg.pinecone.namespace)
    except Exception as e:
        logger.error(
            "[INGEST][EMAIL] Pinecone API error upserting vectors for email %s: %r",
            page_id,
            e,
        )
        raise

    logger.info("[INGEST][EMAIL] Updating Email page status for email %s...", page_id)
    try:
        update_email_page(page_id, now)
    except Exception as e:
        logger.error(
            "[INGEST][EMAIL] Notion API error updating email page %s: %r",
            page_id,
            e,
        )
        raise

    # Extract and link metabolite features from email content
    try:
        # Use full_text (includes header) for feature extraction
        feature_names = extract_features_from_text(full_text)
        if feature_names:
            logger.info(
                "[INGEST][EMAIL] Extracted %d metabolite feature(s) from email %s",
                len(feature_names),
                page_id,
            )
            # Use page ID with dashes for Notion API
            page_id_with_dashes = email_page["id"]
            link_features_to_notion_items(
                feature_names=feature_names,
                item_page_id=page_id_with_dashes,
                item_type="email",
            )
        else:
            logger.debug(
                "[INGEST][EMAIL] No metabolite features found in email %s", page_id
            )
    except Exception as e:
        # Log but don't raise - feature extraction is non-critical
        logger.warning(
            "[INGEST][EMAIL] Feature extraction error for email %s: %r",
            page_id,
            e,
        )

    # Detect and ingest signatures from email content
    try:
        # Get source metadata for signature inference
        source_metadata = {
            "diseases": doc_meta.get("diseases", []),
            "matrix": doc_meta.get("matrix", []),
            "model_systems": doc_meta.get("model_systems", []),
        }

        # No attached files for emails typically
        attachment_paths: List[Path] = []

        page_id_with_dashes = email_page["id"]
        detect_and_ingest_signatures_from_content(
            all_text_content=full_text,
            attachment_paths=attachment_paths,
            source_page_id=page_id_with_dashes,
            source_type="email",
            source_metadata=source_metadata,
            source_name=title or f"Email {page_id[:8]}",
        )
    except Exception as e:
        logger.warning(
            "[INGEST][EMAIL] Error detecting/ingesting signatures for email %s: %r",
            page_id,
            e,
        )
        # Non-blocking - continue

    logger.info("[INGEST][EMAIL] Email ingestion complete for email %s", page_id)


def batch_ingest_emails(parent_type: str = "Email") -> None:
    """
    Batch process all not-embedded emails from the Email DB.

    This function:
    1. Cleans up orphaned chunks (chunks whose parent emails no longer exist)
    2. Fetches all emails with "Embedding Status" = "Not Embedded"
    3. Processes each email through ingest_email()
    4. Handles errors gracefully, continuing with remaining emails

    Args:
        parent_type: Type label to use for ingested emails ("Email" or "Note")
    """
    cleanup_orphaned_chunks()

    logger.info("[INGEST][EMAIL] %s", "=" * 60)

    try:
        emails = fetch_not_embedded_emails()
    except Exception as e:
        logger.error("[INGEST][EMAIL] Error fetching not-embedded emails: %r", e)
        raise

    if not emails:
        logger.info("[INGEST][EMAIL] All emails are already embedded!")
        return

    total = len(emails)
    logger.info("[INGEST][EMAIL] Starting batch ingestion of %d email(s)...", total)

    for idx, email_page in enumerate(emails, 1):
        props = email_page.get("properties", {})
        title_prop = props.get("Title", {}).get("title", [])
        title = title_prop[0]["plain_text"] if title_prop else "Untitled"

        logger.info("[INGEST][EMAIL] [%d/%d] %s", idx, total, title[:80])

        try:
            ingest_email(email_page, parent_type=parent_type)
        except Exception as e:
            logger.error(
                "[INGEST][EMAIL] Error processing email [%d/%d] %s: %r",
                idx,
                total,
                title[:80],
                e,
            )

        # Rate limiting
        import time
        time.sleep(0.5)  # Simple rate limiting

    logger.info("[INGEST][EMAIL] Batch complete: processed %d email(s)", total)
