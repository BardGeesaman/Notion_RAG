"""
RAG embedding for lipidomics datasets.

This module provides functions for embedding lipidomics datasets into Pinecone
for semantic search and retrieval.
"""

from __future__ import annotations

from typing import Any, List, Optional, Set

from amprenta_rag.clients.vector_store import get_vector_store
from amprenta_rag.config import get_config
# DEPRECATED: Notion imports removed - Postgres is now source of truth
# from amprenta_rag.ingestion.dataset_notion_utils import update_dataset_embedding_metadata
from amprenta_rag.ingestion.pinecone_utils import sanitize_metadata
from amprenta_rag.ingestion.text_embedding_utils import chunk_text, embed_texts
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def embed_lipidomics_dataset(
    page_id: str,
    dataset_name: str,
    species: Set[str],
    signature_matches: Optional[List[Any]] = None,
) -> None:
    """
    Embed lipidomics dataset into Pinecone for RAG.

    Args:
        page_id: Notion page ID
        dataset_name: Dataset name
        species: Set of normalized species
        signature_matches: Optional list of signature match results
    """
    cfg = get_config()

    try:
        # Build text representation
        text_parts = [
            f"Internal Lipidomics Dataset: {dataset_name}",
            "Data Origin: Internal – Amprenta",
            "Dataset Source Type: Processed table",
            "",
            f"Contains {len(species)} normalized lipid species:",
        ]

        # Add species list (truncate if too long)
        species_list = sorted(list(species))
        if len(species_list) > 50:
            text_parts.append(", ".join(species_list[:50]))
            text_parts.append(f"... and {len(species_list) - 50} more species")
        else:
            text_parts.append(", ".join(species_list))

        # Add signature matches if available
        if signature_matches:
            text_parts.append("")
            text_parts.append("Signature Matches:")
            for match in signature_matches[:5]:  # Top 5
                match_name = getattr(match, "signature_name", "Unknown")
                match_score = getattr(match, "score", 0.0)
                match_overlap = getattr(match, "overlap_fraction", 0.0)
                text_parts.append(
                    f"- {match_name}: score {match_score:.3f}, overlap {match_overlap:.2f}"
                )

        dataset_text = "\n".join(text_parts)

        # Chunk and embed
        chunks = chunk_text(dataset_text, max_chars=2000)
        if not chunks:
            logger.warning(
                "[INGEST][LIPIDOMICS] No chunks generated for dataset %s",
                page_id,
            )
            return

        embeddings = embed_texts(chunks)

        # Upsert to vector store
        store = get_vector_store()
        vectors = []
        for order, (chunk, emb) in enumerate(zip(chunks, embeddings)):
            chunk_id = f"{page_id.replace('-', '')}_lipidomics_chunk_{order:03d}"

            meta = {
                "source": "Dataset",
                "source_type": "Dataset",
                "dataset_page_id": page_id.replace("-", ""),
                "dataset_name": dataset_name,
                "data_origin": "Internal – Amprenta",
                "snippet": chunk[:300],
            }

            vectors.append(
                {
                    "id": chunk_id,
                    "values": emb,
                    "metadata": sanitize_metadata(meta),
                }
            )

        # Batch upsert
        batch_size = 100
        for i in range(0, len(vectors), batch_size):
            batch = vectors[i : i + batch_size]
            store.upsert(vectors=batch, namespace=cfg.pinecone.namespace)

        # DEPRECATED: Notion metadata update removed - Postgres is now source of truth
        # Embedding metadata is stored in Pinecone

        logger.info(
            "[INGEST][LIPIDOMICS] Embedded dataset %s to Pinecone (%d vectors)",
            page_id,
            len(vectors),
        )

    except Exception as e:
        logger.warning(
            "[INGEST][LIPIDOMICS] Error embedding dataset %s: %r",
            page_id,
            e,
        )
        # Don't raise - embedding is non-critical

