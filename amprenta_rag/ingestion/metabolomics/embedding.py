"""
RAG embedding for metabolomics datasets.

This module provides functions for embedding metabolomics datasets into Pinecone
for semantic search and retrieval.
"""

from __future__ import annotations

from typing import List, Optional, Set

from amprenta_rag.clients.pinecone_client import get_pinecone_index
from amprenta_rag.config import get_config
# DEPRECATED: Notion imports removed - Postgres is now source of truth
# from amprenta_rag.ingestion.dataset_notion_utils import update_dataset_embedding_metadata
from amprenta_rag.ingestion.pinecone_utils import sanitize_metadata
from amprenta_rag.ingestion.text_embedding_utils import chunk_text, embed_texts
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def embed_metabolomics_dataset(
    page_id: str,
    dataset_name: str,
    metabolites: Set[str],
    program_ids: Optional[List[str]] = None,
    experiment_ids: Optional[List[str]] = None,
) -> None:
    """
    Embed metabolomics dataset into Pinecone for RAG.

    Args:
        page_id: Notion page ID
        dataset_name: Dataset name
        metabolites: Set of normalized metabolites
        program_ids: Optional list of Program page IDs
        experiment_ids: Optional list of Experiment page IDs
    """
    cfg = get_config()

    try:
        # Build text representation
        text_parts = [
            f"Dataset: {dataset_name}",
            "Type: Metabolomics (internal)",
            "Data Origin: Internal – Amprenta",
        ]

        if program_ids:
            text_parts.append(f"Programs: {len(program_ids)} program(s) linked")

        if experiment_ids:
            text_parts.append(f"Experiments: {len(experiment_ids)} experiment(s) linked")

        text_parts.append("")
        text_parts.append(f"Metabolites ({len(metabolites)} unique):")

        # Add metabolite list (truncate if too long)
        metabolite_list = sorted(list(metabolites))
        if len(metabolite_list) > 100:
            text_parts.append(", ".join(metabolite_list[:100]))
            text_parts.append(f"... and {len(metabolite_list) - 100} more metabolites")
        else:
            text_parts.append(", ".join(metabolite_list))

        dataset_text = "\n".join(text_parts)

        # Chunk and embed
        chunks = chunk_text(dataset_text, max_chars=2000)
        if not chunks:
            logger.warning(
                "[INGEST][METABOLOMICS] No chunks generated for dataset %s",
                page_id,
            )
            return

        embeddings = embed_texts(chunks)

        # Upsert to Pinecone
        index = get_pinecone_index()
        vectors = []
        for order, (chunk, emb) in enumerate(zip(chunks, embeddings)):
            chunk_id = f"{page_id.replace('-', '')}_metab_chunk_{order:03d}"

            meta = {
                "source": "Dataset",
                "source_type": "Dataset",
                "omics_type": "Metabolomics",
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
            index.upsert(vectors=batch, namespace=cfg.pinecone.namespace)

        # DEPRECATED: Notion metadata update removed - Postgres is now source of truth
        # Embedding metadata is stored in Pinecone

        logger.info(
            "[INGEST][METABOLOMICS] Generated %d chunk(s)",
            len(chunks),
        )
        logger.info(
            "[INGEST][METABOLOMICS] Upserted %d vectors to Pinecone",
            len(vectors),
        )
        logger.info(
            "[INGEST][METABOLOMICS] Updated Embedding IDs and Last Embedded for dataset %s",
            page_id,
        )

    except Exception as e:
        logger.warning(
            "[INGEST][METABOLOMICS] Error embedding dataset %s: %r",
            page_id,
            e,
        )
        # Don't raise - embedding is non-critical

