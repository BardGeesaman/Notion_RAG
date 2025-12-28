"""
RAG embedding for transcriptomics datasets.

This module provides functions for embedding transcriptomics datasets into Pinecone
for semantic search and retrieval.
"""

from __future__ import annotations

from typing import List, Optional, Set

import pandas as pd

from amprenta_rag.clients.vector_store import get_vector_store
from amprenta_rag.config import get_config
# DEPRECATED: Notion imports removed - Postgres is now source of truth
# from amprenta_rag.ingestion.dataset_notion_utils import update_dataset_embedding_metadata
from amprenta_rag.utils.metadata import sanitize_metadata
from amprenta_rag.ingestion.text_embedding_utils import chunk_text, embed_texts
from amprenta_rag.ingestion.transcriptomics.text_building import build_dge_text_representation
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def embed_transcriptomics_dataset(
    page_id: str,
    dataset_name: str,
    genes: Set[str],
    df: pd.DataFrame,
    gene_column: str,
    program_ids: Optional[List[str]] = None,
    experiment_ids: Optional[List[str]] = None,
) -> None:
    """
    Embed transcriptomics dataset into Pinecone for RAG.

    Args:
        page_id: Notion page ID
        dataset_name: Dataset name
        genes: Set of normalized genes
        df: Full DataFrame with DGE data
        gene_column: Name of the gene column
        program_ids: Optional list of Program page IDs
        experiment_ids: Optional list of Experiment page IDs
    """
    cfg = get_config()

    try:
        # Build text representation
        dataset_text = build_dge_text_representation(
            dataset_name=dataset_name,
            genes=genes,
            df=df,
            gene_column=gene_column,
            program_ids=program_ids,
            experiment_ids=experiment_ids,
        )

        # Chunk and embed
        chunks = chunk_text(dataset_text, max_chars=2000)
        if not chunks:
            logger.warning(
                "[INGEST][TRANSCRIPTOMICS] No chunks generated for dataset %s",
                page_id,
            )
            return

        embeddings = embed_texts(chunks)

        # Upsert to vector store
        store = get_vector_store()
        vectors = []
        for order, (chunk, emb) in enumerate(zip(chunks, embeddings)):
            chunk_id = f"{page_id.replace('-', '')}_trans_chunk_{order:03d}"

            meta = {
                "source": "Dataset",
                "source_type": "Dataset",
                "omics_type": "Transcriptomics",
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
            "[INGEST][TRANSCRIPTOMICS] Generated %d chunk(s) for transcriptomics dataset %s",
            len(chunks),
            page_id,
        )
        logger.info(
            "[INGEST][TRANSCRIPTOMICS] Upserted %d vectors to Pinecone",
            len(vectors),
        )
        logger.info(
            "[INGEST][TRANSCRIPTOMICS] Updated Embedding IDs and Last Embedded for dataset %s",
            page_id,
        )

    except Exception as e:
        logger.warning(
            "[INGEST][TRANSCRIPTOMICS] Error embedding dataset %s: %r",
            page_id,
            e,
        )
        # Don't raise - embedding is non-critical

