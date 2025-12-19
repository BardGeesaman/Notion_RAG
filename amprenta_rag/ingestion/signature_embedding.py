"""
Signature embedding utilities.

This module handles embedding signatures into Pinecone for RAG queries.
"""

from __future__ import annotations

import textwrap
from typing import Any, Dict, List

from amprenta_rag.clients.pinecone_client import get_pinecone_index
from amprenta_rag.config import get_config
from amprenta_rag.ingestion.pinecone_utils import sanitize_metadata
from amprenta_rag.ingestion.text_embedding_utils import chunk_text, embed_texts
from amprenta_rag.logging_utils import get_logger
from amprenta_rag.signatures.signature_loader import Signature

logger = get_logger(__name__)

__all__ = [
    "embed_signature",
]


def embed_signature(
    signature_page_id: str,
    signature: Signature,
) -> None:
    """
    Embed a signature into Pinecone for RAG queries.

    Creates text representation of signature, chunks it, embeds, and upserts to Pinecone.

    Args:
        signature_page_id: Notion page ID of signature (with dashes)
        signature: Signature object with components
    """
    cfg = get_config()

    try:
        # Build text representation of signature (multi-omics support)
        signature_type = "Multi-Omics Signature" if (signature.modalities and len(signature.modalities) > 1) else "Signature"
        text_parts = [
            f"{signature_type}: {signature.name}",
        ]

        if signature.modalities:
            modalities_str = ", ".join(mod.title() for mod in signature.modalities)
            text_parts.append(f"Modalities: {modalities_str}")

        if signature.description:
            text_parts.append(f"Description: {signature.description}")

        text_parts.append("\nComponents:")
        for comp in signature.components:
            feature_name = getattr(comp, "feature_name", comp.species)
            feature_type = getattr(comp, "feature_type", "lipid")

            comp_line = f"- {feature_name} ({feature_type})"
            if comp.direction:
                comp_line += f" [{comp.direction}]"
            if comp.weight:
                comp_line += f" [weight: {comp.weight}]"
            text_parts.append(comp_line)

        signature_text = "\n".join(text_parts)

        # Chunk (1-2 chunks should be enough for a signature)
        chunks = chunk_text(signature_text, max_chars=2000)
        if not chunks:
            logger.debug(
                "[INGEST][SIGNATURES] No chunks generated for signature %s",
                signature.name,
            )
            return

        logger.info(
            "[INGEST][SIGNATURES] Embedding signature '%s' (%d chunk(s))",
            signature.name,
            len(chunks),
        )

        # Embed chunks
        embeddings = embed_texts(chunks)

        # Upsert to Pinecone
        index = get_pinecone_index()

        vectors: List[Dict[str, Any]] = []
        for order, (chunk, emb) in enumerate(zip(chunks, embeddings)):
            chunk_id = f"{signature_page_id.replace('-', '')}_sig_chunk_{order:03d}"

            meta: Dict[str, Any] = {
                "source": "Signature",
                "source_type": "Signature",
                "signature_page_id": signature_page_id.replace("-", ""),
                "signature_name": signature.name,
                "snippet": textwrap.shorten(chunk, width=300),
            }

            vectors.append(
                {
                    "id": chunk_id,
                    "values": emb,
                    "metadata": sanitize_metadata(meta),
                }
            )

        index.upsert(vectors=vectors, namespace=cfg.pinecone.namespace)

        logger.info(
            "[INGEST][SIGNATURES] Embedded signature '%s' to Pinecone (%d vectors)",
            signature.name,
            len(vectors),
        )
    except Exception as e:
        logger.warning(
            "[INGEST][SIGNATURES] Error embedding signature %s: %r",
            signature.name,
            e,
        )
        # Don't raise - embedding is non-critical

