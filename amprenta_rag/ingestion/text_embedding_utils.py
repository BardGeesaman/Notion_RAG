"""
Text chunking and embedding utilities.

This module provides shared functions for text chunking and embedding
that are used across multiple ingestion modules.
"""

from __future__ import annotations

from typing import List

from amprenta_rag.clients.openai_client import get_default_models, get_openai_client
from amprenta_rag.logging_utils import get_logger
from amprenta_rag.utils.error_handling import RetryConfig, retry_with_backoff

logger = get_logger(__name__)

__all__ = [
    "chunk_text",
    "embed_texts",
]


@retry_with_backoff(RetryConfig(max_attempts=3, initial_delay=1.0, max_delay=8.0))
def _embed_batch(client, model: str, batch: List[str]) -> List[List[float]]:
    resp = client.embeddings.create(
        model=model,
        input=batch,
    )
    return [d.embedding for d in resp.data]


def embed_texts(texts: List[str], batch_size: int = 64) -> List[List[float]]:
    """
    Embed texts in batches using the configured embedding model.

    Args:
        texts: List of text strings to embed
        batch_size: Number of texts to embed per batch

    Returns:
        List of embedding vectors (one per input text)
    """
    if not texts:
        return []

    client = get_openai_client()
    _, embed_model = get_default_models()
    all_embeddings: List[List[float]] = []

    for i in range(0, len(texts), batch_size):
        batch = texts[i : i + batch_size]
        try:
            all_embeddings.extend(_embed_batch(client, embed_model, batch))
        except Exception as e:
            logger.error(
                "[INGEST][EMBED] OpenAI embedding error for batch %d-%d: %r",
                i,
                min(i + batch_size, len(texts)),
                e,
            )
            raise

    return all_embeddings


def chunk_text(text: str, max_chars: int = 2000) -> List[str]:
    """
    Paragraph-aware, character-limited chunking.

    Splits text into chunks that respect paragraph boundaries while
    staying within the character limit.

    Args:
        text: Text to chunk
        max_chars: Maximum characters per chunk

    Returns:
        List of text chunks
    """
    text = text.strip()
    if not text:
        return []

    paras = [p.strip() for p in text.split("\n") if p.strip()]
    chunks: List[str] = []
    buf = ""

    for p in paras:
        if buf and len(buf) + len(p) + 1 > max_chars:
            chunks.append(buf.strip())
            buf = p
        else:
            if buf:
                buf += "\n" + p
            else:
                buf = p

    if buf:
        chunks.append(buf.strip())

    return chunks

