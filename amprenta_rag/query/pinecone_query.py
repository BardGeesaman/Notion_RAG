# amprenta_rag/query/pinecone_query.py

"""
Low-level Pinecone query operations.

This module handles:
- Embedding user queries
- Building metadata filters
- Executing Pinecone queries
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional

from amprenta_rag.clients.openai_client import (get_default_models,
                                                get_openai_client)
from amprenta_rag.clients.pinecone_client import get_pinecone_index
from amprenta_rag.config import get_config
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def embed_query(text: str) -> List[float]:
    """
    Embed a user query using the configured OpenAI embedding model.

    Args:
        text: User query text to embed

    Returns:
        Embedding vector as a list of floats
    """
    client = get_openai_client()
    _, embedding_model = get_default_models()
    try:
        resp = client.embeddings.create(
            model=embedding_model,
            input=text,
        )
    except Exception as e:
        logger.error("[PINECONE] OpenAI API error embedding query: %r", e)
        raise
    return resp.data[0].embedding


def build_meta_filter(
    disease: Optional[str] = None,
    target: Optional[str] = None,
    lipid: Optional[str] = None,
    signature: Optional[str] = None,
) -> Optional[Dict[str, Any]]:
    """
    Build a Pinecone metadata filter from high-level query constraints.

    Each argument is a SINGLE value (for now). If provided, they are applied
    as "$in" filters against the corresponding metadata fields:
      - diseases
      - targets
      - lipids (canonical IDs) OR lipids_raw (fall back)
      - lipid_signatures

    Args:
        disease: Optional disease name to filter by
        target: Optional molecular target to filter by
        lipid: Optional lipid (canonical ID or raw label) to filter by
        signature: Optional lipid signature to filter by

    Returns:
        Pinecone metadata filter dict, or None if no filters specified
    """
    filt: Dict[str, Any] = {}

    if disease:
        filt["diseases"] = {"$in": [disease]}
    if target:
        filt["targets"] = {"$in": [target]}
    if signature:
        filt["lipid_signatures"] = {"$in": [signature]}

    # For lipid, try canonical IDs first ("lipids"), then fall back to raw names ("lipids_raw")
    if lipid:
        filt["$or"] = [
            {"lipids": {"$in": [lipid]}},
            {"lipids_raw": {"$in": [lipid]}},
        ]

    return filt or None


def query_pinecone(
    user_query: str,
    top_k: int = 10,
    meta_filter: Optional[Dict[str, Any]] = None,
    source_types: Optional[List[str]] = None,
) -> List[Dict[str, Any]]:
    """
    Execute a raw Pinecone query and return the matches list.

    This function:
    1. Embeds the user query
    2. Queries Pinecone with the embedding, top_k, and optional metadata filter
    3. Returns the raw matches list from the API response

    Args:
        user_query: User query text
        top_k: Number of results to retrieve
        meta_filter: Optional Pinecone metadata filter, e.g.:
            {"diseases": {"$in": ["ALS"]}, "lipid_signatures": {"$in": ["ALS-CSF-Core-6Ceramides"]}}
        source_types: Optional list of source types to filter by (e.g., ["Literature", "Experiment"])

    Returns:
        List of match dictionaries with metadata from Pinecone
    """
    index = get_pinecone_index()
    vector = embed_query(user_query)
    cfg = get_config()

    # Build combined filter
    combined_filter: Dict[str, Any] = {}
    if meta_filter:
        combined_filter.update(meta_filter)

    # Add source type filtering if specified
    if source_types:
        combined_filter["source_type"] = {"$in": source_types}

    kwargs: Dict[str, Any] = {
        "vector": vector,
        "top_k": top_k,
        "include_metadata": True,
        "namespace": cfg.pinecone.namespace,  # use same namespace as upserts
    }
    if combined_filter:
        kwargs["filter"] = combined_filter

    try:
        res = index.query(**kwargs)
    except Exception as e:
        logger.error(
            "[PINECONE] Pinecone API error querying (top_k=%d, filter=%s): %r",
            top_k,
            meta_filter,
            e,
        )
        raise
    matches = getattr(res, "matches", None) or res.get("matches", [])
    return matches
