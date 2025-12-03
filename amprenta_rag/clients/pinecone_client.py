# amprenta_rag/clients/pinecone_client.py

from __future__ import annotations

from typing import Any

from pinecone import Pinecone

from amprenta_rag.config import get_config

_pc: Pinecone | None = None
_index: Any | None = None


def get_pinecone_client() -> Pinecone:
    global _pc
    if _pc is None:
        cfg = get_config().pinecone
        _pc = Pinecone(api_key=cfg.api_key)
    return _pc


def get_pinecone_index() -> Any:
    """
    Return the Pinecone index handle using the modern Pinecone client (v8+).
    We avoid importing Index directly because it's not exported as a symbol.
    """
    global _index
    if _index is None:
        cfg = get_config().pinecone
        pc = get_pinecone_client()
        _index = pc.Index(cfg.index_name)
    return _index
