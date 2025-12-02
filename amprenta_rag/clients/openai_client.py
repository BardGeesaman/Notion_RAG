# amprenta_rag/clients/openai_client.py

from __future__ import annotations

from openai import OpenAI

from amprenta_rag.config import get_config


_client: OpenAI | None = None


def get_openai_client() -> OpenAI:
    global _client
    if _client is None:
        cfg = get_config().openai
        _client = OpenAI(api_key=cfg.api_key)
    return _client


def get_default_models() -> tuple[str, str]:
    cfg = get_config().openai
    return cfg.model, cfg.embedding_model