"""Cross-encoder reranker for RAG chunks."""
from __future__ import annotations

import time
from typing import List, Dict, Any

from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)

try:
    from sentence_transformers import CrossEncoder
except ImportError:  # pragma: no cover
    CrossEncoder = None  # type: ignore


class CrossEncoderReranker:
    def __init__(self, model_name: str = "cross-encoder/ms-marco-MiniLM-L-6-v2") -> None:
        self.model_name = model_name
        self._model: Any = None

    def _load_model(self):
        if self._model is None:
            if CrossEncoder is None:
                raise ImportError("sentence_transformers not installed")
            self._model = CrossEncoder(self.model_name)

    def rerank(self, query: str, documents: List[Dict[str, Any]], top_k: int = 10) -> List[Dict[str, Any]]:
        if not documents or not query:
            return []
        self._load_model()

        pairs = []
        texts = []
        ids = []
        for doc in documents:
            doc_text = doc.get("chunk_text") or doc.get("content") or ""
            pairs.append((query, doc_text))
            texts.append(doc_text)
            ids.append(doc.get("id") or doc.get("chunk_id"))

        start = time.time()
        scores = self._model.predict(pairs)
        logger.debug("[RERANK] Scored %d docs in %.3f s", len(documents), time.time() - start)

        # Attach scores and sort
        reranked = []
        for doc, score in zip(documents, scores):
            out = dict(doc)
            out["rerank_score"] = float(score)
            reranked.append(out)

        reranked.sort(key=lambda x: x.get("rerank_score", 0), reverse=True)
        return reranked[:top_k]


_singleton: CrossEncoderReranker | None = None


def get_reranker() -> CrossEncoderReranker:
    global _singleton
    if _singleton is None:
        _singleton = CrossEncoderReranker()
    return _singleton
