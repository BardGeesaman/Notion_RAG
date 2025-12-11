"""Semantic cache for RAG query results."""
import time
import logging
from typing import Optional, Dict, Any, Tuple
import numpy as np

logger = logging.getLogger(__name__)


class SemanticCache:
    """Cache RAG results keyed by query embedding similarity."""
    
    def __init__(self, ttl_seconds: int = 3600, similarity_threshold: float = 0.92):
        self.ttl_seconds = ttl_seconds
        self.similarity_threshold = similarity_threshold
        self._cache: Dict[str, Tuple[np.ndarray, Any, float]] = {}  # key -> (embedding, result, timestamp)
    
    def _get_embedding(self, query: str) -> np.ndarray:
        """Get embedding for query using OpenAI."""
        from amprenta_rag.clients.openai_client import get_openai_client
        client = get_openai_client()
        response = client.embeddings.create(
            model="text-embedding-3-small",
            input=query
        )
        return np.array(response.data[0].embedding)
    
    def _cosine_similarity(self, a: np.ndarray, b: np.ndarray) -> float:
        """Compute cosine similarity between two vectors."""
        return float(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b)))
    
    def get(self, query: str) -> Optional[Any]:
        """Find cached result for semantically similar query."""
        query_embedding = self._get_embedding(query)
        now = time.time()
        
        for cache_key, (cached_embedding, result, timestamp) in list(self._cache.items()):
            # Check TTL
            if now - timestamp > self.ttl_seconds:
                del self._cache[cache_key]
                continue
            
            similarity = self._cosine_similarity(query_embedding, cached_embedding)
            if similarity >= self.similarity_threshold:
                logger.info("[CACHE] Hit (similarity=%.3f) for query: %s", similarity, query[:50])
                return result
        
        logger.info("[CACHE] Miss for query: %s", query[:50])
        return None
    
    def set(self, query: str, result: Any) -> None:
        """Store result in cache."""
        embedding = self._get_embedding(query)
        self._cache[query] = (embedding, result, time.time())
        logger.info("[CACHE] Stored result for query: %s", query[:50])
    
    def clear(self) -> None:
        """Clear all cached entries."""
        self._cache.clear()
        logger.info("[CACHE] Cleared")


# Module-level singleton
_cache_instance: Optional[SemanticCache] = None

def get_semantic_cache(ttl_seconds: int = 3600, similarity_threshold: float = 0.92) -> SemanticCache:
    """Get or create the semantic cache singleton."""
    global _cache_instance
    if _cache_instance is None:
        _cache_instance = SemanticCache(ttl_seconds, similarity_threshold)
    return _cache_instance
