"""Cache management utilities for admin operations."""

from typing import Dict, Any, List
import time
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def get_cache_stats() -> Dict[str, Dict[str, Any]]:
    """Get statistics for all caches."""
    stats = {}
    
    # Semantic Cache
    try:
        from amprenta_rag.query.semantic_cache import get_semantic_cache
        semantic_cache = get_semantic_cache()
        stats["semantic_cache"] = {
            "type": "SemanticCache",
            "entries": len(semantic_cache._cache),
            "ttl_seconds": semantic_cache.ttl_seconds,
            "similarity_threshold": semantic_cache.similarity_threshold,
        }
    except Exception as e:
        logger.warning("Failed to get semantic cache stats: %s", e)
        stats["semantic_cache"] = {"type": "SemanticCache", "error": str(e)}
    
    # Dataset Feature Cache
    try:
        from amprenta_rag.ingestion.dataset_feature_cache import get_feature_cache
        feature_cache = get_feature_cache()
        stats["dataset_feature_cache"] = {
            "type": "DatasetFeatureCache",
            "entries": len(feature_cache.cache),
            "ttl_seconds": feature_cache.default_ttl,
            "max_size": feature_cache.max_size,
            "hits": feature_cache.hits,
            "misses": feature_cache.misses,
            "evictions": feature_cache.evictions,
        }
    except Exception as e:
        logger.warning("Failed to get dataset feature cache stats: %s", e)
        stats["dataset_feature_cache"] = {"type": "DatasetFeatureCache", "error": str(e)}
    
    # Enhanced Feature Cache
    try:
        from amprenta_rag.ingestion.enhanced_feature_cache import get_enhanced_feature_cache
        enhanced_cache = get_enhanced_feature_cache()
        stats["enhanced_feature_cache"] = {
            "type": "EnhancedDatasetFeatureCache",
            "entries": len(enhanced_cache.cache),
            "default_ttl_seconds": enhanced_cache.default_ttl.total_seconds(),
            "max_size": enhanced_cache.max_size,
            "hits": enhanced_cache.hits,
            "misses": enhanced_cache.misses,
            "evictions": enhanced_cache.evictions,
            "enable_persistence": enhanced_cache.enable_persistence,
        }
    except Exception as e:
        logger.warning("Failed to get enhanced feature cache stats: %s", e)
        stats["enhanced_feature_cache"] = {"type": "EnhancedDatasetFeatureCache", "error": str(e)}
    
    return stats


def clear_cache(cache_name: str) -> bool:
    """Clear a specific cache by name."""
    try:
        if cache_name == "semantic_cache":
            from amprenta_rag.query.semantic_cache import get_semantic_cache
            get_semantic_cache().clear()
            logger.info("Cleared semantic cache")
            return True
        
        elif cache_name == "dataset_feature_cache":
            from amprenta_rag.ingestion.dataset_feature_cache import clear_feature_cache
            clear_feature_cache()
            logger.info("Cleared dataset feature cache")
            return True
        
        elif cache_name == "enhanced_feature_cache":
            from amprenta_rag.ingestion.enhanced_feature_cache import get_enhanced_feature_cache
            enhanced_cache = get_enhanced_feature_cache()
            enhanced_cache.clear()
            logger.info("Cleared enhanced feature cache")
            return True
        
        else:
            logger.warning("Unknown cache name: %s", cache_name)
            return False
            
    except Exception as e:
        logger.error("Failed to clear cache %s: %s", cache_name, e)
        return False


def clear_all_caches() -> Dict[str, bool]:
    """Clear all caches."""
    results = {}
    cache_names = ["semantic_cache", "dataset_feature_cache", "enhanced_feature_cache"]
    
    for cache_name in cache_names:
        results[cache_name] = clear_cache(cache_name)
    
    logger.info("Cleared all caches: %s", results)
    return results


def get_cache_summary() -> Dict[str, Any]:
    """Get a high-level summary of cache usage."""
    stats = get_cache_stats()
    
    summary = {
        "total_caches": len(stats),
        "total_entries": 0,
        "caches_with_errors": [],
        "cache_names": list(stats.keys()),
    }
    
    for cache_name, cache_stats in stats.items():
        if "error" in cache_stats:
            summary["caches_with_errors"].append(cache_name)
        else:
            summary["total_entries"] += cache_stats.get("entries", 0)
    
    return summary
