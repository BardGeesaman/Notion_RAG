"""LLM-assisted relevance and novelty scoring for datasets and papers."""
from __future__ import annotations

import hashlib
import json
import logging
import time
from dataclasses import dataclass
from typing import Any, Dict, List, Optional

import numpy as np
from sklearn.metrics.pairwise import cosine_similarity

from amprenta_rag.clients.openai_client import get_openai_client
from amprenta_rag.config import get_config

logger = logging.getLogger(__name__)


@dataclass
class RelevanceScore:
    """Result of relevance scoring."""
    
    item_id: str
    overall_score: float  # 0-1 scale
    disease_match: float  # 0-1 scale
    target_overlap: float  # 0-1 scale
    data_quality: float  # 0-1 scale
    explanation: str
    processing_time_seconds: float
    cached: bool = False


@dataclass
class NoveltyScore:
    """Result of novelty scoring."""
    
    item_id: str
    novelty_score: float  # 0-1 scale (1 = most novel)
    max_similarity: float  # Highest similarity to existing items
    most_similar_item_id: Optional[str]  # ID of most similar existing item
    explanation: str  # LLM explanation of WHY it's novel/not novel
    processing_time_seconds: float
    cached: bool = False


@dataclass
class ScoredItem:
    """Combined scoring result for an item."""
    
    item_id: str
    relevance_score: Optional[RelevanceScore] = None
    novelty_score: Optional[NoveltyScore] = None


# Cache for storing scoring results (in-memory for MVP)
_SCORE_CACHE: Dict[str, Dict[str, Any]] = {}
CACHE_TTL_SECONDS = 24 * 60 * 60  # 24 hours


def _get_cache_key(item_id: str, context_hash: str, score_type: str) -> str:
    """Generate cache key for scoring results."""
    return f"{score_type}:{item_id}:{context_hash}"


def _get_from_cache(cache_key: str) -> Optional[Dict[str, Any]]:
    """Get scoring result from cache if not expired."""
    if cache_key not in _SCORE_CACHE:
        return None
    
    cached_data = _SCORE_CACHE[cache_key]
    if time.time() - cached_data["timestamp"] > CACHE_TTL_SECONDS:
        # Expired, remove from cache
        del _SCORE_CACHE[cache_key]
        return None
    
    return cached_data["result"]


def _set_cache(cache_key: str, result: Dict[str, Any]) -> None:
    """Store scoring result in cache."""
    _SCORE_CACHE[cache_key] = {
        "result": result,
        "timestamp": time.time(),
    }


def _hash_context(context: Dict[str, Any]) -> str:
    """Create hash of context for caching."""
    # Sort keys to ensure consistent hashing
    context_str = json.dumps(context, sort_keys=True)
    return hashlib.md5(context_str.encode()).hexdigest()


def score_relevance(
    item: Dict[str, Any],
    context: Dict[str, Any],
    criteria: Optional[Dict[str, Any]] = None,
) -> RelevanceScore:
    """
    Score relevance of an item using LLM with structured output.
    
    Args:
        item: Item to score (dataset, paper, etc.) with fields like title, description, etc.
        context: Research context including diseases of interest, targets, etc.
        criteria: Optional custom criteria weights and thresholds
        
    Returns:
        RelevanceScore with detailed breakdown and explanation
    """
    start_time = time.time()
    item_id = item.get("id", item.get("accession", "unknown"))
    
    # Check cache first
    context_hash = _hash_context(context)
    cache_key = _get_cache_key(item_id, context_hash, "relevance")
    cached_result = _get_from_cache(cache_key)
    
    if cached_result:
        logger.debug(f"[RELEVANCE] Cache hit for item {item_id}")
        cached_result["cached"] = True
        cached_result["processing_time_seconds"] = time.time() - start_time
        return RelevanceScore(**cached_result)
    
    # Default criteria if not provided
    if criteria is None:
        criteria = {
            "disease_weight": 0.4,
            "target_weight": 0.3,
            "data_quality_weight": 0.3,
        }
    
    try:
        cfg = get_config()
        client = get_openai_client()
        
        # Prepare prompt for LLM
        prompt = _build_relevance_prompt(item, context, criteria)
        
        logger.debug(f"[RELEVANCE] Scoring item {item_id} with LLM")
        
        response = client.chat.completions.create(
            model=cfg.openai.model,
            messages=[
                {
                    "role": "system",
                    "content": "You are a biomedical research relevance scoring assistant. Score items based on their relevance to the given research context and return structured JSON.",
                },
                {"role": "user", "content": prompt},
            ],
            temperature=0.1,  # Low temperature for consistent scoring
            max_tokens=800,
            response_format={"type": "json_object"},
        )
        
        response_text = response.choices[0].message.content
        if not response_text:
            raise ValueError("Empty response from LLM")
        
        # Parse JSON response
        llm_result = json.loads(response_text)
        
        # Validate and extract scores
        overall_score = float(llm_result.get("overall_score", 0.0))
        disease_match = float(llm_result.get("disease_match", 0.0))
        target_overlap = float(llm_result.get("target_overlap", 0.0))
        data_quality = float(llm_result.get("data_quality", 0.0))
        explanation = llm_result.get("explanation", "No explanation provided")
        
        # Ensure scores are in valid range
        overall_score = max(0.0, min(1.0, overall_score))
        disease_match = max(0.0, min(1.0, disease_match))
        target_overlap = max(0.0, min(1.0, target_overlap))
        data_quality = max(0.0, min(1.0, data_quality))
        
        processing_time = time.time() - start_time
        
        result = RelevanceScore(
            item_id=item_id,
            overall_score=overall_score,
            disease_match=disease_match,
            target_overlap=target_overlap,
            data_quality=data_quality,
            explanation=explanation,
            processing_time_seconds=processing_time,
            cached=False,
        )
        
        # Cache the result
        cache_data = {
            "item_id": item_id,
            "overall_score": overall_score,
            "disease_match": disease_match,
            "target_overlap": target_overlap,
            "data_quality": data_quality,
            "explanation": explanation,
            "cached": False,
        }
        _set_cache(cache_key, cache_data)
        
        logger.info(f"[RELEVANCE] Scored item {item_id}: {overall_score:.3f} overall ({processing_time:.2f}s)")
        return result
        
    except Exception as e:
        logger.error(f"[RELEVANCE] Error scoring item {item_id}: {e}", exc_info=True)
        processing_time = time.time() - start_time
        
        # Return default low score on error
        return RelevanceScore(
            item_id=item_id,
            overall_score=0.1,
            disease_match=0.0,
            target_overlap=0.0,
            data_quality=0.0,
            explanation=f"Error during scoring: {str(e)}",
            processing_time_seconds=processing_time,
            cached=False,
        )


def score_novelty(
    item: Dict[str, Any],
    existing_items: List[Dict[str, Any]],
) -> NoveltyScore:
    """
    Score novelty using hybrid approach: embeddings for similarity + LLM for explanation.
    
    Args:
        item: Item to score for novelty
        existing_items: List of existing items to compare against
        
    Returns:
        NoveltyScore with similarity-based score and LLM explanation
    """
    start_time = time.time()
    item_id = item.get("id", item.get("accession", "unknown"))
    
    # Check cache (use existing items hash as context)
    existing_items_hash = _hash_context({"existing_items": [i.get("id", str(i)) for i in existing_items]})
    cache_key = _get_cache_key(item_id, existing_items_hash, "novelty")
    cached_result = _get_from_cache(cache_key)
    
    if cached_result:
        logger.debug(f"[NOVELTY] Cache hit for item {item_id}")
        cached_result["cached"] = True
        cached_result["processing_time_seconds"] = time.time() - start_time
        return NoveltyScore(**cached_result)
    
    try:
        # Step 1: Compute embeddings for similarity calculation
        item_embedding = _get_item_embedding(item)
        existing_embeddings = [_get_item_embedding(existing_item) for existing_item in existing_items]
        
        if not existing_embeddings:
            # No existing items to compare against - highly novel
            novelty_score = 1.0
            max_similarity = 0.0
            most_similar_item_id = None
        else:
            # Step 2: Calculate cosine similarities (cheap computation)
            similarities = []
            for existing_embedding in existing_embeddings:
                similarity = cosine_similarity([item_embedding], [existing_embedding])[0][0]
                similarities.append(similarity)
            
            # Find maximum similarity and most similar item
            max_similarity = float(np.max(similarities))
            most_similar_idx = int(np.argmax(similarities))
            most_similar_item_id = existing_items[most_similar_idx].get("id", 
                                                                       existing_items[most_similar_idx].get("accession", "unknown"))
            
            # Score = 1 - max_similarity (higher similarity = lower novelty)
            novelty_score = 1.0 - max_similarity
        
        # Step 3: LLM for explanation only (not score computation)
        explanation = _get_novelty_explanation(item, existing_items, novelty_score, most_similar_item_id)
        
        processing_time = time.time() - start_time
        
        result = NoveltyScore(
            item_id=item_id,
            novelty_score=novelty_score,
            max_similarity=max_similarity,
            most_similar_item_id=most_similar_item_id,
            explanation=explanation,
            processing_time_seconds=processing_time,
            cached=False,
        )
        
        # Cache the result
        cache_data = {
            "item_id": item_id,
            "novelty_score": novelty_score,
            "max_similarity": max_similarity,
            "most_similar_item_id": most_similar_item_id,
            "explanation": explanation,
            "cached": False,
        }
        _set_cache(cache_key, cache_data)
        
        logger.info(f"[NOVELTY] Scored item {item_id}: {novelty_score:.3f} novelty, {max_similarity:.3f} max_sim ({processing_time:.2f}s)")
        return result
        
    except Exception as e:
        logger.error(f"[NOVELTY] Error scoring item {item_id}: {e}", exc_info=True)
        processing_time = time.time() - start_time
        
        # Return default medium novelty on error
        return NoveltyScore(
            item_id=item_id,
            novelty_score=0.5,
            max_similarity=0.5,
            most_similar_item_id=None,
            explanation=f"Error during novelty scoring: {str(e)}",
            processing_time_seconds=processing_time,
            cached=False,
        )


def batch_score(
    items: List[Dict[str, Any]],
    context: Dict[str, Any],
    score_relevance_flag: bool = True,
    score_novelty_flag: bool = True,
) -> List[ScoredItem]:
    """
    Score multiple items in batch for relevance and/or novelty.
    
    Args:
        items: List of items to score
        context: Research context for relevance scoring
        score_relevance_flag: Whether to compute relevance scores
        score_novelty_flag: Whether to compute novelty scores
        
    Returns:
        List of ScoredItem objects with requested scores
    """
    logger.info(f"[BATCH] Starting batch scoring of {len(items)} items (relevance={score_relevance_flag}, novelty={score_novelty_flag})")
    start_time = time.time()
    
    results = []
    
    for i, item in enumerate(items):
        item_id = item.get("id", item.get("accession", f"item_{i}"))
        scored_item = ScoredItem(item_id=item_id)
        
        try:
            # Score relevance if requested
            if score_relevance_flag:
                scored_item.relevance_score = score_relevance(item, context)
                
                # Brief pause to respect API rate limits
                time.sleep(0.1)
            
            # Score novelty if requested
            if score_novelty_flag:
                # For novelty, compare against all other items in the batch
                other_items = [other_item for j, other_item in enumerate(items) if j != i]
                scored_item.novelty_score = score_novelty(item, other_items)
                
                # Brief pause to respect API rate limits
                time.sleep(0.1)
            
            results.append(scored_item)
            
        except Exception as e:
            logger.error(f"[BATCH] Error scoring item {item_id}: {e}")
            # Add item with error scores
            if score_relevance_flag:
                scored_item.relevance_score = RelevanceScore(
                    item_id=item_id,
                    overall_score=0.0,
                    disease_match=0.0,
                    target_overlap=0.0,
                    data_quality=0.0,
                    explanation=f"Batch scoring error: {str(e)}",
                    processing_time_seconds=0.0,
                    cached=False,
                )
            if score_novelty_flag:
                scored_item.novelty_score = NoveltyScore(
                    item_id=item_id,
                    novelty_score=0.5,
                    max_similarity=0.5,
                    most_similar_item_id=None,
                    explanation=f"Batch scoring error: {str(e)}",
                    processing_time_seconds=0.0,
                    cached=False,
                )
            results.append(scored_item)
    
    processing_time = time.time() - start_time
    logger.info(f"[BATCH] Completed batch scoring in {processing_time:.2f}s")
    
    return results


def _build_relevance_prompt(item: Dict[str, Any], context: Dict[str, Any], criteria: Dict[str, Any]) -> str:
    """Build prompt for LLM relevance scoring."""
    prompt = f"""Score the relevance of this research item to the given context.

ITEM TO SCORE:
Title: {item.get('title', 'N/A')}
Description: {item.get('description', 'N/A')}
Species: {item.get('species', 'N/A')}
Assay Type: {item.get('assay_type', 'N/A')}
Sample Count: {item.get('sample_count', 'N/A')}

RESEARCH CONTEXT:
Diseases of Interest: {context.get('diseases', 'N/A')}
Target Molecules: {context.get('targets', 'N/A')}
Preferred Species: {context.get('species', 'N/A')}
Required Assay Types: {context.get('assay_types', 'N/A')}
Minimum Sample Size: {context.get('min_sample_size', 'N/A')}

SCORING CRITERIA:
Disease Match Weight: {criteria.get('disease_weight', 0.4)}
Target Overlap Weight: {criteria.get('target_weight', 0.3)}
Data Quality Weight: {criteria.get('data_quality_weight', 0.3)}

Please score the item on a 0-1 scale for each criterion and provide an overall score.

Return your response as JSON with this exact structure:
{{
  "overall_score": 0.85,
  "disease_match": 0.9,
  "target_overlap": 0.8,
  "data_quality": 0.85,
  "explanation": "This dataset is highly relevant because..."
}}

Consider:
- Disease match: How well do the item's diseases align with the research context?
- Target overlap: How much overlap is there between the item's targets and context targets?
- Data quality: Sample size, assay quality, species relevance, experimental design quality
- Overall score should be a weighted combination based on the criteria weights
"""
    return prompt


def _get_item_embedding(item: Dict[str, Any]) -> np.ndarray:
    """
    Get embedding for an item using OpenAI embeddings API.
    
    For MVP, we'll use a simple text-based approach.
    In production, this could be cached or use a dedicated embedding service.
    """
    try:
        client = get_openai_client()
        
        # Combine key fields into text for embedding
        text_parts = []
        if item.get("title"):
            text_parts.append(f"Title: {item['title']}")
        if item.get("description"):
            text_parts.append(f"Description: {item['description'][:500]}")  # Truncate long descriptions
        if item.get("species"):
            text_parts.append(f"Species: {item['species']}")
        if item.get("assay_type"):
            text_parts.append(f"Assay: {item['assay_type']}")
        
        text = " ".join(text_parts)
        if not text.strip():
            # Return zero vector for empty items
            return np.zeros(1536)  # Default OpenAI embedding size
        
        response = client.embeddings.create(
            model="text-embedding-ada-002",
            input=text,
        )
        
        embedding = np.array(response.data[0].embedding)
        return embedding
        
    except Exception as e:
        logger.error(f"Error getting embedding: {e}")
        # Return random vector as fallback
        return np.random.random(1536)


def _get_novelty_explanation(
    item: Dict[str, Any],
    existing_items: List[Dict[str, Any]],
    novelty_score: float,
    most_similar_item_id: Optional[str],
) -> str:
    """Get LLM explanation for novelty score."""
    try:
        cfg = get_config()
        client = get_openai_client()
        
        # Find the most similar item for context
        most_similar_item = None
        if most_similar_item_id:
            for existing_item in existing_items:
                if existing_item.get("id") == most_similar_item_id or existing_item.get("accession") == most_similar_item_id:
                    most_similar_item = existing_item
                    break
        
        prompt = f"""Explain why this research item has a novelty score of {novelty_score:.2f} (scale 0-1, where 1 is most novel).

ITEM BEING EVALUATED:
Title: {item.get('title', 'N/A')}
Description: {item.get('description', 'N/A')[:300]}...

MOST SIMILAR EXISTING ITEM:
{f"Title: {most_similar_item.get('title', 'N/A')}" if most_similar_item else "No similar items found"}
{f"Description: {most_similar_item.get('description', 'N/A')[:300]}..." if most_similar_item else ""}

Provide a brief explanation (2-3 sentences) of what makes this item novel or not novel compared to existing research.
Focus on unique aspects like novel targets, diseases, methodologies, or populations studied.
"""
        
        response = client.chat.completions.create(
            model=cfg.openai.model,
            messages=[
                {
                    "role": "system",
                    "content": "You are a research novelty assessment assistant. Explain what makes research novel or similar to existing work.",
                },
                {"role": "user", "content": prompt},
            ],
            temperature=0.3,
            max_tokens=200,
        )
        
        explanation = response.choices[0].message.content
        return explanation or "No explanation available."
        
    except Exception as e:
        logger.error(f"Error getting novelty explanation: {e}")
        if novelty_score > 0.7:
            return "High novelty - appears to be substantially different from existing items."
        elif novelty_score > 0.3:
            return "Moderate novelty - some similarities to existing items but with unique aspects."
        else:
            return "Low novelty - very similar to existing items in the dataset."


__all__ = [
    "RelevanceScore",
    "NoveltyScore",
    "ScoredItem",
    "score_relevance",
    "score_novelty",
    "batch_score",
]
