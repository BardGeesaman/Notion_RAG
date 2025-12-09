"""
Helper functions for cross-omics reasoning.

Shared utilities for extracting properties, retrieving chunks,
and grouping data by omics type.

Notion support has been removed - Postgres is now the source of truth.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Set

from amprenta_rag.config import get_config
from amprenta_rag.logging_utils import get_logger
from amprenta_rag.query.pinecone_query import query_pinecone

logger = get_logger(__name__)


def get_page_text(page_id: str) -> str:
    """Stub: Notion support removed. Returns empty string."""
    logger.debug("[CROSS-OMICS][HELPERS] get_page_text() is a no-op (Notion removed)")
    return ""


def fetch_notion_page(page_id: str) -> Optional[Dict[str, Any]]:
    """Stub: Notion support removed. Returns None."""
    logger.debug("[CROSS-OMICS][HELPERS] fetch_notion_page() is a no-op (Notion removed)")
    return None


def extract_relation_ids(page: Dict[str, Any], property_name: str) -> List[str]:
    """Extract page IDs from a relation property."""
    props = page.get("properties", {}) or {}
    relation_prop = props.get(property_name, {})
    if relation_prop.get("type") == "relation":
        relations = relation_prop.get("relation", []) or []
        return [r.get("id", "") for r in relations if r.get("id")]
    return []


def extract_select_values(page: Dict[str, Any], property_name: str) -> List[str]:
    """Extract values from a select or multi_select property."""
    props = page.get("properties", {}) or {}
    select_prop = props.get(property_name, {})
    
    if select_prop.get("type") == "select":
        select_val = select_prop.get("select")
        if select_val:
            return [select_val.get("name", "")]
    
    elif select_prop.get("type") == "multi_select":
        multi_select = select_prop.get("multi_select", []) or []
        return [item.get("name", "") for item in multi_select if item.get("name")]
    
    return []


def extract_text_property(page: Dict[str, Any], property_name: str) -> Optional[str]:
    """Extract text from a title or rich_text property."""
    props = page.get("properties", {}) or {}
    text_prop = props.get(property_name, {})
    
    if text_prop.get("type") == "title":
        title_parts = text_prop.get("title", []) or []
        if title_parts:
            return title_parts[0].get("plain_text", "")
    
    elif text_prop.get("type") == "rich_text":
        rich_text = text_prop.get("rich_text", []) or []
        if rich_text:
            return "".join(rt.get("plain_text", "") for rt in rich_text)
    
    return None


def get_chunk_text(chunk: Dict[str, Any]) -> Optional[str]:
    """Get full chunk text from Notion, fallback to snippet."""
    meta = chunk.get("metadata", {}) or getattr(chunk, "metadata", {})
    
    # Try to get full chunk text from Notion
    chunk_page_id = meta.get("notion_chunk_page_id")
    if chunk_page_id:
        try:
            full_text = get_page_text(chunk_page_id)
            if full_text:
                return full_text
        except Exception:
            pass
    
    # Fallback to snippet
    snippet = meta.get("snippet", "")
    if snippet:
        return snippet
    
    return None


def retrieve_chunks_for_objects(
    object_ids: List[str],
    object_type: str,
    top_k_per_object: int = 20,
) -> List[Dict[str, Any]]:
    """
    Retrieve Pinecone chunks for a list of Notion objects.
    
    Args:
        object_ids: List of Notion page IDs
        object_type: Type of objects ("dataset", "experiment", "signature", etc.)
        top_k_per_object: Maximum chunks per object
        
    Returns:
        List of match dictionaries from Pinecone
    """
    if not object_ids:
        return []
    
    all_matches: List[Dict[str, Any]] = []
    
    # Query Pinecone for chunks associated with these objects
    # Use metadata filters to find chunks with matching page IDs
    for obj_id in object_ids:
        try:
            # Remove dashes from ID for metadata matching (Pinecone stores IDs without dashes)
            obj_id_clean = obj_id.replace("-", "")
            
            # Build filter based on object type
            meta_filter: Dict[str, Any] = {}
            
            if object_type == "dataset":
                meta_filter["dataset_page_id"] = obj_id_clean
            elif object_type == "experiment":
                meta_filter["experiment_page_id"] = obj_id_clean
            elif object_type == "signature":
                meta_filter["signature_page_id"] = obj_id_clean
            elif object_type == "program":
                meta_filter["program_page_id"] = obj_id_clean
            
            # Query Pinecone with a generic query to get all chunks for this object
            query_text = f"{object_type} data analysis results"
            matches = query_pinecone(
                user_query=query_text,
                top_k=top_k_per_object,
                meta_filter=meta_filter,
                source_types=None,
            )
            
            all_matches.extend(matches)
            
        except Exception as e:
            logger.warning(
                "[RAG][CROSS-OMICS] Error retrieving chunks for %s %s: %r",
                object_type,
                obj_id,
                e,
            )
            continue
    
    # Deduplicate by chunk ID
    seen_ids: Set[str] = set()
    unique_matches: List[Dict[str, Any]] = []
    for match in all_matches:
        match_id = match.get("id") or getattr(match, "id", None)
        if match_id and match_id not in seen_ids:
            seen_ids.add(match_id)
            unique_matches.append(match)
    
    logger.info(
        "[RAG][CROSS-OMICS] Retrieved %d unique chunks for %d %s objects",
        len(unique_matches),
        len(object_ids),
        object_type,
    )
    
    return unique_matches


def group_chunks_by_omics_type(chunks: List[Dict[str, Any]]) -> Dict[str, List[Dict[str, Any]]]:
    """Group chunks by omics type."""
    grouped: Dict[str, List[Dict[str, Any]]] = {
        "Lipidomics": [],
        "Metabolomics": [],
        "Proteomics": [],
        "Transcriptomics": [],
        "Other": [],
    }
    
    for chunk in chunks:
        meta = chunk.get("metadata", {}) or getattr(chunk, "metadata", {})
        omics_type = meta.get("omics_type", "Other")
        
        if omics_type in grouped:
            grouped[omics_type].append(chunk)
        else:
            grouped["Other"].append(chunk)
    
    return grouped

