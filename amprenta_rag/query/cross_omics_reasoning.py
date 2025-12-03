"""
Cross-omics RAG reasoning capabilities.

High-level reasoning functions that generate multi-omics summaries across:
- Programs
- Signatures
- Features (genes, proteins, metabolites, lipids)
- Datasets
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Set

import requests

from amprenta_rag.clients.notion_client import notion_headers
from amprenta_rag.clients.notion_client import get_page_text
from amprenta_rag.clients.openai_client import get_default_models, get_openai_client
from amprenta_rag.config import get_config
from amprenta_rag.logging_utils import get_logger
from amprenta_rag.query.pinecone_query import build_meta_filter, query_pinecone

logger = get_logger(__name__)

# Maximum context chunks to send to LLM
MAX_CONTEXT_CHUNKS = 50
MAX_CHUNK_LENGTH = 2000


def _fetch_notion_page(page_id: str) -> Optional[Dict[str, Any]]:
    """Helper to fetch a Notion page."""
    cfg = get_config()
    try:
        url = f"{cfg.notion.base_url}/pages/{page_id}"
        resp = requests.get(url, headers=notion_headers(), timeout=30)
        resp.raise_for_status()
        return resp.json()
    except Exception as e:
        logger.warning(
            "[RAG][CROSS-OMICS] Error fetching Notion page %s: %r",
            page_id,
            e,
        )
        return None


def _extract_relation_ids(page: Dict[str, Any], property_name: str) -> List[str]:
    """Extract page IDs from a relation property."""
    props = page.get("properties", {}) or {}
    relation_prop = props.get(property_name, {})
    if relation_prop.get("type") == "relation":
        relations = relation_prop.get("relation", []) or []
        return [r.get("id", "") for r in relations if r.get("id")]
    return []


def _extract_select_values(page: Dict[str, Any], property_name: str) -> List[str]:
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


def _extract_text_property(page: Dict[str, Any], property_name: str) -> Optional[str]:
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


def _get_chunk_text(chunk: Dict[str, Any]) -> Optional[str]:
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


def _retrieve_chunks_for_objects(
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


def _group_chunks_by_omics_type(chunks: List[Dict[str, Any]]) -> Dict[str, List[Dict[str, Any]]]:
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


def _synthesize_cross_omics_summary(
    prompt: str,
    context_chunks: List[str],
    max_chunks: int = MAX_CONTEXT_CHUNKS,
) -> str:
    """Use OpenAI to synthesize a cross-omics summary."""
    client = get_openai_client()
    chat_model, _ = get_default_models()
    
    # Truncate chunks to fit token budget
    truncated_chunks = context_chunks[:max_chunks]
    
    # Limit each chunk length
    limited_chunks = [
        chunk[:MAX_CHUNK_LENGTH] + ("..." if len(chunk) > MAX_CHUNK_LENGTH else "")
        for chunk in truncated_chunks
    ]
    
    context = "\n\n---\n\n".join(limited_chunks)
    
    system_prompt = (
        "You are an expert in multi-omics analysis (lipidomics, metabolomics, proteomics, transcriptomics). "
        "You are given snippets of information about multi-omics evidence.\n\n"
        "Summarize the cross-omics evidence with the following structure:\n\n"
        "1. High-level context\n"
        "2. Per-omics findings:\n"
        "   - Lipidomics\n"
        "   - Metabolomics\n"
        "   - Proteomics\n"
        "   - Transcriptomics\n"
        "3. Cross-omics convergence:\n"
        "   - Features (genes/proteins/metabolites/lipids) that consistently change across multiple omics.\n"
        "4. Cross-omics divergence:\n"
        "   - Conflicting signals or modality-specific changes.\n"
        "5. Disease, model system, and matrix context.\n"
        "6. Key open questions and next experimental steps.\n\n"
        "Only use information provided in the context. Do not hallucinate external facts. "
        "Label modality-specific findings clearly. Be concise but comprehensive."
    )
    
    user_content = f"{prompt}\n\nRelevant context chunks:\n{context}"
    
    try:
        resp = client.chat.completions.create(
            model=chat_model,
            messages=[
                {"role": "system", "content": system_prompt},
                {"role": "user", "content": user_content},
            ],
            temperature=0.3,
        )
        return resp.choices[0].message.content.strip()  # type: ignore[union-attr]
    except Exception as e:
        logger.error(
            "[RAG][CROSS-OMICS] OpenAI API error synthesizing summary: %r",
            e,
        )
        raise


def cross_omics_program_summary(
    program_page_id: str,
    top_k_per_omics: int = 20,
) -> str:
    """
    Generate a cross-omics summary for a Program.

    Args:
        program_page_id: Notion page ID of Program (with dashes)
        top_k_per_omics: Maximum chunks to retrieve per omics type

    Returns:
        Textual summary (markdown format)
    """
    logger.info(
        "[RAG][CROSS-OMICS] Generating cross-omics summary for program %s",
        program_page_id,
    )
    
    # Fetch program page
    program_page = _fetch_notion_page(program_page_id)
    if not program_page:
        return f"Error: Could not fetch Program page {program_page_id}."
    
    # Try "Program" property first (actual schema), fallback to "Name"
    program_name = _extract_text_property(program_page, "Program") or _extract_text_property(program_page, "Name") or "Unknown Program"
    
    # Find related Experiments (schema shows "Experiments" relation)
    # Try "Experiments" first, then fallback to "Related Experiments"
    experiment_ids = _extract_relation_ids(program_page, "Experiments")
    if not experiment_ids:
        experiment_ids = _extract_relation_ids(program_page, "Related Experiments")
    
    # Datasets might be linked via Experiments, or check for "Program Datasets"
    dataset_ids = _extract_relation_ids(program_page, "Program Datasets")
    
    logger.info(
        "[RAG][CROSS-OMICS] Found %d experiments, %d datasets for program %s",
        len(experiment_ids),
        len(dataset_ids),
        program_name,
    )
    
    if not experiment_ids and not dataset_ids:
        return (
            f"No sufficient multi-omics context found for program '{program_name}'. "
            "No experiments or datasets are linked to this program."
        )
    
    # Retrieve chunks from datasets
    dataset_chunks = _retrieve_chunks_for_objects(
        dataset_ids,
        "dataset",
        top_k_per_object=top_k_per_omics,
    )
    
    # Retrieve chunks from experiments
    experiment_chunks = _retrieve_chunks_for_objects(
        experiment_ids,
        "experiment",
        top_k_per_object=top_k_per_omics,
    )
    
    all_chunks = dataset_chunks + experiment_chunks
    
    if not all_chunks:
        return (
            f"No chunks found for program '{program_name}'. "
            "Linked datasets and experiments may not have been ingested into RAG."
        )
    
    # Group chunks by omics type
    chunks_by_omics = _group_chunks_by_omics_type(all_chunks)
    
    # Get chunk texts - try to get full text from Notion, fallback to snippet
    context_chunks: List[str] = []
    for chunk in all_chunks:
        meta = chunk.get("metadata", {}) or getattr(chunk, "metadata", {})
        title = meta.get("title", "")
        source = meta.get("source_type", "")
        
        # Try to get full chunk text from Notion
        chunk_text = ""
        chunk_page_id = meta.get("notion_chunk_page_id")
        if chunk_page_id:
            try:
                full_text = get_page_text(chunk_page_id)
                if full_text:
                    chunk_text = full_text
            except Exception as e:
                logger.debug(
                    "[RAG][CROSS-OMICS] Could not fetch full text for chunk %s: %r",
                    chunk_page_id,
                    e,
                )
        
        # Fallback to snippet if no full text
        if not chunk_text:
            snippet = meta.get("snippet", "")
            chunk_text = snippet
        
        if chunk_text:
            context_chunks.append(f"[{source}] {title}\n{chunk_text}")
    
    omics_counts = {
        omics: len(chunks)
        for omics, chunks in chunks_by_omics.items()
        if chunks
    }
    
    logger.info(
        "[RAG][CROSS-OMICS] Retrieved %d total chunks for program %s: %s",
        len(all_chunks),
        program_name,
        omics_counts,
    )
    
    # Build prompt
    prompt = f"""Generate a cross-omics summary for the program: {program_name}

The program is linked to:
- {len(experiment_ids)} experiment(s)
- {len(dataset_ids)} dataset(s)

Chunks retrieved:
- {omics_counts.get('Lipidomics', 0)} lipidomics chunks
- {omics_counts.get('Metabolomics', 0)} metabolomics chunks
- {omics_counts.get('Proteomics', 0)} proteomics chunks
- {omics_counts.get('Transcriptomics', 0)} transcriptomics chunks
"""
    
    # Synthesize summary
    summary = _synthesize_cross_omics_summary(prompt, context_chunks)
    
    return summary


def cross_omics_signature_summary(
    signature_page_id: str,
    top_k_datasets: int = 20,
    top_k_chunks: int = 100,
) -> str:
    """
    Generate a cross-omics summary for a Signature.

    Args:
        signature_page_id: Notion page ID of Signature (with dashes)
        top_k_datasets: Maximum datasets to include
        top_k_chunks: Maximum chunks to retrieve

    Returns:
        Textual summary (markdown format)
    """
    logger.info(
        "[RAG][CROSS-OMICS] Generating cross-omics summary for signature %s",
        signature_page_id,
    )
    
    # Fetch signature page
    signature_page = _fetch_notion_page(signature_page_id)
    if not signature_page:
        return f"Error: Could not fetch Signature page {signature_page_id}."
    
    signature_name = _extract_text_property(signature_page, "Name") or "Unknown Signature"
    modalities = _extract_select_values(signature_page, "Modalities")
    disease = _extract_select_values(signature_page, "Disease")
    matrix = _extract_select_values(signature_page, "Matrix")
    
    # Find datasets with this signature
    # Query Experimental Data Assets database for datasets with this signature in Related Signature(s)
    cfg = get_config()
    dataset_ids: List[str] = []
    
    try:
        # Query for datasets with this signature
        # This requires querying the Experimental Data Assets database
        # For now, we'll query Pinecone chunks to find datasets
        # A more robust approach would query Notion directly
        
        # Query Pinecone for chunks with this signature
        signature_id_clean = signature_page_id.replace("-", "")
        meta_filter = {"signature_page_id": signature_id_clean}
        
        signature_chunks = query_pinecone(
            user_query=f"signature {signature_name} multi-omics analysis",
            top_k=top_k_chunks,
            meta_filter=meta_filter,
            source_types=["Signature", "Dataset"],
        )
        
        # Extract unique dataset IDs from chunks
        for chunk in signature_chunks:
            meta = chunk.get("metadata", {}) or getattr(chunk, "metadata", {})
            dataset_id = meta.get("dataset_page_id")
            if dataset_id:
                # Add dashes back
                dataset_id_with_dashes = f"{dataset_id[:8]}-{dataset_id[8:12]}-{dataset_id[12:16]}-{dataset_id[16:20]}-{dataset_id[20:]}"
                if dataset_id_with_dashes not in dataset_ids:
                    dataset_ids.append(dataset_id_with_dashes)
        
        # Also query Experimental Data Assets database directly
        try:
            exp_data_db_id = cfg.notion.experimental_data_assets_db_id
            if exp_data_db_id:
                url = f"{cfg.notion.base_url}/databases/{exp_data_db_id}/query"
                payload = {
                    "filter": {
                        "property": "Related Signature(s)",
                        "relation": {"contains": signature_page_id},
                    },
                    "page_size": 100,
                }
                
                resp = requests.post(
                    url,
                    headers=notion_headers(),
                    json=payload,
                    timeout=30,
                )
                resp.raise_for_status()
                
                data = resp.json()
                db_dataset_ids = [r.get("id", "") for r in data.get("results", []) if r.get("id")]
                dataset_ids.extend(db_dataset_ids)
                dataset_ids = list(set(dataset_ids))  # Deduplicate
                
        except Exception as e:
            logger.debug(
                "[RAG][CROSS-OMICS] Could not query Experimental Data Assets DB: %r",
                e,
            )
    
    except Exception as e:
        logger.warning(
            "[RAG][CROSS-OMICS] Error finding datasets for signature: %r",
            e,
        )
    
    logger.info(
        "[RAG][CROSS-OMICS] Found %d datasets for signature %s",
        len(dataset_ids),
        signature_name,
    )
    
    # Limit to top_k_datasets
    dataset_ids = dataset_ids[:top_k_datasets]
    
    # Retrieve chunks from datasets and signature
    all_chunks: List[Dict[str, Any]] = []
    
    # Get signature chunks
    signature_id_clean = signature_page_id.replace("-", "")
    signature_chunks = query_pinecone(
        user_query=f"signature {signature_name}",
        top_k=20,
        meta_filter={"signature_page_id": signature_id_clean},
        source_types=["Signature"],
    )
    all_chunks.extend(signature_chunks)
    
    # Get dataset chunks
    if dataset_ids:
        dataset_chunks = _retrieve_chunks_for_objects(
            dataset_ids,
            "dataset",
            top_k_per_object=top_k_chunks // max(len(dataset_ids), 1),
        )
        all_chunks.extend(dataset_chunks)
    
    if not all_chunks:
        return (
            f"No sufficient multi-omics context found for signature '{signature_name}'. "
            "No chunks were found for this signature or its matched datasets."
        )
    
    # Group by omics type
    chunks_by_omics = _group_chunks_by_omics_type(all_chunks)
    
    # Get chunk texts
    context_chunks: List[str] = []
    for chunk in all_chunks[:top_k_chunks]:
        meta = chunk.get("metadata", {}) or getattr(chunk, "metadata", {})
        title = meta.get("title", "")
        source = meta.get("source_type", "")
        
        chunk_text = _get_chunk_text(chunk)
        if chunk_text:
            context_chunks.append(f"[{source}] {title}\n{chunk_text}")
    
    omics_counts = {
        omics: len(chunks)
        for omics, chunks in chunks_by_omics.items()
        if chunks
    }
    
    # Build prompt
    prompt = f"""Generate a cross-omics summary for the signature: {signature_name}

Signature details:
- Modalities: {', '.join(modalities) if modalities else 'Not specified'}
- Disease: {', '.join(disease) if disease else 'Not specified'}
- Matrix: {', '.join(matrix) if matrix else 'Not specified'}
- Matched datasets: {len(dataset_ids)}

Chunks retrieved:
- {omics_counts.get('Lipidomics', 0)} lipidomics chunks
- {omics_counts.get('Metabolomics', 0)} metabolomics chunks
- {omics_counts.get('Proteomics', 0)} proteomics chunks
- {omics_counts.get('Transcriptomics', 0)} transcriptomics chunks
- {len(signature_chunks)} signature chunks
"""
    
    # Synthesize summary
    summary = _synthesize_cross_omics_summary(prompt, context_chunks)
    
    return summary


def cross_omics_feature_summary(
    feature_name: str,
    feature_type: str,
    top_k_datasets: int = 20,
    top_k_chunks: int = 100,
) -> str:
    """
    Summarize all multi-omics evidence for a single feature.

    Args:
        feature_name: Feature name (e.g., "TP53", "Cer(d18:1/16:0)")
        feature_type: "gene", "protein", "metabolite", or "lipid"
        top_k_datasets: Maximum datasets to include
        top_k_chunks: Maximum chunks to retrieve

    Returns:
        Textual summary (markdown format)
    """
    logger.info(
        "[RAG][CROSS-OMICS] Generating cross-omics summary for %s feature: %s",
        feature_type,
        feature_name,
    )
    
    # Find feature page in appropriate database
    cfg = get_config()
    
    # Map feature type to database ID
    db_map = {
        "gene": cfg.notion.gene_features_db_id if hasattr(cfg.notion, "gene_features_db_id") else None,
        "protein": cfg.notion.protein_features_db_id if hasattr(cfg.notion, "protein_features_db_id") else None,
        "metabolite": cfg.notion.metabolite_features_db_id if hasattr(cfg.notion, "metabolite_features_db_id") else None,
        "lipid": cfg.notion.lipid_species_db_id if hasattr(cfg.notion, "lipid_species_db_id") else None,
    }
    
    db_id = db_map.get(feature_type)
    if not db_id:
        return f"Error: {feature_type} Features database ID not configured."
    
    # Query for feature page
    feature_page_id: Optional[str] = None
    try:
        url = f"{cfg.notion.base_url}/databases/{db_id}/query"
        payload = {
            "filter": {
                "property": "Name",
                "title": {"equals": feature_name},
            },
            "page_size": 1,
        }
        
        resp = requests.post(
            url,
            headers=notion_headers(),
            json=payload,
            timeout=30,
        )
        resp.raise_for_status()
        
        results = resp.json().get("results", [])
        if results:
            feature_page_id = results[0].get("id", "")
    except Exception as e:
        logger.warning(
            "[RAG][CROSS-OMICS] Error finding %s feature page for '%s': %r",
            feature_type,
            feature_name,
            e,
        )
    
    # Find linked datasets, experiments, programs
    dataset_ids: List[str] = []
    experiment_ids: List[str] = []
    program_ids: List[str] = []
    
    if feature_page_id:
        feature_page = _fetch_notion_page(feature_page_id)
        if feature_page:
            # Extract relations based on feature type
            relation_property_candidates = {
                "gene": ["Transcriptomics Datasets", "Datasets", "Related Datasets"],
                "protein": ["Proteomics Datasets", "Datasets", "Related Datasets"],
                "metabolite": ["Metabolomics Datasets", "Datasets", "Related Datasets"],
                "lipid": ["Experimental Data Assets", "Datasets", "Related Datasets"],
            }
            
            candidates = relation_property_candidates.get(feature_type, ["Datasets", "Related Datasets"])
            for prop_name in candidates:
                dataset_ids = _extract_relation_ids(feature_page, prop_name)
                if dataset_ids:
                    break
    
    # Query Pinecone for chunks mentioning this feature
    query_text = f"{feature_type} {feature_name} multi-omics"
    feature_chunks = query_pinecone(
        user_query=query_text,
        top_k=top_k_chunks,
        meta_filter=None,
        source_types=None,
    )
    
    # Filter chunks that actually mention the feature
    filtered_chunks = []
    feature_name_lower = feature_name.lower()
    for chunk in feature_chunks:
        meta = chunk.get("metadata", {}) or getattr(chunk, "metadata", {})
        snippet = meta.get("snippet", "").lower()
        title = meta.get("title", "").lower()
        text = f"{title} {snippet}"
        
        if feature_name_lower in text:
            filtered_chunks.append(chunk)
    
    # Also get chunks from linked datasets
    if dataset_ids:
        dataset_chunks = _retrieve_chunks_for_objects(
            dataset_ids[:top_k_datasets],
            "dataset",
            top_k_per_object=top_k_chunks // max(len(dataset_ids[:top_k_datasets]), 1),
        )
        filtered_chunks.extend(dataset_chunks)
    
    if not filtered_chunks:
        return (
            f"No sufficient multi-omics context found for {feature_type} feature '{feature_name}'. "
            "No relevant chunks or linked datasets were found."
        )
    
    # Group by omics type
    chunks_by_omics = _group_chunks_by_omics_type(filtered_chunks)
    
    # Get chunk texts - try to get full text from Notion, fallback to snippet
    context_chunks: List[str] = []
    for chunk in filtered_chunks[:top_k_chunks]:
        meta = chunk.get("metadata", {}) or getattr(chunk, "metadata", {})
        title = meta.get("title", "")
        source = meta.get("source_type", "")
        
        # Try to get full chunk text from Notion
        chunk_text = ""
        chunk_page_id = meta.get("notion_chunk_page_id")
        if chunk_page_id:
            try:
                full_text = get_page_text(chunk_page_id)
                if full_text:
                    chunk_text = full_text
            except Exception:
                pass
        
        # Fallback to snippet if no full text
        if not chunk_text:
            snippet = meta.get("snippet", "")
            chunk_text = snippet
        
        if chunk_text:
            context_chunks.append(f"[{source}] {title}\n{chunk_text}")
    
    omics_counts = {
        omics: len(chunks)
        for omics, chunks in chunks_by_omics.items()
        if chunks
    }
    
    logger.info(
        "[RAG][CROSS-OMICS] Retrieved %d chunks for %s feature '%s': %s",
        len(filtered_chunks),
        feature_type,
        feature_name,
        omics_counts,
    )
    
    # Build prompt
    prompt = f"""Generate a cross-omics summary for the {feature_type} feature: {feature_name}

This feature appears in:
- {len(dataset_ids)} linked dataset(s)

Chunks retrieved across omics:
- {omics_counts.get('Lipidomics', 0)} lipidomics chunks
- {omics_counts.get('Metabolomics', 0)} metabolomics chunks
- {omics_counts.get('Proteomics', 0)} proteomics chunks
- {omics_counts.get('Transcriptomics', 0)} transcriptomics chunks

Summarize how this feature behaves across different omics modalities and contexts.
"""
    
    # Synthesize summary
    summary = _synthesize_cross_omics_summary(prompt, context_chunks)
    
    return summary


def cross_omics_dataset_summary(
    dataset_page_id: str,
    top_k_chunks: int = 100,
) -> str:
    """
    Summarize cross-omics context for a single dataset.

    Args:
        dataset_page_id: Notion page ID of dataset (with dashes)
        top_k_chunks: Maximum chunks to retrieve

    Returns:
        Textual summary (markdown format)
    """
    logger.info(
        "[RAG][CROSS-OMICS] Generating cross-omics summary for dataset %s",
        dataset_page_id,
    )
    
    # Fetch dataset page
    dataset_page = _fetch_notion_page(dataset_page_id)
    if not dataset_page:
        return f"Error: Could not fetch Dataset page {dataset_page_id}."
    
    dataset_name = _extract_text_property(dataset_page, "Name") or "Unknown Dataset"
    omics_type = _extract_select_values(dataset_page, "Omics Type")
    signature_ids = _extract_relation_ids(dataset_page, "Related Signature(s)")
    experiment_ids = _extract_relation_ids(dataset_page, "Related Experiments")
    program_ids = _extract_relation_ids(dataset_page, "Related Programs")
    
    # Get dataset chunks
    dataset_chunks = _retrieve_chunks_for_objects(
        [dataset_page_id],
        "dataset",
        top_k_per_object=top_k_chunks,
    )
    
    # Get signature chunks
    signature_chunks: List[Dict[str, Any]] = []
    for sig_id in signature_ids[:10]:  # Limit to 10 signatures
        sig_id_clean = sig_id.replace("-", "")
        chunks = query_pinecone(
            user_query=f"signature",
            top_k=5,
            meta_filter={"signature_page_id": sig_id_clean},
            source_types=["Signature"],
        )
        signature_chunks.extend(chunks)
    
    # Get experiment chunks
    experiment_chunks = _retrieve_chunks_for_objects(
        experiment_ids,
        "experiment",
        top_k_per_object=10,
    )
    
    all_chunks = dataset_chunks + signature_chunks + experiment_chunks
    
    if not all_chunks:
        return (
            f"No sufficient multi-omics context found for dataset '{dataset_name}'. "
            "No chunks were found for this dataset or its related signatures/experiments."
        )
    
    # Group by omics type
    chunks_by_omics = _group_chunks_by_omics_type(all_chunks)
    
    # Get chunk texts
    context_chunks: List[str] = []
    for chunk in all_chunks[:top_k_chunks]:
        meta = chunk.get("metadata", {}) or getattr(chunk, "metadata", {})
        title = meta.get("title", "")
        source = meta.get("source_type", "")
        
        chunk_text = _get_chunk_text(chunk)
        if chunk_text:
            context_chunks.append(f"[{source}] {title}\n{chunk_text}")
    
    omics_counts = {
        omics: len(chunks)
        for omics, chunks in chunks_by_omics.items()
        if chunks
    }
    
    # Build prompt
    prompt = f"""Generate a cross-omics summary for the dataset: {dataset_name}

Dataset details:
- Omics Type: {', '.join(omics_type) if omics_type else 'Not specified'}
- Linked signatures: {len(signature_ids)}
- Linked experiments: {len(experiment_ids)}
- Linked programs: {len(program_ids)}

Context chunks retrieved:
- {omics_counts.get('Lipidomics', 0)} lipidomics chunks
- {omics_counts.get('Metabolomics', 0)} metabolomics chunks
- {omics_counts.get('Proteomics', 0)} proteomics chunks
- {omics_counts.get('Transcriptomics', 0)} transcriptomics chunks
- {len(signature_chunks)} signature chunks
- {len(experiment_chunks)} experiment chunks
"""
    
    # Synthesize summary
    summary = _synthesize_cross_omics_summary(prompt, context_chunks)
    
    return summary

