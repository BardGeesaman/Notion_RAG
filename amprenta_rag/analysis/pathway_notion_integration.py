"""
Notion integration for pathway analysis results.

Links pathway enrichment results back to Notion, creating pathway pages
and establishing relations between pathways, features, datasets, and signatures.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Set

import requests

from amprenta_rag.analysis.pathway_analysis import Pathway, PathwayEnrichmentResult
from amprenta_rag.clients.notion_client import notion_headers
from amprenta_rag.config import get_config
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def create_or_update_pathway_page(
    pathway: Pathway,
    enrichment_result: Optional[PathwayEnrichmentResult] = None,
) -> Optional[str]:
    """
    Create or update a pathway page in Notion.
    
    Args:
        pathway: Pathway object
        enrichment_result: Optional enrichment result for this pathway
        
    Returns:
        Notion page ID if created/updated, None otherwise
    """
    cfg = get_config()
    
    if not hasattr(cfg.notion, "pathways_db_id") or not cfg.notion.pathways_db_id:
        logger.warning(
            "[ANALYSIS][PATHWAY-NOTION] Pathways DB ID not configured. "
            "Set NOTION_PATHWAYS_DB_ID in .env"
        )
        return None
    
    db_id = cfg.notion.pathways_db_id
    
    # Check if pathway page already exists
    existing_page_id = _find_pathway_page(pathway.pathway_id, db_id)
    
    # Build page properties
    props: Dict[str, Any] = {
        "Pathway Name": {"title": [{"text": {"content": pathway.name}}]},
        "Pathway ID": {"rich_text": [{"text": {"content": pathway.pathway_id}}]},
        "Source": {"select": {"name": pathway.source}},
    }
    
    if pathway.description:
        props["Description"] = {"rich_text": [{"text": {"content": pathway.description[:2000]}}]}  # Notion limit
    
    # Add enrichment metrics if available
    if enrichment_result:
        if enrichment_result.p_value is not None:
            props["P-value"] = {"number": enrichment_result.p_value}
        if enrichment_result.adjusted_p_value is not None:
            props["FDR P-value"] = {"number": enrichment_result.adjusted_p_value}
        if enrichment_result.enrichment_ratio is not None:
            props["Enrichment Ratio"] = {"number": enrichment_result.enrichment_ratio}
        if enrichment_result.input_features is not None:
            props["Matched Features"] = {"number": enrichment_result.input_features}
        if enrichment_result.pathway_size is not None:
            props["Pathway Size"] = {"number": enrichment_result.pathway_size}
    
    try:
        if existing_page_id:
            # Update existing page
            url = f"{cfg.notion.base_url}/pages/{existing_page_id}"
            payload = {"properties": props}
            
            resp = requests.patch(
                url,
                headers=notion_headers(),
                json=payload,
                timeout=30,
            )
            resp.raise_for_status()
            
            logger.info(
                "[ANALYSIS][PATHWAY-NOTION] Updated pathway page %s for %s",
                existing_page_id,
                pathway.pathway_id,
            )
            return existing_page_id
        else:
            # Create new page
            url = f"{cfg.notion.base_url}/pages"
            payload = {
                "parent": {"database_id": db_id},
                "properties": props,
            }
            
            resp = requests.post(
                url,
                headers=notion_headers(),
                json=payload,
                timeout=30,
            )
            resp.raise_for_status()
            
            page_id = resp.json().get("id", "")
            logger.info(
                "[ANALYSIS][PATHWAY-NOTION] Created pathway page %s for %s",
                page_id,
                pathway.pathway_id,
            )
            return page_id
            
    except Exception as e:
        logger.error(
            "[ANALYSIS][PATHWAY-NOTION] Error creating/updating pathway page for %s: %r",
            pathway.pathway_id,
            e,
        )
        return None


def _find_pathway_page(pathway_id: str, db_id: str) -> Optional[str]:
    """Find existing pathway page by pathway ID."""
    cfg = get_config()
    
    try:
        url = f"{cfg.notion.base_url}/databases/{db_id}/query"
        payload = {
            "filter": {
                "property": "Pathway ID",
                "rich_text": {"equals": pathway_id},
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
            return results[0].get("id", "")
        
        return None
        
    except Exception as e:
        logger.debug(
            "[ANALYSIS][PATHWAY-NOTION] Error finding pathway page: %r",
            e,
        )
        return None


def link_pathway_to_features(
    pathway_page_id: str,
    feature_page_ids: List[str],
    feature_type: str,
) -> None:
    """
    Link a pathway page to feature pages in Notion.
    
    Args:
        pathway_page_id: Notion page ID of pathway
        feature_page_ids: List of feature page IDs to link
        feature_type: Type of features (gene, protein, metabolite, lipid)
    """
    cfg = get_config()
    
    if not feature_page_ids:
        return
    
    # Determine the relation property name based on feature type
    relation_property_map = {
        "gene": "Related Features",  # Or "Related Genes" if exists
        "protein": "Related Features",  # Or "Related Proteins" if exists
        "metabolite": "Related Features",  # Or "Related Metabolites" if exists
        "lipid": "Related Features",  # Or "Related Lipids" if exists
    }
    
    relation_property = relation_property_map.get(feature_type, "Related Features")
    
    try:
        # Fetch current relations
        url = f"{cfg.notion.base_url}/pages/{pathway_page_id}"
        resp = requests.get(url, headers=notion_headers(), timeout=30)
        resp.raise_for_status()
        
        page = resp.json()
        props = page.get("properties", {}) or {}
        
        # Find the relation property (may have different names)
        current_relations = []
        for prop_name, prop_data in props.items():
            if prop_data.get("type") == "relation":
                if "feature" in prop_name.lower() or "gene" in prop_name.lower() or "protein" in prop_name.lower():
                    current_relations = prop_data.get("relation", []) or []
                    relation_property = prop_name
                    break
        
        # Add new relations (avoid duplicates)
        existing_ids = {r.get("id", "") for r in current_relations if r.get("id")}
        new_relations = [{"id": fid} for fid in feature_page_ids if fid not in existing_ids]
        
        if new_relations:
            updated_relations = current_relations + new_relations
            
            # Update pathway page
            update_url = f"{cfg.notion.base_url}/pages/{pathway_page_id}"
            payload = {
                "properties": {
                    relation_property: {"relation": updated_relations},
                },
            }
            
            resp = requests.patch(
                update_url,
                headers=notion_headers(),
                json=payload,
                timeout=30,
            )
            resp.raise_for_status()
            
            logger.info(
                "[ANALYSIS][PATHWAY-NOTION] Linked %d %s features to pathway %s",
                len(new_relations),
                feature_type,
                pathway_page_id,
            )
        
    except Exception as e:
        logger.warning(
            "[ANALYSIS][PATHWAY-NOTION] Error linking features to pathway: %r",
            e,
        )


def link_pathway_to_dataset(
    pathway_page_id: str,
    dataset_page_id: str,
) -> None:
    """
    Link a pathway page to a dataset page in Notion.
    
    Args:
        pathway_page_id: Notion page ID of pathway
        dataset_page_id: Notion page ID of dataset
    """
    cfg = get_config()
    
    try:
        # Fetch current relations
        url = f"{cfg.notion.base_url}/pages/{pathway_page_id}"
        resp = requests.get(url, headers=notion_headers(), timeout=30)
        resp.raise_for_status()
        
        page = resp.json()
        props = page.get("properties", {}) or {}
        
        # Find the relation property for datasets
        relation_property = None
        current_relations = []
        
        for prop_name, prop_data in props.items():
            if prop_data.get("type") == "relation":
                prop_lower = prop_name.lower()
                if "dataset" in prop_lower or "data asset" in prop_lower:
                    current_relations = prop_data.get("relation", []) or []
                    relation_property = prop_name
                    break
        
        if not relation_property:
            logger.debug(
                "[ANALYSIS][PATHWAY-NOTION] No dataset relation property found for pathway %s",
                pathway_page_id,
            )
            return
        
        # Add dataset relation (avoid duplicates)
        existing_ids = {r.get("id", "") for r in current_relations if r.get("id")}
        if dataset_page_id not in existing_ids:
            updated_relations = current_relations + [{"id": dataset_page_id}]
            
            # Update pathway page
            update_url = f"{cfg.notion.base_url}/pages/{pathway_page_id}"
            payload = {
                "properties": {
                    relation_property: {"relation": updated_relations},
                },
            }
            
            resp = requests.patch(
                update_url,
                headers=notion_headers(),
                json=payload,
                timeout=30,
            )
            resp.raise_for_status()
            
            logger.info(
                "[ANALYSIS][PATHWAY-NOTION] Linked dataset %s to pathway %s",
                dataset_page_id,
                pathway_page_id,
            )
        
    except Exception as e:
        logger.warning(
            "[ANALYSIS][PATHWAY-NOTION] Error linking dataset to pathway: %r",
            e,
        )


def link_pathway_to_signature(
    pathway_page_id: str,
    signature_page_id: str,
) -> None:
    """
    Link a pathway page to a signature page in Notion.
    
    Args:
        pathway_page_id: Notion page ID of pathway
        signature_page_id: Notion page ID of signature
    """
    cfg = get_config()
    
    try:
        # Fetch current relations
        url = f"{cfg.notion.base_url}/pages/{pathway_page_id}"
        resp = requests.get(url, headers=notion_headers(), timeout=30)
        resp.raise_for_status()
        
        page = resp.json()
        props = page.get("properties", {}) or {}
        
        # Find the relation property for signatures
        relation_property = None
        current_relations = []
        
        for prop_name, prop_data in props.items():
            if prop_data.get("type") == "relation":
                prop_lower = prop_name.lower()
                if "signature" in prop_lower:
                    current_relations = prop_data.get("relation", []) or []
                    relation_property = prop_name
                    break
        
        if not relation_property:
            logger.debug(
                "[ANALYSIS][PATHWAY-NOTION] No signature relation property found for pathway %s",
                pathway_page_id,
            )
            return
        
        # Add signature relation (avoid duplicates)
        existing_ids = {r.get("id", "") for r in current_relations if r.get("id")}
        if signature_page_id not in existing_ids:
            updated_relations = current_relations + [{"id": signature_page_id}]
            
            # Update pathway page
            update_url = f"{cfg.notion.base_url}/pages/{pathway_page_id}"
            payload = {
                "properties": {
                    relation_property: {"relation": updated_relations},
                },
            }
            
            resp = requests.patch(
                update_url,
                headers=notion_headers(),
                json=payload,
                timeout=30,
            )
            resp.raise_for_status()
            
            logger.info(
                "[ANALYSIS][PATHWAY-NOTION] Linked signature %s to pathway %s",
                signature_page_id,
                pathway_page_id,
            )
        
    except Exception as e:
        logger.warning(
            "[ANALYSIS][PATHWAY-NOTION] Error linking signature to pathway: %r",
            e,
        )


def update_dataset_with_pathway_summary(
    dataset_page_id: str,
    enrichment_results: List[PathwayEnrichmentResult],
) -> None:
    """
    Update a dataset page with pathway enrichment summary.
    
    Args:
        dataset_page_id: Notion page ID of dataset
        enrichment_results: List of pathway enrichment results
    """
    cfg = get_config()
    
    if not enrichment_results:
        return
    
    try:
        # Generate summary text
        summary_parts = [
            f"Found {len(enrichment_results)} significantly enriched pathway(s).\n\n"
        ]
        
        # Top 5 pathways
        top_pathways = enrichment_results[:5]
        for i, result in enumerate(top_pathways, 1):
            pathway = result.pathway
            summary_parts.append(
                f"{i}. {pathway.name} ({pathway.source}) - "
                f"p={result.adjusted_p_value:.4f}, "
                f"enrichment={result.enrichment_ratio:.2f}x\n"
            )
        
        summary_text = "".join(summary_parts)
        
        # Update dataset page
        url = f"{cfg.notion.base_url}/pages/{dataset_page_id}"
        
        # Try to find a pathway summary property
        resp = requests.get(url, headers=notion_headers(), timeout=30)
        resp.raise_for_status()
        page = resp.json()
        props = page.get("properties", {}) or {}
        
        # Look for pathway-related properties
        pathway_property = None
        for prop_name in ["Pathway Summary", "Enriched Pathways", "Pathway Enrichment"]:
            if prop_name in props:
                pathway_property = prop_name
                break
        
        if pathway_property:
            payload = {
                "properties": {
                    pathway_property: {
                        "rich_text": [{"text": {"content": summary_text}}],
                    },
                },
            }
            
            resp = requests.patch(
                url,
                headers=notion_headers(),
                json=payload,
                timeout=30,
            )
            resp.raise_for_status()
            
            logger.info(
                "[ANALYSIS][PATHWAY-NOTION] Updated pathway summary for dataset %s",
                dataset_page_id,
            )
        else:
            logger.debug(
                "[ANALYSIS][PATHWAY-NOTION] No pathway summary property found for dataset %s",
                dataset_page_id,
            )
        
    except Exception as e:
        logger.warning(
            "[ANALYSIS][PATHWAY-NOTION] Error updating dataset pathway summary: %r",
            e,
        )

