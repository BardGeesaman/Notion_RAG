# amprenta_rag/ingestion/metadata_semantic.py

from __future__ import annotations

from typing import Dict, Any, List, Optional

import json
import requests
import re

from amprenta_rag.config import get_config
from amprenta_rag.clients.notion_client import notion_headers
from amprenta_rag.ingestion.zotero_api import ZoteroItem
from amprenta_rag.logging_utils import get_logger

logger = get_logger(__name__)


def _fetch_notion_page(page_id: str) -> Dict[str, Any]:
    cfg = get_config().notion
    url = f"{cfg.base_url}/pages/{page_id}"
    resp = requests.get(url, headers=notion_headers())
    resp.raise_for_status()
    return resp.json()


def _get_select_name(props: Dict[str, Any], name: str) -> Optional[str]:
    sel = props.get(name, {}).get("select")
    if sel and sel.get("name"):
        return sel["name"]
    return None


def _get_multi_names(props: Dict[str, Any], name: str) -> List[str]:
    ms = props.get(name, {}).get("multi_select", []) or []
    return [x["name"] for x in ms if x.get("name")]


def _get_relation_ids(props: Dict[str, Any], name: str) -> List[str]:
    """Extract page IDs from a Notion relation property."""
    rel = props.get(name, {}).get("relation", []) or []
    return [r.get("id", "").replace("-", "") for r in rel if r.get("id")]


def enforce_signature_reverse_link(signature_page_id: str, literature_page_id: str) -> None:
    """
    Create reverse link from signature page to literature/email page.
    Updates the signature page's "Source Papers" relation to include the literature page.
    """
    cfg = get_config().notion
    try:
        # Fetch current signature page
        sig_page = _fetch_notion_page(signature_page_id)
        props = sig_page.get("properties", {}) or {}
        current_rel = props.get("Source Papers", {}).get("relation", []) or []
        existing_ids = {r.get("id", "").replace("-", "") for r in current_rel if r.get("id")}
        
        # Check if literature page is already linked
        lit_page_id_clean = literature_page_id.replace("-", "")
        if lit_page_id_clean in existing_ids:
            return  # Already linked, skip
        
        # Add literature page to relation
        updated_rel = [{"id": lit_page_id_clean}] + current_rel
        
        # Update signature page
        payload = {
            "properties": {
                "Source Papers": {"relation": updated_rel},
            },
        }
        
        url = f"{cfg.base_url}/pages/{signature_page_id}"
        resp = requests.patch(
            url,
            headers=notion_headers(),
            data=json.dumps(payload),
        )
        
        if resp.status_code >= 300:
            logger.error(
                "[METADATA] Failed to update signature page %s reverse link: %s",
                signature_page_id,
                resp.text,
            )
        else:
            logger.info(
                "[METADATA] Added reverse link: signature %s -> literature %s",
                signature_page_id,
                literature_page_id,
            )
            
    except Exception as e:
        logger.error(
            "[METADATA] Error creating reverse link for signature %s -> literature %s: %r",
            signature_page_id,
            literature_page_id,
            e,
        )


def _collect_signature_metadata(signature_ids: List[str]) -> Dict[str, Any]:
    """
    Fetch Lipid Signature pages and collect their metadata.
    
    For each signature page ID, fetches the page and extracts:
    - Short ID (rich_text or title)
    - Biomarker Role (multi_select)
    - Phenotype Axes (multi_select)
    - Data Ownership (select)
    
    Args:
        signature_ids: List of Notion page IDs (without dashes) for signature pages
        
    Returns:
        Dictionary with:
        - sig_short_ids: List[str] - Short IDs from signature pages
        - sig_roles: List[str] - Biomarker roles (aggregated)
        - sig_axes: List[str] - Phenotype axes (aggregated)
        - sig_ownership: List[str] - Data ownership values (aggregated)
    """
    sig_short_ids: List[str] = []
    sig_roles: List[str] = []
    sig_axes: List[str] = []
    sig_ownership: List[str] = []
    
    if not signature_ids:
        return {
            "sig_short_ids": [],
            "sig_roles": [],
            "sig_axes": [],
            "sig_ownership": [],
        }
    
    for sig_id in signature_ids:
        try:
            sig_page = _fetch_notion_page(sig_id)
            sig_props = sig_page.get("properties", {}) or {}
            
            # Extract Short ID from rich_text or title, fallback to Name property
            short_id = None
            short_id_prop = sig_props.get("Short ID", {})
            if "rich_text" in short_id_prop:
                rich_text = short_id_prop.get("rich_text", []) or []
                if rich_text:
                    short_id = rich_text[0].get("plain_text", "").strip()
            elif "title" in short_id_prop:
                title = short_id_prop.get("title", []) or []
                if title:
                    short_id = title[0].get("plain_text", "").strip()
            
            # Fallback to Name property if Short ID is missing/empty
            if not short_id:
                name_prop = sig_props.get("Name", {})
                if "title" in name_prop:
                    title = name_prop.get("title", []) or []
                    if title:
                        short_id = title[0].get("plain_text", "").strip()
                elif "rich_text" in name_prop:
                    rich_text = name_prop.get("rich_text", []) or []
                    if rich_text:
                        short_id = rich_text[0].get("plain_text", "").strip()
            
            if short_id:
                sig_short_ids.append(short_id)
            
            # Extract Biomarker Role (multi-select)
            roles = _get_multi_names(sig_props, "Biomarker Role")
            sig_roles.extend(roles)
            
            # Extract Phenotype Axes (multi-select)
            axes = _get_multi_names(sig_props, "Phenotype Axes")
            sig_axes.extend(axes)
            
            # Extract Data Ownership (select)
            ownership = _get_select_name(sig_props, "Data Ownership")
            if ownership:
                sig_ownership.append(ownership)
                
        except Exception as e:
            logger.warning(
                "[METADATA] Error fetching signature page %s: %r",
                sig_id,
                e,
            )
            continue
    
    return {
        "sig_short_ids": sig_short_ids,
        "sig_roles": list(set(sig_roles)),  # Deduplicate
        "sig_axes": list(set(sig_axes)),  # Deduplicate
        "sig_ownership": list(set(sig_ownership)),  # Deduplicate
    }


def get_literature_semantic_metadata(parent_page_id: str, item: ZoteroItem) -> Dict[str, Any]:
    """
    Read semantic + lipid metadata from the Literature DB page and build
    a doc-level + lipid-level metadata dict that we will mirror into Pinecone.

    Expects (but does not require) these Notion properties on the Literature DB:
      - Disease (multi-select)
      - Targets (multi-select)
      - Modality (multi-select)
      - Stage (select)
      - Model Systems (multi-select)
      - Biomarker Role (multi-select)
      - Importance (number)
      - Lipid Species (raw) (multi-select or text)
      - Canonical Lipid Species (relation) – optional
      - Lipid Signatures (relation) – links to Lipid Signatures DB pages
      - Phenotype Axes (multi-select)
      - Matrix (multi-select)
      - Treatment Arms (multi-select)
      
    For Lipid Signatures relation: Fetches linked signature pages and extracts:
      - Short ID (rich_text or title)
      - Biomarker Role (multi-select)
      - Phenotype Axes (multi-select)
      - Data Ownership (select)
    """
    page = _fetch_notion_page(parent_page_id)
    props = page.get("properties", {}) or {}

    # --- Basic semantic fields ---
    diseases = _get_multi_names(props, "Disease")
    targets = _get_multi_names(props, "Targets")
    modality = _get_multi_names(props, "Modality")
    stage = _get_select_name(props, "Stage")
    model_systems = _get_multi_names(props, "Model Systems")
    biomarker_role = _get_multi_names(props, "Biomarker Role")

    importance = props.get("Importance", {}).get("number")

    phenotype_axes = _get_multi_names(props, "Phenotype Axes")
    matrix = _get_multi_names(props, "Matrix")
    treatment_arms = _get_multi_names(props, "Treatment Arms")

    # --- Lipid-related fields ---
    lipids_raw = _get_multi_names(props, "Lipid Species (raw)")

    canonical_rel = props.get("Canonical Lipid Species", {}).get("relation", []) or []
    canonical_lipids = [r.get("id", "").replace("-", "") for r in canonical_rel if r.get("id")]

    # --- Lipid Signatures: Read from relation and fetch signature metadata ---
    signature_ids = _get_relation_ids(props, "Lipid Signatures")
    
    # Create reverse links from signature pages to this literature page
    for sig_id in signature_ids:
        enforce_signature_reverse_link(sig_id, parent_page_id)
    
    sig_meta = _collect_signature_metadata(signature_ids)
    
    # Get legacy multi-select values for backward compatibility (if any)
    lipid_signatures_legacy = _get_multi_names(props, "Lipid Signatures")
    lipid_signature_role_legacy = _get_multi_names(props, "Lipid Signature Role")
    
    # Use signature Short IDs (prefer relation-based, fall back to legacy multi-select)
    lipid_signatures = sig_meta["sig_short_ids"] if sig_meta["sig_short_ids"] else lipid_signatures_legacy
    lipid_signature_role = sig_meta["sig_roles"] if sig_meta["sig_roles"] else lipid_signature_role_legacy

    # --- Simple doc_type inference from Zotero item_type and tags ---
    doc_type = "JournalArticle"
    if item.item_type.lower() == "book":
        doc_type = "Book"
    elif item.item_type.lower() == "thesis":
        doc_type = "Thesis"
    elif any("review" in t.lower() for t in item.tags):
        doc_type = "Review"
    elif any("trial" in t.lower() for t in item.tags):
        doc_type = "ClinicalTrialDesign"

    # Extract year from date string if available
    year_val = None
    if item.date:
        m = re.search(r"(19|20)\d{2}", item.date)
        if m:
            year_val = int(m.group(0))

    doc_meta: Dict[str, Any] = {
        "doc_id": f"ZOTERO:{item.key}",
        "doc_source": "Zotero",
        "doc_source_subtype": "Attachment",  # notes will override if needed
        "doc_type": doc_type,
        "diseases": diseases,
        "targets": targets,
        "modality": modality,
        "stage": stage,
        "model_systems": model_systems,
        "biomarker_role": biomarker_role,
        "year": year_val,
        "journal": item.journal,
        "importance": importance,
        "manual_tags": item.tags,
    }

    # Extend phenotype_axes with signature axes
    all_phenotype_axes = sorted(
        set(phenotype_axes + sig_meta["sig_axes"])
    )
    
    lipid_meta: Dict[str, Any] = {
        "lipids_raw": lipids_raw,
        "lipids": canonical_lipids,
        "lipid_classes": [],  # can be filled later by deref'ing Lipid Species DB
        "lipid_signatures": lipid_signatures,  # List of Short IDs from signature pages
        "lipid_signature_role": lipid_signature_role,  # Aggregated biomarker roles
        "phenotype_axes": all_phenotype_axes,  # Extended with signature axes
        "signature_ownership": sig_meta["sig_ownership"],  # Aggregated ownership values
        "matrix": matrix,
        "treatment_arms": treatment_arms,
    }

    return {**doc_meta, **lipid_meta}


def get_email_semantic_metadata(email_page: Dict[str, Any]) -> Dict[str, Any]:
    """
    Read semantic + lipid metadata from an Email DB page and build a
    doc-level + lipid-level metadata dict.

    Expects (but does not require) these Notion properties on the Email DB:
      - Disease (multi-select)
      - Targets (multi-select)
      - Modality (multi-select)
      - Stage (select)
      - Model Systems (multi-select)
      - Biomarker Role (multi-select)
      - Importance (number)
      - Lipid Signatures (relation) – links to Lipid Signatures DB pages
      - Phenotype Axes (multi-select)
      - Matrix (multi-select)
      - Treatment Arms (multi-select)
      - Lipid Species (raw) (multi-select)
      - Canonical Lipid Species (relation) – optional
      
    For Lipid Signatures relation: Fetches linked signature pages and extracts:
      - Short ID (rich_text or title)
      - Biomarker Role (multi-select)
      - Phenotype Axes (multi-select)
      - Data Ownership (select)
      
    Args:
        email_page: Notion page object from Email DB query
        
    Returns:
        Dictionary with semantic and lipid metadata fields
    """
    props = email_page.get("properties", {}) or {}

    diseases = _get_multi_names(props, "Disease")
    targets = _get_multi_names(props, "Targets")
    modality = _get_multi_names(props, "Modality")
    stage = _get_select_name(props, "Stage")
    model_systems = _get_multi_names(props, "Model Systems")
    biomarker_role = _get_multi_names(props, "Biomarker Role")

    importance = props.get("Importance", {}).get("number")

    phenotype_axes = _get_multi_names(props, "Phenotype Axes")
    matrix = _get_multi_names(props, "Matrix")
    treatment_arms = _get_multi_names(props, "Treatment Arms")

    lipids_raw = _get_multi_names(props, "Lipid Species (raw)")
    canonical_rel = props.get("Canonical Lipid Species", {}).get("relation", []) or []
    canonical_lipids = [r.get("id", "").replace("-", "") for r in canonical_rel if r.get("id")]

    # --- Lipid Signatures: Read from relation and fetch signature metadata ---
    signature_ids = _get_relation_ids(props, "Lipid Signatures")
    
    # Create reverse links from signature pages to this email page
    email_page_id = email_page.get("id", "").replace("-", "")
    if email_page_id:
        for sig_id in signature_ids:
            enforce_signature_reverse_link(sig_id, email_page_id)
    
    sig_meta = _collect_signature_metadata(signature_ids)
    
    # Get legacy multi-select values for backward compatibility (if any)
    lipid_signatures_legacy = _get_multi_names(props, "Lipid Signatures")
    lipid_signature_role_legacy = _get_multi_names(props, "Lipid Signature Role")
    
    # Use signature Short IDs (prefer relation-based, fall back to legacy multi-select)
    lipid_signatures = sig_meta["sig_short_ids"] if sig_meta["sig_short_ids"] else lipid_signatures_legacy
    lipid_signature_role = sig_meta["sig_roles"] if sig_meta["sig_roles"] else lipid_signature_role_legacy

    # Extend phenotype_axes with signature axes
    all_phenotype_axes = sorted(
        set(phenotype_axes + sig_meta["sig_axes"])
    )

    # doc-level identity & type will be filled in ingest_email()
    base_meta: Dict[str, Any] = {
        "diseases": diseases,
        "targets": targets,
        "modality": modality,
        "stage": stage,
        "model_systems": model_systems,
        "biomarker_role": biomarker_role,
        "importance": importance,
        "lipids_raw": lipids_raw,
        "lipids": canonical_lipids,
        "lipid_classes": [],
        "lipid_signatures": lipid_signatures,
        "lipid_signature_role": lipid_signature_role,
        "phenotype_axes": all_phenotype_axes,
        "signature_ownership": sig_meta["sig_ownership"],
        "matrix": matrix,
        "treatment_arms": treatment_arms,
    }

    return base_meta


def get_experiment_semantic_metadata(exp_page: Dict[str, Any]) -> Dict[str, Any]:
    """
    Extract semantic + lipid + signature metadata from an Experiments DB page.

    Expected properties on Experiments DB:
      - Disease (multi-select)
      - Targets (multi-select)
      - Modality (multi-select)
      - Stage (select)
      - Model Systems (multi-select)
      - Biomarker Role (multi-select)
      - Phenotype Axes (multi-select)
      - Matrix (multi-select)
      - Treatment Arms (multi-select)
      - Lipid Signatures (relation -> Lipid Signatures DB)

    For Lipid Signatures relation: Fetches linked signature pages and extracts:
      - Short ID (rich_text or title)
      - Biomarker Role (multi-select)
      - Phenotype Axes (multi-select)
      - Data Ownership (select)

    Args:
        exp_page: Notion page object from Experiments DB query

    Returns:
        Dictionary with semantic and lipid signature metadata fields
    """
    props = exp_page.get("properties", {}) or {}

    diseases = _get_multi_names(props, "Disease")
    targets = _get_multi_names(props, "Targets")
    modality = _get_multi_names(props, "Modality")
    stage = _get_select_name(props, "Stage")
    model_systems = _get_multi_names(props, "Model Systems")
    biomarker_role = _get_multi_names(props, "Biomarker Role")

    phenotype_axes = _get_multi_names(props, "Phenotype Axes")
    matrix = _get_multi_names(props, "Matrix")
    treatment_arms = _get_multi_names(props, "Treatment Arms")

    # --- Lipid Signatures: Read from relation and fetch signature metadata ---
    signature_ids = _get_relation_ids(props, "Lipid Signatures")
    
    # Create reverse links from signature pages to this experiment page
    exp_page_id = exp_page.get("id", "").replace("-", "")
    if exp_page_id:
        for sig_id in signature_ids:
            enforce_signature_reverse_link(sig_id, exp_page_id)
    
    sig_meta = _collect_signature_metadata(signature_ids)
    
    # Get legacy multi-select values for backward compatibility (if any)
    lipid_signatures_legacy = _get_multi_names(props, "Lipid Signatures")
    lipid_signature_role_legacy = _get_multi_names(props, "Lipid Signature Role")
    
    # Use signature Short IDs (prefer relation-based, fall back to legacy multi-select)
    lipid_signatures = sig_meta["sig_short_ids"] if sig_meta["sig_short_ids"] else lipid_signatures_legacy
    lipid_signature_role = sig_meta["sig_roles"] if sig_meta["sig_roles"] else lipid_signature_role_legacy

    # Merge axes from signatures
    all_axes = sorted(
        set(phenotype_axes + sig_meta["sig_axes"])
    )

    base_meta: Dict[str, Any] = {
        "diseases": diseases,
        "targets": targets,
        "modality": modality,
        "stage": stage,
        "model_systems": model_systems,
        "biomarker_role": biomarker_role,
        "phenotype_axes": all_axes,
        "matrix": matrix,
        "treatment_arms": treatment_arms,
        "lipid_signatures": lipid_signatures,
        "lipid_signature_role": lipid_signature_role,
        "signature_ownership": sig_meta["sig_ownership"],
    }

    return base_meta


def get_dataset_semantic_metadata(dataset_page: Dict[str, Any]) -> Dict[str, Any]:
    """
    Extract semantic + lipid + signature metadata from a Dataset DB page.

    Expected properties on Dataset DB:
      - Disease (multi-select)
      - Matrix (multi-select)
      - Model Systems (multi-select)
      - Dataset Source Type (select)
      - Data Origin (select)
      - Related Signature(s) (relation -> Lipid Signatures DB)

    For Related Signature(s) relation: Fetches linked signature pages and extracts:
      - Short ID (rich_text or title)
      - Biomarker Role (multi-select)
      - Phenotype Axes (multi-select)
      - Data Ownership (select)

    Args:
        dataset_page: Notion page object from Dataset DB query

    Returns:
        Dictionary with semantic and lipid signature metadata fields
    """
    props = dataset_page.get("properties", {}) or {}

    diseases = _get_multi_names(props, "Disease")
    matrix = _get_multi_names(props, "Matrix")
    model_systems = _get_multi_names(props, "Model Systems")
    dataset_type = _get_select_name(props, "Dataset Source Type")
    dataset_origin = _get_select_name(props, "Data Origin")

    # --- Lipid Signatures: Read from relation and fetch signature metadata ---
    signature_ids = _get_relation_ids(props, "Related Signature(s)")
    
    # Create reverse links from signature pages to this dataset page
    dataset_page_id = dataset_page.get("id", "").replace("-", "")
    if dataset_page_id:
        for sig_id in signature_ids:
            enforce_signature_reverse_link(sig_id, dataset_page_id)
    
    sig_meta = _collect_signature_metadata(signature_ids)
    
    # Get legacy multi-select values for backward compatibility (if any)
    lipid_signatures_legacy = _get_multi_names(props, "Lipid Signatures")
    lipid_signature_role_legacy = _get_multi_names(props, "Lipid Signature Role")
    
    # Use signature Short IDs (prefer relation-based, fall back to legacy multi-select)
    lipid_signatures = sig_meta["sig_short_ids"] if sig_meta["sig_short_ids"] else lipid_signatures_legacy
    lipid_signature_role = sig_meta["sig_roles"] if sig_meta["sig_roles"] else lipid_signature_role_legacy

    # Merge axes from signatures
    phenotype_axes = sig_meta["sig_axes"]  # Datasets may not have their own axes field

    base_meta: Dict[str, Any] = {
        "diseases": diseases,
        "matrix": matrix,
        "model_systems": model_systems,
        "dataset_origin": dataset_origin,
        "dataset_type": dataset_type,
        "lipid_signatures": lipid_signatures,
        "lipid_signature_role": lipid_signature_role,
        "phenotype_axes": sorted(set(phenotype_axes)) if phenotype_axes else [],
        "signature_ownership": sig_meta["sig_ownership"],
    }

    return base_meta