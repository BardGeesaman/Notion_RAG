# amprenta_rag/ingestion/metadata_semantic.py

from __future__ import annotations

from typing import Dict, Any, List, Optional

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
      - Lipid Signatures (multi-select)
      - Lipid Signature Role (multi-select)
      - Phenotype Axes (multi-select)
      - Matrix (multi-select)
      - Treatment Arms (multi-select)
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

    lipid_signatures = _get_multi_names(props, "Lipid Signatures")
    lipid_signature_role = _get_multi_names(props, "Lipid Signature Role")

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

    lipid_meta: Dict[str, Any] = {
        "lipids_raw": lipids_raw,
        "lipids": canonical_lipids,
        "lipid_classes": [],  # can be filled later by deref'ing Lipid Species DB
        "lipid_signatures": lipid_signatures,
        "lipid_signature_role": lipid_signature_role,
        "phenotype_axes": phenotype_axes,
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
      - Lipid Signatures (multi-select)
      - Lipid Signature Role (multi-select)
      - Phenotype Axes (multi-select)
      - Matrix (multi-select)
      - Treatment Arms (multi-select)
      - Lipid Species (raw) (multi-select)
      - Canonical Lipid Species (relation) – optional
      
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

    lipid_signatures = _get_multi_names(props, "Lipid Signatures")
    lipid_signature_role = _get_multi_names(props, "Lipid Signature Role")

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
        "phenotype_axes": phenotype_axes,
        "matrix": matrix,
        "treatment_arms": treatment_arms,
    }

    return base_meta