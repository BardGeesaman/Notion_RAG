"""
Literature metadata extraction.

Extracts semantic and lipid metadata from Literature DB pages.
"""

from __future__ import annotations

import re
from typing import Any, Dict

from amprenta_rag.ingestion.metadata.helpers import (
    fetch_notion_page,
    get_multi_names,
    get_relation_ids,
    get_select_name,
)
from amprenta_rag.ingestion.metadata.signature_metadata import (
    collect_signature_metadata,
    enforce_signature_reverse_link,
)
from amprenta_rag.ingestion.zotero_api import ZoteroItem


def get_literature_semantic_metadata(
    parent_page_id: str, item: ZoteroItem
) -> Dict[str, Any]:
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
    page = fetch_notion_page(parent_page_id)
    props = page.get("properties", {}) or {}

    # --- Basic semantic fields ---
    diseases = get_multi_names(props, "Disease")
    targets = get_multi_names(props, "Targets")
    modality = get_multi_names(props, "Modality")
    stage = get_select_name(props, "Stage")
    model_systems = get_multi_names(props, "Model Systems")
    biomarker_role = get_multi_names(props, "Biomarker Role")

    importance = props.get("Importance", {}).get("number")

    phenotype_axes = get_multi_names(props, "Phenotype Axes")
    matrix = get_multi_names(props, "Matrix")
    treatment_arms = get_multi_names(props, "Treatment Arms")

    # --- Lipid-related fields ---
    lipids_raw = get_multi_names(props, "Lipid Species (raw)")

    canonical_rel = props.get("Canonical Lipid Species", {}).get("relation", []) or []
    canonical_lipids = [
        r.get("id", "").replace("-", "") for r in canonical_rel if r.get("id")
    ]

    # --- Lipid Signatures: Read from relation and fetch signature metadata ---
    signature_ids = get_relation_ids(props, "Lipid Signatures")

    # Create reverse links from signature pages to this literature page
    for sig_id in signature_ids:
        enforce_signature_reverse_link(sig_id, parent_page_id)

    sig_meta = collect_signature_metadata(signature_ids)

    # Get legacy multi-select values for backward compatibility (if any)
    lipid_signatures_legacy = get_multi_names(props, "Lipid Signatures")
    lipid_signature_role_legacy = get_multi_names(props, "Lipid Signature Role")

    # Use signature Short IDs (prefer relation-based, fall back to legacy multi-select)
    lipid_signatures = (
        sig_meta["sig_short_ids"]
        if sig_meta["sig_short_ids"]
        else lipid_signatures_legacy
    )
    lipid_signature_role = (
        sig_meta["sig_roles"] if sig_meta["sig_roles"] else lipid_signature_role_legacy
    )

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
    all_phenotype_axes = sorted(set(phenotype_axes + sig_meta["sig_axes"]))

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

