"""
Email metadata extraction.

Extracts semantic and lipid metadata from Email DB pages.
"""

from __future__ import annotations

from typing import Any, Dict

from amprenta_rag.ingestion.metadata.helpers import (
    get_multi_names,
    get_relation_ids,
    get_select_name,
)
from amprenta_rag.ingestion.metadata.signature_metadata import (
    collect_signature_metadata,
    enforce_signature_reverse_link,
)


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

    lipids_raw = get_multi_names(props, "Lipid Species (raw)")
    canonical_rel = props.get("Canonical Lipid Species", {}).get("relation", []) or []
    canonical_lipids = [
        r.get("id", "").replace("-", "") for r in canonical_rel if r.get("id")
    ]

    # --- Lipid Signatures: Read from relation and fetch signature metadata ---
    signature_ids = get_relation_ids(props, "Lipid Signatures")

    # Create reverse links from signature pages to this email page
    email_page_id = email_page.get("id", "").replace("-", "")
    if email_page_id:
        for sig_id in signature_ids:
            enforce_signature_reverse_link(sig_id, email_page_id)

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

    # Extend phenotype_axes with signature axes
    all_phenotype_axes = sorted(set(phenotype_axes + sig_meta["sig_axes"]))

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

