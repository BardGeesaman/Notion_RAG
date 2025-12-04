"""
Dataset metadata extraction.

Extracts semantic and lipid metadata from Dataset DB pages.
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

    diseases = get_multi_names(props, "Disease")
    matrix = get_multi_names(props, "Matrix")
    model_systems = get_multi_names(props, "Model Systems")
    dataset_type = get_select_name(props, "Dataset Source Type")
    dataset_origin = get_select_name(props, "Data Origin")

    # --- Lipid Signatures: Read from relation and fetch signature metadata ---
    signature_ids = get_relation_ids(props, "Related Signature(s)")

    # Create reverse links from signature pages to this dataset page
    dataset_page_id = dataset_page.get("id", "").replace("-", "")
    if dataset_page_id:
        for sig_id in signature_ids:
            enforce_signature_reverse_link(sig_id, dataset_page_id)

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

