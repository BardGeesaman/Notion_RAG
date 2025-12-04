"""
Experiment metadata extraction.

Extracts semantic and lipid metadata from Experiments DB pages.
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

    diseases = get_multi_names(props, "Disease")
    targets = get_multi_names(props, "Targets")
    modality = get_multi_names(props, "Modality")
    stage = get_select_name(props, "Stage")
    model_systems = get_multi_names(props, "Model Systems")
    biomarker_role = get_multi_names(props, "Biomarker Role")

    phenotype_axes = get_multi_names(props, "Phenotype Axes")
    matrix = get_multi_names(props, "Matrix")
    treatment_arms = get_multi_names(props, "Treatment Arms")

    # --- Lipid Signatures: Read from relation and fetch signature metadata ---
    signature_ids = get_relation_ids(props, "Lipid Signatures")

    # Create reverse links from signature pages to this experiment page
    exp_page_id = exp_page.get("id", "").replace("-", "")
    if exp_page_id:
        for sig_id in signature_ids:
            enforce_signature_reverse_link(sig_id, exp_page_id)

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
    all_axes = sorted(set(phenotype_axes + sig_meta["sig_axes"]))

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

