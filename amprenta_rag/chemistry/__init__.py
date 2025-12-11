"""
Chemistry and HTS integration module.

This package provides functionality for managing chemical compounds,
HTS campaigns, and biochemical assay results.
"""

from __future__ import annotations

from amprenta_rag.chemistry.database import (
    get_hits_for_campaign,
    initialize_chemistry_database,
    find_compound_by_inchi_key,
    insert_biochemical_results,
    insert_compound,
    insert_hts_campaign,
    insert_hts_results,
    insert_compound_signature_link,
    get_compounds_for_signature,
    get_signatures_for_compound,
)
from amprenta_rag.chemistry.normalization import (
    compute_molecular_descriptors,
    generate_compound_id,
    normalize_smiles,
)
from amprenta_rag.chemistry.schema import (
    BiochemicalResult,
    Compound,
    HTSCampaign,
    HTSResult,
    CompoundSignatureLink,
)

__all__ = [
    "Compound",
    "HTSCampaign",
    "HTSResult",
    "BiochemicalResult",
    "CompoundSignatureLink",
    "initialize_chemistry_database",
    "insert_compound",
    "find_compound_by_inchi_key",
    "insert_hts_campaign",
    "insert_hts_results",
    "insert_biochemical_results",
    "get_hits_for_campaign",
    "insert_compound_signature_link",
    "get_compounds_for_signature",
    "get_signatures_for_compound",
    "normalize_smiles",
    "compute_molecular_descriptors",
    "generate_compound_id",
]
