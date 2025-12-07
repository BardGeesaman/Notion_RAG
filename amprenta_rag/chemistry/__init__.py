"""
Chemistry and HTS integration module.

This package provides functionality for managing chemical compounds,
HTS campaigns, and biochemical assay results.
"""

from __future__ import annotations

from amprenta_rag.chemistry.database import (
    find_compound_by_inchi_key,
    get_chemistry_db_path,
    get_hits_for_campaign,
    initialize_chemistry_database,
    insert_biochemical_results,
    insert_compound,
    insert_hts_campaign,
    insert_hts_results,
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
)

__all__ = [
    "Compound",
    "HTSCampaign",
    "HTSResult",
    "BiochemicalResult",
    "initialize_chemistry_database",
    "get_chemistry_db_path",
    "insert_compound",
    "find_compound_by_inchi_key",
    "insert_hts_campaign",
    "insert_hts_results",
    "insert_biochemical_results",
    "get_hits_for_campaign",
    "normalize_smiles",
    "compute_molecular_descriptors",
    "generate_compound_id",
]
