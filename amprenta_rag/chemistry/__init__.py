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
from amprenta_rag.chemistry.notion_integration import (
    create_biochemical_hit_page,
    create_compound_feature_page,
    create_hts_campaign_page,
    find_or_create_compound_page,
    find_or_create_hts_campaign_page,
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
    "create_compound_feature_page",
    "find_or_create_compound_page",
    "create_hts_campaign_page",
    "find_or_create_hts_campaign_page",
    "create_biochemical_hit_page",
]
