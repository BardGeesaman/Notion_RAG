"""
Chemistry and HTS integration module.

This package provides functionality for managing chemical compounds,
HTS campaigns, and biochemical assay results.
"""

from __future__ import annotations

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
    "normalize_smiles",
    "compute_molecular_descriptors",
    "generate_compound_id",
]
