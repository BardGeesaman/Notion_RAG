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

__all__ = [
    "normalize_smiles",
    "compute_molecular_descriptors",
    "generate_compound_id",
]
