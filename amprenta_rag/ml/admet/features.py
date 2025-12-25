"""Feature name mapping for ADMET models.

ADMET feature vector schema (dim=2054):
- 0..2047: Morgan fingerprint bit positions
- 2048..2053: RDKit descriptors (6)
"""

from __future__ import annotations

from typing import List


DESCRIPTOR_NAMES: List[str] = [
    "MolWt",
    "MolLogP",
    "TPSA",
    "NumHDonors",
    "NumHAcceptors",
    "NumRotatableBonds",
]

FEATURE_NAMES: List[str] = [f"MorganBit_{i:04d}" for i in range(2048)] + list(DESCRIPTOR_NAMES)


def get_feature_names() -> List[str]:
    """Return the full ordered list of feature names (len=2054)."""
    return FEATURE_NAMES


def get_feature_name(index: int) -> str:
    """Return the feature name for a single index (0 <= index < 2054)."""
    if index < 0 or index >= len(FEATURE_NAMES):
        raise IndexError(f"feature index out of range: {index}")
    return FEATURE_NAMES[index]


__all__ = ["FEATURE_NAMES", "DESCRIPTOR_NAMES", "get_feature_names", "get_feature_name"]


