"""
Automatic signature matching and scoring for datasets.

This package provides:
- Signature scoring against datasets
- Automatic signature detection during dataset ingestion
- Notion writebacks for signature matches
- Lipidomics-aware species mapping

Maintains backward compatibility by re-exporting all public functions.
"""

from __future__ import annotations

from amprenta_rag.ingestion.signature_matching.matching import (
    find_matching_signatures_for_dataset,
    score_signature_against_dataset,
)
from amprenta_rag.ingestion.signature_matching.models import SignatureMatchResult
from amprenta_rag.ingestion.signature_matching.species_mapping import (
    map_raw_lipid_to_canonical_species,
)
from amprenta_rag.ingestion.signature_matching.writeback import (
    update_dataset_with_signature_matches,
)

__all__ = [
    "SignatureMatchResult",
    "map_raw_lipid_to_canonical_species",
    "find_matching_signatures_for_dataset",
    "score_signature_against_dataset",
    "update_dataset_with_signature_matches",
]

