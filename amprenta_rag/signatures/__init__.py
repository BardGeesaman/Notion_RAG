"""
Signature scoring engine for lipid signatures.

This module provides tools for:
- Matching lipid species between datasets and signatures
- Loading signature definitions from Notion or TSV
- Computing direction-aware, weight-aware signature scores
"""

from .signature_loader import SignatureComponent, load_signature_from_tsv
from .signature_scoring import SignatureScoreResult, score_signature
from .species_matching import match_species, normalize_species_name

__all__ = [
    "match_species",
    "normalize_species_name",
    "load_signature_from_tsv",
    "SignatureComponent",
    "score_signature",
    "SignatureScoreResult",
]
