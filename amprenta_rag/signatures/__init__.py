"""
Signature scoring engine for lipid signatures.

This module provides tools for:
- Matching lipid species between datasets and signatures
- Loading signature definitions from Notion or TSV
- Computing direction-aware, weight-aware signature scores
"""

from .species_matching import match_species, normalize_species_name
from .signature_loader import load_signature_from_tsv, SignatureComponent
from .signature_scoring import score_signature, SignatureScoreResult

__all__ = [
    "match_species",
    "normalize_species_name",
    "load_signature_from_tsv",
    "SignatureComponent",
    "score_signature",
    "SignatureScoreResult",
]

