"""
Data models for signature matching.

Contains dataclasses and types used throughout the signature matching system.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List

from amprenta_rag.signatures.signature_scoring import SignatureScoreResult


@dataclass
class SignatureMatchResult:
    """Result of matching a signature against a dataset."""

    signature_page_id: str
    signature_name: str
    score: float
    overlap_fraction: float
    matched_components: List[str]
    missing_components: List[str]
    conflicting_components: List[str]
    score_result: SignatureScoreResult

