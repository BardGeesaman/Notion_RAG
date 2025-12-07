"""
Common feature types for ingestion.
"""

from __future__ import annotations

from enum import Enum


class FeatureType(str, Enum):
    """Enum for feature types."""
    GENE = "gene"
    PROTEIN = "protein"
    METABOLITE = "metabolite"
    LIPID = "lipid"
    OTHER = "other"

