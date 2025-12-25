"""Biomarker discovery module with stable feature selection methods."""

from amprenta_rag.ml.biomarker.datasets import load_biomarker_dataset
from amprenta_rag.ml.biomarker.discovery import BiomarkerDiscoveryService
from amprenta_rag.ml.biomarker.importance import CVFeatureImportance
from amprenta_rag.ml.biomarker.stability import StabilitySelector
from amprenta_rag.ml.biomarker.statistical import (
    t_test_selection,
    mann_whitney_selection,
    anova_selection,
    fdr_correction,
)

__all__ = [
    "StabilitySelector",
    "CVFeatureImportance",
    "BiomarkerDiscoveryService",
    "load_biomarker_dataset",
    "t_test_selection",
    "mann_whitney_selection",
    "anova_selection",
    "fdr_correction",
]

