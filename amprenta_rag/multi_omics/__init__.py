"""Multi-omics latent factor utilities (data extraction + MOFA runner)."""

from .data_extractor import align_samples, extract_feature_matrix
from .mofa_runner import run_mofa

__all__ = ["extract_feature_matrix", "align_samples", "run_mofa"]


