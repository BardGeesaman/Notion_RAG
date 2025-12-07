"""
Program-level signature mapping analysis.

This package contains functions for computing program-signature scoring matrices,
coverage analysis, and convergence indicators.
"""

from amprenta_rag.analysis.program_maps.convergence import compute_convergence_indicators
from amprenta_rag.analysis.program_maps.coverage import compute_program_omics_coverage
from amprenta_rag.analysis.program_maps.models import ProgramOmicsCoverage, ProgramSignatureMap, ProgramSignatureScore
from amprenta_rag.analysis.program_maps.reporting import (
    generate_program_map_report,
    generate_program_signature_map,
    update_notion_with_program_map,
)
from amprenta_rag.analysis.program_maps.scoring import compute_program_signature_scores, get_program_datasets

__all__ = [
    "ProgramSignatureScore",
    "ProgramOmicsCoverage",
    "ProgramSignatureMap",
    "get_program_datasets",
    "compute_program_signature_scores",
    "compute_program_omics_coverage",
    "compute_convergence_indicators",
    "generate_program_signature_map",
    "generate_program_map_report",
    "update_notion_with_program_map",
]
