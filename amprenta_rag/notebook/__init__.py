"""Notebook copilot module for AI-assisted cell generation."""

from amprenta_rag.notebook.context import AnalysisContext, generate_context_cell
from amprenta_rag.notebook.helpers import (
    load_dataset,
    load_experiment,
    load_campaign,
    load_compound,
    run_hts_qc,
    fit_dose_response,
    publish_to_rag,
)

__all__ = [
    "AnalysisContext",
    "generate_context_cell",
    "load_dataset",
    "load_experiment",
    "load_campaign",
    "load_compound",
    "run_hts_qc",
    "fit_dose_response",
    "publish_to_rag",
]


