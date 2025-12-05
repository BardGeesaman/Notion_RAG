"""
Evidence reporting module.

This package provides functionality for generating automated evidence reports
for Programs, Experiments, Datasets, Signatures, and Features.
"""

from __future__ import annotations

from amprenta_rag.reporting.evidence_report import (
    EvidenceReport,
    format_evidence_report,
    generate_dataset_evidence_report,
    generate_experiment_evidence_report,
    generate_feature_evidence_report,
    generate_program_evidence_report,
    generate_signature_evidence_report,
    write_evidence_report_to_notion,
)

__all__ = [
    "EvidenceReport",
    "generate_program_evidence_report",
    "generate_experiment_evidence_report",
    "generate_dataset_evidence_report",
    "generate_signature_evidence_report",
    "generate_feature_evidence_report",
    "format_evidence_report",
    "write_evidence_report_to_notion",
]

