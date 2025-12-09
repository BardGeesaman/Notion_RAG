"""
Evidence Report Engine.

Generates comprehensive cross-omics evidence reports for Programs, Experiments,
Datasets, Signatures, and Features using RAG and cross-omics reasoning.

Postgres is now the source of truth - Notion support has been removed.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime
from typing import Any, Dict, List, Optional
from uuid import UUID

from amprenta_rag.database.base import get_db
from amprenta_rag.database.models import Dataset, Program, Signature
from amprenta_rag.domain.analytics import EvidenceReport, EvidenceSection
from amprenta_rag.logging_utils import get_logger
from amprenta_rag.query.cross_omics import (
    cross_omics_dataset_summary_postgres,
    cross_omics_feature_summary_postgres,
    cross_omics_program_summary_postgres,
    cross_omics_signature_summary_postgres,
)
from amprenta_rag.signatures.signature_validation import validate_signature_against_all_datasets

logger = get_logger(__name__)


@dataclass
class EvidenceReportLegacy:
    """
    Legacy evidence report structure for backward compatibility.

    Attributes:
        entity_type: Type of entity ("program", "experiment", "dataset", "signature", "feature")
        entity_id: ID of the entity
        entity_name: Name of the entity
        summary: Generated summary text (markdown format)
        metadata: Additional metadata (disease, model system, matrix, etc.)
    """

    entity_type: str
    entity_id: str
    entity_name: str
    summary: str
    metadata: Dict[str, Any] = field(default_factory=dict)


def generate_program_report(program_id: UUID) -> EvidenceReport:
    """Generate an evidence report for a Program from Postgres."""
    db = next(get_db())
    program = db.query(Program).filter(Program.id == program_id).first()
    datasets = getattr(program, "datasets", [])
    sigs = getattr(program, "signatures", []) if hasattr(program, "signatures") else []
    sections = []
    
    # Section: Overview
    overview = EvidenceSection(
        title="Overview",
        summary_text=f"Program: {program.name if program else program_id}",
        supporting_datasets=[d.id for d in datasets],
        key_features=None,
        signatures=[s.id for s in sigs],
        references=None,
    )
    sections.append(overview)
    
    return EvidenceReport(
        entity_id=program_id,
        entity_type="program",
        generated_at=datetime.utcnow(),
        sections=sections,
    )


def generate_dataset_report(dataset_id: UUID) -> EvidenceReport:
    """Generate an evidence report for a Dataset from Postgres."""
    db = next(get_db())
    d = db.query(Dataset).filter(Dataset.id == dataset_id).first()
    
    sections = [
        EvidenceSection(
            title="Overview",
            summary_text=f"Dataset: {d.name if d else dataset_id}",
            supporting_datasets=None,
            key_features=None,
            signatures=None,
            references=None,
        )
    ]
    return EvidenceReport(
        entity_id=dataset_id,
        entity_type="dataset",
        generated_at=datetime.utcnow(),
        sections=sections,
    )


def generate_signature_report(signature_id: UUID) -> EvidenceReport:
    """Generate an evidence report for a Signature from Postgres."""
    db = next(get_db())
    s = db.query(Signature).filter(Signature.id == signature_id).first()
    
    sections = [
        EvidenceSection(
            title="Overview",
            summary_text=f"Signature: {getattr(s, 'name', signature_id)}",
            supporting_datasets=None,
            key_features=None,
            signatures=None,
            references=None,
        )
    ]
    
    val = validate_signature_against_all_datasets(signature_id)
    sections.append(
        EvidenceSection(
            title="Validation Metrics",
            summary_text=val.summary,
            supporting_datasets=val.matched_dataset_ids,
            key_features=None,
            signatures=None,
            references=None,
        )
    )
    return EvidenceReport(
        entity_id=signature_id,
        entity_type="signature",
        generated_at=datetime.utcnow(),
        sections=sections,
    )


def generate_program_evidence_report(
    program_page_id: str,
    top_k_per_omics: int = 20,
) -> EvidenceReportLegacy:
    """
    Generate an evidence report for a Program.

    Args:
        program_page_id: Page ID of program
        top_k_per_omics: Number of top chunks per omics type

    Returns:
        EvidenceReportLegacy object
    """
    logger.info(
        "[REPORTING][EVIDENCE] Generating evidence report for program %s",
        program_page_id,
    )

    # Use existing cross-omics summary function
    summary = cross_omics_program_summary_postgres(
        program_page_id=program_page_id,
        top_k_per_omics=top_k_per_omics,
    )

    program_name = f"Program {program_page_id[:8]}"

    return EvidenceReportLegacy(
        entity_type="program",
        entity_id=program_page_id,
        entity_name=program_name,
        summary=summary,
        metadata={},
    )


def generate_experiment_evidence_report(
    experiment_page_id: str,
    top_k_chunks: int = 100,
) -> EvidenceReportLegacy:
    """
    Generate an evidence report for an Experiment.

    Args:
        experiment_page_id: Page ID of experiment
        top_k_chunks: Number of top chunks to retrieve

    Returns:
        EvidenceReportLegacy object
    """
    logger.info(
        "[REPORTING][EVIDENCE] Generating evidence report for experiment %s",
        experiment_page_id,
    )

    summary = cross_omics_dataset_summary_postgres(
        dataset_page_id=experiment_page_id,
        top_k_chunks=top_k_chunks,
    )

    experiment_name = f"Experiment {experiment_page_id[:8]}"

    return EvidenceReportLegacy(
        entity_type="experiment",
        entity_id=experiment_page_id,
        entity_name=experiment_name,
        summary=summary,
        metadata={},
    )


def generate_dataset_evidence_report(
    dataset_page_id: str,
    top_k_chunks: int = 100,
) -> EvidenceReportLegacy:
    """
    Generate an evidence report for a Dataset.

    Args:
        dataset_page_id: Page ID of dataset
        top_k_chunks: Number of top chunks to retrieve

    Returns:
        EvidenceReportLegacy object
    """
    logger.info(
        "[REPORTING][EVIDENCE] Generating evidence report for dataset %s",
        dataset_page_id,
    )

    summary = cross_omics_dataset_summary_postgres(
        dataset_page_id=dataset_page_id,
        top_k_chunks=top_k_chunks,
    )

    dataset_name = f"Dataset {dataset_page_id[:8]}"

    return EvidenceReportLegacy(
        entity_type="dataset",
        entity_id=dataset_page_id,
        entity_name=dataset_name,
        summary=summary,
        metadata={},
    )


def generate_signature_evidence_report(
    signature_page_id: str,
    top_k_datasets: int = 20,
    top_k_chunks: int = 100,
) -> EvidenceReportLegacy:
    """
    Generate an evidence report for a Signature.

    Args:
        signature_page_id: Page ID of signature
        top_k_datasets: Number of top matching datasets to include
        top_k_chunks: Number of top chunks per dataset

    Returns:
        EvidenceReportLegacy object
    """
    logger.info(
        "[REPORTING][EVIDENCE] Generating evidence report for signature %s",
        signature_page_id,
    )

    summary = cross_omics_signature_summary_postgres(
        signature_page_id=signature_page_id,
        top_k_datasets=top_k_datasets,
        top_k_chunks=top_k_chunks,
    )

    signature_name = f"Signature {signature_page_id[:8]}"

    return EvidenceReportLegacy(
        entity_type="signature",
        entity_id=signature_page_id,
        entity_name=signature_name,
        summary=summary,
        metadata={},
    )


def generate_feature_evidence_report(
    feature_name: str,
    feature_type: str,
    top_k_datasets: int = 20,
    top_k_chunks: int = 100,
) -> EvidenceReportLegacy:
    """
    Generate an evidence report for a Feature.

    Args:
        feature_name: Name of the feature
        feature_type: Type of feature ("gene", "protein", "metabolite", "lipid")
        top_k_datasets: Number of top datasets to include
        top_k_chunks: Number of top chunks per dataset

    Returns:
        EvidenceReportLegacy object
    """
    logger.info(
        "[REPORTING][EVIDENCE] Generating evidence report for %s feature '%s'",
        feature_type,
        feature_name,
    )

    summary = cross_omics_feature_summary_postgres(
        feature_name=feature_name,
        feature_type=feature_type,
        top_k_datasets=top_k_datasets,
        top_k_chunks=top_k_chunks,
    )

    return EvidenceReportLegacy(
        entity_type="feature",
        entity_id=feature_name,
        entity_name=feature_name,
        summary=summary,
        metadata={"feature_type": feature_type},
    )


def format_evidence_report(
    report: EvidenceReportLegacy,
    include_metadata: bool = True,
) -> str:
    """
    Format an evidence report as a complete markdown document.

    Args:
        report: EvidenceReportLegacy object
        include_metadata: Whether to include metadata section

    Returns:
        Formatted markdown document
    """
    doc = f"# Evidence Report: {report.entity_name}\n\n"
    doc += f"**Entity Type**: {report.entity_type.capitalize()}\n"
    doc += f"**Entity ID**: `{report.entity_id}`\n\n"

    if include_metadata and report.metadata:
        doc += "## Metadata\n\n"
        for key, value in report.metadata.items():
            if isinstance(value, list):
                doc += f"- **{key.replace('_', ' ').title()}**: {', '.join(str(v) for v in value)}\n"
            else:
                doc += f"- **{key.replace('_', ' ').title()}**: {value}\n"
        doc += "\n"

    doc += "## Summary\n\n"
    doc += report.summary

    return doc


def write_evidence_report_to_notion(
    report: EvidenceReportLegacy,
    target_page_id: Optional[str] = None,
    property_name: str = "Evidence Report",
) -> bool:
    """
    DEPRECATED: Notion support removed. Returns False.
    
    Previously wrote an evidence report to a Notion page.
    """
    logger.debug(
        "[REPORTING][EVIDENCE] write_evidence_report_to_notion() is a no-op (Notion removed)"
    )
    return False
