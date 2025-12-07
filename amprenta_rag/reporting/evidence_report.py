"""
Evidence Report Engine.

Generates comprehensive cross-omics evidence reports for Programs, Experiments,
Datasets, Signatures, and Features using RAG and cross-omics reasoning.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from typing import Optional
from uuid import UUID

from amprenta_rag.database.base import get_db
from amprenta_rag.database.models import Dataset, Program, Signature
from amprenta_rag.domain.analytics import EvidenceReport, EvidenceSection
from amprenta_rag.logging_utils import get_logger
from amprenta_rag.query.cross_omics import (
    cross_omics_dataset_summary,
    cross_omics_feature_summary,
    cross_omics_program_summary,
    cross_omics_signature_summary,
)
from amprenta_rag.signatures.signature_validation import validate_signature_against_all_datasets

# Optionally import or stub cross-omics/RAG helpers here

logger = get_logger(__name__)


@dataclass
class EvidenceReport:
    """
    Represents a generated evidence report.

    Attributes:
        entity_type: Type of entity ("program", "experiment", "dataset", "signature", "feature")
        entity_id: Notion page ID of the entity
        entity_name: Name of the entity
        summary: Generated summary text (markdown format)
        metadata: Additional metadata (disease, model system, matrix, etc.)
    """

    entity_type: str
    entity_id: str
    entity_name: str
    summary: str
    metadata: dict = None


def generate_program_report(program_id: UUID) -> EvidenceReport:
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
    # Optionally use cross-omics summary and RAG helpers for richer sections.
    # sections.append(...)
    return EvidenceReport(
        entity_id=program_id,
        entity_type="program",
        generated_at=datetime.utcnow(),
        sections=sections,
    )


def generate_dataset_report(dataset_id: UUID) -> EvidenceReport:
    db = next(get_db())
    d = db.query(Dataset).filter(Dataset.id == dataset_id).first()
    # Populate basic summary section
    # Optionally add: linked signatures, key features, RAG summaries
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
        entity_id=dataset_id, entity_type="dataset", generated_at=datetime.utcnow(), sections=sections
    )


def generate_signature_report(signature_id: UUID) -> EvidenceReport:
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
        entity_id=signature_id, entity_type="signature", generated_at=datetime.utcnow(), sections=sections
    )


def generate_program_evidence_report(
    program_page_id: str,
    top_k_per_omics: int = 20,
) -> EvidenceReport:
    """
    Generate an evidence report for a Program.

    Args:
        program_page_id: Notion page ID of program (with dashes)
        top_k_per_omics: Number of top chunks per omics type

    Returns:
        EvidenceReport object
    """
    logger.info(
        "[REPORTING][EVIDENCE] Generating evidence report for program %s",
        program_page_id,
    )

    # Use existing cross-omics summary function
    summary = cross_omics_program_summary(
        program_page_id=program_page_id,
        top_k_per_omics=top_k_per_omics,
    )

    # Extract program name
    program_name = f"Program {program_page_id[:8]}"
    try:
        from amprenta_rag.query.cross_omics.helpers import extract_text_property, fetch_notion_page

        page = fetch_notion_page(program_page_id)
        program_name = extract_text_property(page, "Program") or extract_text_property(page, "Name") or program_name
    except Exception as e:
        logger.debug("[REPORTING][EVIDENCE] Could not fetch program name: %r", e)

    # Extract metadata
    metadata = {}
    try:
        from amprenta_rag.query.cross_omics.helpers import extract_select_values, fetch_notion_page

        page = fetch_notion_page(program_page_id)
        props = page.get("properties", {}) or {}

        if "Disease" in props:
            metadata["disease"] = extract_select_values(props["Disease"])
        if "Model Systems" in props:
            metadata["model_systems"] = extract_select_values(props["Model Systems"])
    except Exception as e:
        logger.debug("[REPORTING][EVIDENCE] Could not extract metadata: %r", e)

    return EvidenceReport(
        entity_type="program",
        entity_id=program_page_id,
        entity_name=program_name,
        summary=summary,
        metadata=metadata,
    )


def generate_experiment_evidence_report(
    experiment_page_id: str,
    top_k_chunks: int = 100,
) -> EvidenceReport:
    """
    Generate an evidence report for an Experiment.

    Args:
        experiment_page_id: Notion page ID of experiment (with dashes)
        top_k_chunks: Number of top chunks to retrieve

    Returns:
        EvidenceReport object
    """
    logger.info(
        "[REPORTING][EVIDENCE] Generating evidence report for experiment %s",
        experiment_page_id,
    )

    # Use dataset summary function (experiments are similar to datasets)
    summary = cross_omics_dataset_summary(
        dataset_page_id=experiment_page_id,
        top_k_chunks=top_k_chunks,
    )

    # Extract experiment name
    experiment_name = f"Experiment {experiment_page_id[:8]}"
    try:
        from amprenta_rag.query.cross_omics.helpers import extract_text_property, fetch_notion_page

        page = fetch_notion_page(experiment_page_id)
        experiment_name = extract_text_property(page, "Title") or extract_text_property(page, "Name") or experiment_name
    except Exception as e:
        logger.debug("[REPORTING][EVIDENCE] Could not fetch experiment name: %r", e)

    # Extract metadata
    metadata = {}
    try:
        from amprenta_rag.query.cross_omics.helpers import extract_select_values, fetch_notion_page

        page = fetch_notion_page(experiment_page_id)
        props = page.get("properties", {}) or {}

        if "Disease" in props:
            metadata["disease"] = extract_select_values(props["Disease"])
        if "Matrix" in props:
            metadata["matrix"] = extract_select_values(props["Matrix"])
        if "Model Systems" in props:
            metadata["model_systems"] = extract_select_values(props["Model Systems"])
    except Exception as e:
        logger.debug("[REPORTING][EVIDENCE] Could not extract metadata: %r", e)

    return EvidenceReport(
        entity_type="experiment",
        entity_id=experiment_page_id,
        entity_name=experiment_name,
        summary=summary,
        metadata=metadata,
    )


def generate_dataset_evidence_report(
    dataset_page_id: str,
    top_k_chunks: int = 100,
) -> EvidenceReport:
    """
    Generate an evidence report for a Dataset.

    Args:
        dataset_page_id: Notion page ID of dataset (with dashes)
        top_k_chunks: Number of top chunks to retrieve

    Returns:
        EvidenceReport object
    """
    logger.info(
        "[REPORTING][EVIDENCE] Generating evidence report for dataset %s",
        dataset_page_id,
    )

    # Use existing cross-omics summary function
    summary = cross_omics_dataset_summary(
        dataset_page_id=dataset_page_id,
        top_k_chunks=top_k_chunks,
    )

    # Extract dataset name
    dataset_name = f"Dataset {dataset_page_id[:8]}"
    try:
        from amprenta_rag.query.cross_omics.helpers import extract_text_property, fetch_notion_page

        page = fetch_notion_page(dataset_page_id)
        dataset_name = (
            extract_text_property(page, "Experiment Name") or extract_text_property(page, "Name") or dataset_name
        )
    except Exception as e:
        logger.debug("[REPORTING][EVIDENCE] Could not fetch dataset name: %r", e)

    # Extract metadata
    metadata = {}
    try:
        from amprenta_rag.query.cross_omics.helpers import extract_select_values, fetch_notion_page

        page = fetch_notion_page(dataset_page_id)
        props = page.get("properties", {}) or {}

        if "Disease" in props:
            metadata["disease"] = extract_select_values(props["Disease"])
        if "Matrix" in props:
            metadata["matrix"] = extract_select_values(props["Matrix"])
        if "Model Systems" in props:
            metadata["model_systems"] = extract_select_values(props["Model Systems"])
        if "Omics Type" in props:
            metadata["omics_type"] = extract_select_values(props["Omics Type"])
    except Exception as e:
        logger.debug("[REPORTING][EVIDENCE] Could not extract metadata: %r", e)

    return EvidenceReport(
        entity_type="dataset",
        entity_id=dataset_page_id,
        entity_name=dataset_name,
        summary=summary,
        metadata=metadata,
    )


def generate_signature_evidence_report(
    signature_page_id: str,
    top_k_datasets: int = 20,
    top_k_chunks: int = 100,
) -> EvidenceReport:
    """
    Generate an evidence report for a Signature.

    Args:
        signature_page_id: Notion page ID of signature (with dashes)
        top_k_datasets: Number of top matching datasets to include
        top_k_chunks: Number of top chunks per dataset

    Returns:
        EvidenceReport object
    """
    logger.info(
        "[REPORTING][EVIDENCE] Generating evidence report for signature %s",
        signature_page_id,
    )

    # Use existing cross-omics summary function
    summary = cross_omics_signature_summary(
        signature_page_id=signature_page_id,
        top_k_datasets=top_k_datasets,
        top_k_chunks=top_k_chunks,
    )

    # Extract signature name
    signature_name = f"Signature {signature_page_id[:8]}"
    try:
        from amprenta_rag.query.cross_omics.helpers import extract_text_property, fetch_notion_page

        page = fetch_notion_page(signature_page_id)
        signature_name = extract_text_property(page, "Name") or signature_name
    except Exception as e:
        logger.debug("[REPORTING][EVIDENCE] Could not fetch signature name: %r", e)

    # Extract metadata
    metadata = {}
    try:
        from amprenta_rag.query.cross_omics.helpers import extract_select_values, fetch_notion_page

        page = fetch_notion_page(signature_page_id)
        props = page.get("properties", {}) or {}

        if "Modalities" in props:
            metadata["modalities"] = extract_select_values(props["Modalities"])
    except Exception as e:
        logger.debug("[REPORTING][EVIDENCE] Could not extract metadata: %r", e)

    return EvidenceReport(
        entity_type="signature",
        entity_id=signature_page_id,
        entity_name=signature_name,
        summary=summary,
        metadata=metadata,
    )


def generate_feature_evidence_report(
    feature_name: str,
    feature_type: str,
    top_k_datasets: int = 20,
    top_k_chunks: int = 100,
) -> EvidenceReport:
    """
    Generate an evidence report for a Feature.

    Args:
        feature_name: Name of the feature
        feature_type: Type of feature ("gene", "protein", "metabolite", "lipid")
        top_k_datasets: Number of top datasets to include
        top_k_chunks: Number of top chunks per dataset

    Returns:
        EvidenceReport object
    """
    logger.info(
        "[REPORTING][EVIDENCE] Generating evidence report for %s feature '%s'",
        feature_type,
        feature_name,
    )

    # Use existing cross-omics summary function
    summary = cross_omics_feature_summary(
        feature_name=feature_name,
        feature_type=feature_type,
        top_k_datasets=top_k_datasets,
        top_k_chunks=top_k_chunks,
    )

    return EvidenceReport(
        entity_type="feature",
        entity_id=feature_name,  # Features don't have Notion page IDs
        entity_name=feature_name,
        summary=summary,
        metadata={"feature_type": feature_type},
    )


def format_evidence_report(
    report: EvidenceReport,
    include_metadata: bool = True,
) -> str:
    """
    Format an evidence report as a complete markdown document.

    Args:
        report: EvidenceReport object
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
    report: EvidenceReport,
    target_page_id: Optional[str] = None,
    property_name: str = "Evidence Report",
) -> bool:
    """
    Write an evidence report to a Notion page.

    Args:
        report: EvidenceReport object
        target_page_id: Notion page ID to write to (if None, uses report.entity_id)
        property_name: Name of the property to write to

    Returns:
        True if successful, False otherwise
    """
    if target_page_id is None:
        target_page_id = report.entity_id

    logger.info(
        "[REPORTING][EVIDENCE] Writing evidence report to Notion page %s (property: %s)",
        target_page_id,
        property_name,
    )

    try:
        import requests

        # DEPRECATED: Notion imports removed
        # from amprenta_rag.clients.notion_client import notion_headers
        from amprenta_rag.config import get_config
        
        def notion_headers() -> Dict[str, str]:
            """DEPRECATED: Notion support removed. Returns empty headers dict."""
            logger.debug("[REPORTING][EVIDENCE] notion_headers() deprecated - Notion support removed")
            return {}

        cfg = get_config()

        # Format report
        formatted_report = format_evidence_report(report, include_metadata=True)

        # Truncate if too long (Notion rich text limit is ~2000 characters per block)
        if len(formatted_report) > 2000:
            formatted_report = formatted_report[:1997] + "..."

        # Update Notion page
        url = f"{cfg.notion.base_url}/pages/{target_page_id}"
        payload = {
            "properties": {
                property_name: {
                    "rich_text": [
                        {
                            "type": "text",
                            "text": {"content": formatted_report},
                        }
                    ],
                },
            },
        }

        resp = requests.patch(
            url,
            headers=notion_headers(),
            json=payload,
            timeout=30,
        )
        resp.raise_for_status()

        logger.info(
            "[REPORTING][EVIDENCE] Successfully wrote evidence report to %s",
            target_page_id,
        )
        return True

    except Exception as e:
        logger.warning(
            "[REPORTING][EVIDENCE] Error writing evidence report to Notion: %r",
            e,
        )
        return False
